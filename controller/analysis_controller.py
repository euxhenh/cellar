import dash
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .cellar.core import ttest, enrich, get_heatmap, get_violin_plot


# DE genes
def get_update_de_table_func(prefix):
    actp = 1 if prefix == 'main' else 2
    an = 'a1' if prefix == 'main' else 'a2'

    def _func(n1, cluster_id):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        if an not in dbroot.adatas:
            raise PreventUpdate

        if 'labels' not in dbroot.adatas[an]['adata'].obs:
            raise PreventUpdate

        cluster_id = int(cluster_id[len(prefix + "-cluster"):])

        title = f"PLOT {actp}: DE genes for cluster {cluster_id} vs rest"

        test = ttest(dbroot.adatas[an]['adata'], cluster_id)
        test['id'] = test['gene'].copy()

        return [{"name": i.upper(), "id": i}
                for i in test.columns if i != 'id'],\
            test.to_dict('records'), 'csv', title

    return _func


for prefix in ["main", "side"]:
    app.callback(
        Output(prefix + "-de-table", "columns"),
        Output(prefix + "-de-table", "data"),
        Output(prefix + "-de-table", "export_format"),
        Output(prefix + "-de-analysis-title", "children"),
        Input(prefix + "-find-de-genes-btn", "n_clicks"),
        State(prefix + "-de-cluster-select", "value"),
        prevent_initial_call=True
    )(get_update_de_table_func(prefix))


# Enrichment
def get_enrichment_table_func(prefix):
    actp = 1 if prefix == 'main' else 2
    an = 'a1' if prefix == 'main' else 'a2'

    def _func(n1, gene_set, de_genes_no, de_data):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        if de_data is None or len(de_data) == 0:
            raise PreventUpdate
        de_gene_list = [r['gene'] for r in de_data][:de_genes_no]

        gene_set = gene_set[len('gene-set-'):]
        logger.info(f"Running enrichment analysis for {gene_set} gene set " +
                    f"with {de_genes_no} top genes.")

        title = f"PLOT {actp}: {gene_set}"

        enr = enrich(dbroot.adatas[an]['adata'], gene_set, de_gene_list)

        return [{"name": i.upper(), "id": i} for i in enr.columns],\
            enr.to_dict('records'), 'csv', title

    return _func


for prefix in ["main", "side"]:
    app.callback(
        Output(prefix + "-enrich-table", "columns"),
        Output(prefix + "-enrich-table", "data"),
        Output(prefix + "-enrich-table", "export_format"),
        Output(prefix + "-enrich-title", "children"),

        Input(prefix + "-run-enrich-btn", "n_clicks"),
        State(prefix + "-gene-set-dropdown", "value"),
        State(prefix + "-de-genes-enrich-no", "value"),
        State(prefix + "-de-table", "data"),
        prevent_initial_call=True
    )(get_enrichment_table_func(prefix))


# Analysis plots
def get_plot_analysis_func(prefix):
    actp = 1 if prefix == 'main' else 2
    an = 'a1' if prefix == 'main' else 'a2'

    def _func(n1, n2, feature_list):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        if 'labels' not in dbroot.adatas[an]['adata'].obs:
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if button_id == prefix + '-heatmap':
            fig = get_heatmap(dbroot.adatas[an]['adata'], feature_list)
        elif button_id == prefix + '-violin-plot':
            if len(feature_list) > 1:
                raise PreventUpdate
            fig = get_violin_plot(
                dbroot.adatas[an]['adata'], feature_list[0], plot1=(actp == 1))

        return fig

    return _func


for prefix in ["main", "side"]:
    app.callback(
        Output(prefix + "-analysis-plot", "figure"),

        Input(prefix + "-heatmap", "n_clicks"),
        Input(prefix + "-violin-plot", "n_clicks"),
        State(prefix + "-feature-list", "value"),
        prevent_initial_call=True
    )(get_plot_analysis_func(prefix))
