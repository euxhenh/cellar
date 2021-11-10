import dash
import os
import pandas as pd
import numpy as np
from zope.interface.declarations import ProvidesClass
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .cellar.core import ttest, enrich, get_heatmap, get_violin_plot
from .cellar.core._tools import _collect_x_from_vars, _collect_x_from_other
from .cellar.utils.exceptions import InternalError, UserError
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification
from layout.misc import empty_analysis_figure


def get_de_cluster_update_func(prefix, an):
    def _func(s1, s2, s3):
        """
        Populates the cluster/subset list used to search for DE genes
        and the list used for annotations. These lists differ if any
        used-defined subsets are found in adata. Only the DE list
        will show these subsets. A third list equivalent to clusters_DE
        is returned for the Freeze functionality of Constrained Leiden.

        First we determine all unique cluster ID's, and then search
        for any user-defined subsets under adata.uns['subsets'].

        The Select menu value for the clusters is of the form:
            prefix  + '-cluster' + ID
        where prefix is in ['main', 'side'], while the value for the
        subsets is of the form
            prefix + '-subset' + s
        where s is the subset name.

        Called upon 1) cluster change, 2) subset definition, 3) subset merge.
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        clusters_DE = []  # The DE cluster list (includes subsets)
        clusters_ann = []  # The annotation cluster list

        value_DE = None  # The default value for the DE list
        value_ann = None  # The default value for the annotations list

        if 'labels' in dbroot.adatas[an]['adata'].obs:
            unq_labels = np.unique(
                dbroot.adatas[an]['adata'].obs['labels'].to_numpy())
            # Set the default value to the first cluster ID
            value_DE = value_ann = prefix + "-cluster" + str(unq_labels[0])

            # Populate the list of cluster values
            clusters_DE = [{
                "label": "Cluster " + str(i),
                "value": prefix + "-cluster" + str(i)
            } for i in unq_labels]

            clusters_ann = clusters_DE.copy()

        # Add subsets to the DE cluster list only
        if 'subsets' in dbroot.adatas[an]['adata'].uns:
            # Get the subset list
            subsets = list(dbroot.adatas[an]['adata'].uns['subsets'].keys())

            # In case no labels are found, but there are subsets
            if len(clusters_DE) == 0:
                value_DE = prefix + "-subset-" + subsets[0]

            subsets = [{
                "label": s,
                "value": prefix + "-subset-" + s
            } for s in subsets]

            # Append subsets to the cluster list
            clusters_DE = clusters_DE + subsets

        clusters_DE2 = [{'label': 'rest', 'value': 'rest'}] + clusters_DE
        value_DE2 = 'rest'

        return clusters_DE, value_DE, clusters_DE2, value_DE2, \
            clusters_ann, value_ann, clusters_DE, []

    return _func


# Loop over both plots
for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + '-de-cluster-select', "options"),
        Output(prefix + '-de-cluster-select', "value"),
        Output(prefix + '-de-cluster-select2', "options"),
        Output(prefix + '-de-cluster-select2', "value"),
        Output(prefix + '-annotation-select', "options"),
        Output(prefix + '-annotation-select', "value"),
        Output('ssclu-Leiden-' + prefix + '-checklist', "options"),
        Output('ssclu-Leiden-' + prefix + '-checklist', "value"),

        Input(prefix + '-cluster-list-signal', "data"),
        Input(prefix + '-subset-list-signal', "data"),
        Input(prefix + '-cluster-merge-signal', "data"),
        prevent_initial_call=True
    )(get_de_cluster_update_func(prefix, an))


# DE genes
def get_update_de_table_func(prefix, an):
    def _func(n1, clean, cluster_id, cluster_id2, alpha, fc_thresh, is_logged):
        """
        Given a cluster id or subset name, find the DE genes for that cluster.

        First parse the id or the subset name, then run t-test.
        TODO: Add the other tests.

        Returns a list of columnts, a dictionary, and enables csv download.
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate
        if cluster_id is None:
            raise PreventUpdate

        if ctx.triggered[0]["prop_id"].split(".")[0] == \
                prefix + "-data-load-clean":
            return None, None, None, None, dash.no_update

        # If cluster_id starts with prefix + '-subset', than it is a subset,
        # otherwise it is considered a cluster ID and is converted to int.
        if cluster_id.startswith(prefix + '-cluster'):
            try:
                cluster_id = int(cluster_id[len(prefix + "-cluster"):])
            except Exception as e:
                logger.error(str(e))
                logger.error(f"Invalid cluster ID found {cluster_id}.")
                raise PreventUpdate
        else:
            cluster_id = str(cluster_id[len(prefix + "-subset-"):])

        # same thing for cluster_id2
        if cluster_id2.startswith(prefix + '-cluster'):
            try:
                cluster_id2 = int(cluster_id2[len(prefix + "-cluster"):])
            except Exception as e:
                logger.error(str(e))
                logger.error(f"Invalid cluster2 ID found {cluster_id2}.")
                raise PreventUpdate
        elif cluster_id2 == 'rest':
            pass
        else:
            cluster_id2 = str(cluster_id2[len(prefix + "-subset-"):])

        if 'labels' not in dbroot.adatas[an]['adata'].obs:
            if isinstance(cluster_id, int) or isinstance(cluster_id2, int):
                logger.error("No labels found. Cannot find DE genes.")
                raise PreventUpdate

        title = f"DE genes for cluster {cluster_id} vs {cluster_id2}"

        logger.info("Running t-Test.")

        if fc_thresh <= 0:
            fc_thresh = None

        # Run tests
        try:
            test = ttest(
                dbroot.adatas[an]['adata'],
                cluster_id, cluster_id2, alpha, fc_thresh, is_logged)
        except UserError as ue:
            logger.error(str(ue))
            return [dash.no_update] * 4 + [_prep_notification(
                str(ue), "warning")]
        except Exception as e:
            logger.error(str(e))
            error_msg = "An error occurred while running t-Test."
            logger.error(error_msg)
            return [dash.no_update] * 4 + [_prep_notification(
                error_msg, "danger")]

        test['id'] = test['gene'].copy()

        return [{
            "name": i.upper(),
            "id": i
        } for i in test.columns if i != 'id'],\
            test.to_dict('records'), 'csv', title, dash.no_update

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-de-table", "columns"),
        Output(prefix + "-de-table", "data"),
        Output(prefix + "-de-table", "export_format"),
        Output(prefix + "-de-analysis-title", "children"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-find-de-genes-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),

        State(prefix + "-de-cluster-select", "value"),
        State(prefix + "-de-cluster-select2", "value"),
        State(prefix + "-de-analysis-alpha", "value"),
        State(prefix + "-de-analysis-fc", "value"),
        State(prefix + "-de-analysis-logged", "value"),
        prevent_initial_call=True
    )(get_update_de_table_func(prefix, an))


# Enrichment
def get_enrichment_table_func(prefix, an):
    def _func(n1, clean, gene_set, de_genes_no, de_data):
        """
        Runs enrichment analysis for gene_set given a list of DE genes.
        Since we obtain the full list of DE genes, an additional
        integer 'de_genes_no' slices the top n genes for the analysis.
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        if ctx.triggered[0]["prop_id"].split(".")[0] == \
                prefix + "-data-load-clean":
            return None, None, None, None, dash.no_update

        if de_data is None or len(de_data) == 0:
            error_msg = "Empty DE list encountered."
            logger.info(error_msg)
            return [dash.no_update] * 4 + [_prep_notification(
                error_msg, "warning")]

        # Make sure de_genes_no is an int
        try:
            de_genes_no = int(de_genes_no)
        except Exception as e:
            logger.error(str(e))
            error_msg = "Number of DE genes to use was not an integer."
            logger.error(error_msg)
            return [dash.no_update] * 4 + [_prep_notification(
                error_msg, "warning")]

        if de_genes_no < 1:
            error_msg = "No. genes to use is less than 1."
            logger.error(error_msg)
            return [dash.no_update] * 4 + [_prep_notification(
                error_msg, "warning")]

        # Gets the gene symbol for every gene and takes top de_genes_no
        de_gene_list = [r['gene'] for r in de_data][:de_genes_no]

        # Parse the name of the gene set
        gene_set = gene_set[len('gene-set-'):]
        logger.info(f"Running enrichment analysis for {gene_set} gene set " +
                    f"with {de_genes_no} top genes.")

        title = f"{gene_set}"

        # Run gseapy's enrichment
        try:
            enr = enrich(dbroot.adatas[an]['adata'], gene_set, de_gene_list)
        except Exception as e:
            logger.error(str(e))
            error_msg = "Error while running enrichment analysis for gene " + \
                f"set {gene_set}."
            logger.error(error_msg)
            return [dash.no_update] * 4 + [_prep_notification(
                error_msg, "danger")]

        return [{
            "name": i.upper(), "id": i
        } for i in enr.columns], enr.to_dict('records'), 'csv', title,\
            dash.no_update

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-enrich-table", "columns"),
        Output(prefix + "-enrich-table", "data"),
        Output(prefix + "-enrich-table", "export_format"),
        Output(prefix + "-enrich-title", "children"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-run-enrich-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),

        State(prefix + "-gene-set-dropdown", "value"),
        State(prefix + "-de-genes-enrich-no", "value"),
        State(prefix + "-de-table", "data"),
        prevent_initial_call=True
    )(get_enrichment_table_func(prefix, an))


# Paste active DE genes
def get_paste_de_genes_func(prefix, an):
    def _func(n1, cells, clean, de_data, cur_page, page_size, options):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        vals = {opt['label']: i for i, opt in enumerate(options)}

        if cur_page is None:
            cur_page = 0

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return [], [], dash.no_update
        elif button_id == prefix + "-paste-de-genes":
            if de_data is None or len(de_data) == 0:
                return dash.no_update, dash.no_update, _prep_notification(
                    "DE table is empty", "warning"
                )
            selected = [
                options[vals[row['gene']]]['value'] for row in
                de_data[cur_page * page_size: (cur_page + 1) * page_size]
            ]
            return selected, selected, dash.no_update
        elif button_id == prefix + "-de-table":
            labels = [
                de_data[cur_page * page_size + cell['row']]['gene']
                for cell in cells
            ]
            labels = np.unique(labels)
            selected = [
                options[vals[label]]['value']
                for label in labels
            ]
            return selected, selected, dash.no_update

        return [], [], dash.no_update
    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-feature-list", "value"),
        Output(prefix + "-feature-list-spatial", "value"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-paste-de-genes", "n_clicks"),
        Input(prefix + "-de-table", "selected_cells"),
        Input(prefix + "-data-load-clean", "data"),

        State(prefix + "-de-table", "data"),
        State(prefix + '-de-table', "page_current"),
        State(prefix + "-de-table", "page_size"),
        State(prefix + "-feature-list", "options"),
        prevent_initial_call=True
    )(get_paste_de_genes_func(prefix, an))


# Analysis plots
def get_plot_analysis_func(prefix, an):
    def _func(n1, n2, clean, feature_list, feature_range, auto_scale):
        """
        Given a single or list of features, return a heatmap or violin plot.
        The violin plot can only be run when a single feature is selected.
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate
        adata = dbroot.adatas[an]['adata']

        if ctx.triggered[0][
                "prop_id"].split(".")[0].endswith("-data-load-clean"):
            return empty_analysis_figure, dash.no_update

        if 'labels' not in adata.obs:
            error_msg = "No labels found; cannot visualize features. Run " +\
                "clustering first."
            logger.info(error_msg)
            return dash.no_update, _prep_notification(error_msg, "warning")

        if feature_list is None or len(feature_list) == 0:
            error_msg = "No features selected."
            logger.info(error_msg)
            return dash.no_update, _prep_notification(error_msg, "warning")

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if button_id == prefix + '-heatmap':
            try:
                if prefix.endswith('other'):
                    fig = get_heatmap(adata, other_list=feature_list)
                else:
                    fig = get_heatmap(adata, feature_list)
            except UserError as ue:
                logger.error(str(ue))
                return dash.no_update, _prep_notification(str(ue), "warning")
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred while running heatmap."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        elif button_id == prefix + '-violin-plot':
            try:
                if prefix.endswith('other'):
                    fig = get_violin_plot(
                        dbroot.adatas[an]['adata'], other_values=feature_list,
                        feature_range=feature_range,
                        palette=dbroot.palettes[prefix[:4]],
                        auto_scale=auto_scale)
                else:
                    fig = get_violin_plot(
                        dbroot.adatas[an]['adata'], feature_list, feature_range,
                        palette=dbroot.palettes[prefix],
                        auto_scale=auto_scale)
            except UserError as ue:
                logger.error(str(ue))
                return dash.no_update, _prep_notification(str(ue), "warning")
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred while running violin plot."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")

        return fig, dash.no_update

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-analysis-plot", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-heatmap", "n_clicks"),
        Input(prefix + "-violin-plot", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        State(prefix + "-feature-list", "value"),
        State(prefix + "-feature-rangeslider", "value"),
        State(prefix + "-auto-scale-expression", "checked"),
        prevent_initial_call=True
    )(get_plot_analysis_func(prefix, an))


for prefix, an in zip(['main-other', 'side-other'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-analysis-plot", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-heatmap", "n_clicks"),
        Input(prefix + "-violin-plot", "n_clicks"),
        Input(prefix[:4] + "-data-load-clean", "data"),
        State(prefix + "-feature-list", "value"),
        State(prefix + "-feature-rangeslider", "value"),
        State(prefix + "-auto-scale-expression", "checked"),
        prevent_initial_call=True
    )(get_plot_analysis_func(prefix, an))


def get_feature_range_disable_func(prefix):
    def _func(auto_scale):
        if auto_scale:
            return True
        return False
    return _func


for prefix in ['main', 'side', 'main-other', 'side-other']:
    app.callback(
        Output(prefix + "-feature-rangeslider", "disabled"),
        Input(prefix + "-auto-scale-expression", "checked")
    )(get_feature_range_disable_func(prefix))


def get_feature_range_func(an, prefix):
    def _func(feature_list):
        """
        Updates the range slider to the range of the first feature selected.
        """
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate
        adata = dbroot.adatas[an]['adata']

        # More than 1 or 0 features selected
        if len(feature_list) > 1 or len(feature_list) < 1:
            return 0, 0, 0, {}, [0, 0]

        feature = feature_list[0]  # Only apply range slider to first feature
        if prefix.endswith('other'):
            feature_vec = _collect_x_from_other([feature], adata).flatten()
        else:
            feature_vec = _collect_x_from_vars([feature], adata).flatten()

        vmin, vmax = float(feature_vec.min()), float(feature_vec.max())
        rang = vmax - vmin
        eps = 1e-3
        marks_step = rang / 5  # 6 marks in total
        step = rang / 100
        disp_vmin, disp_vmax = vmin - eps, vmax + eps

        vmins, vmaxs = "{:.2e}".format(disp_vmin), "{:.2e}".format(disp_vmax)
        marks = {vmin: vmins}
        for i in range(1, 5):
            mark_val = vmin + i * marks_step
            marks[mark_val] = "{:.2e}".format(mark_val)
        marks[vmax] = vmaxs

        value = [disp_vmin, disp_vmax]

        return disp_vmin, disp_vmax, step, marks, value
    return _func


for prefix, an in zip(
        ['main', 'side', 'main-other', 'side-other'],
        ['a1', 'a2', 'a1', 'a2']):
    app.callback(
        Output(prefix + "-feature-rangeslider", "min"),
        Output(prefix + "-feature-rangeslider", "max"),
        Output(prefix + "-feature-rangeslider", "step"),
        Output(prefix + "-feature-rangeslider", "marks"),
        Output(prefix + "-feature-rangeslider", "value"),

        Input(prefix + "-feature-list", "value"),
        prevent_initial_call=True
    )(get_feature_range_func(an, prefix))


def get_asct_table_func():
    def _func(val):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        df = pd.read_csv(
            os.path.join('data/ASCT-B_tables', val),
            skiprows=10,
            encoding="ISO-8859-1")

        # df.drop([
        #     'AS/1/LABEL', 'cell type general label', 'AS/3/LABEL',
        #     'CT/1/LABEL'
        # ], axis=1, inplace=True, errors='ignore')
        # for i in range(15):
        #     df.drop([
        #         f"BG/{i}", f"BG/{i}/LABEL", f"BG/{i}/ID"
        #     ], axis=1, inplace=True, errors='ignore')
        #     df.drop([
        #         f"BP/{i}", f"BP/{i}/LABEL", f"BP/{i}/ID"
        #     ], axis=1, inplace=True, errors='ignore')
        # df.drop([
        #     'FTU', 'REF/1', 'REF/1/DOI', 'REF/1/NOTES',
        #     'REF/2', 'REF/2/DOI', 'REF/2/NOTES'
        # ], axis=1, inplace=True, errors='ignore')

        columns = [{"name": i, "id": i} for i in df.columns]
        data = df.to_dict('records')

        return columns, data

    return _func


for prefix in ['main', 'side']:
    app.callback(
        Output(prefix + "-asct-table", "columns"),
        Output(prefix + "-asct-table", "data"),

        Input(prefix + "-asct-select", "value"),
        prevent_initial_call=True
    )(get_asct_table_func())
