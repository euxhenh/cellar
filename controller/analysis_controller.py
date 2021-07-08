import dash
import numpy as np
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .cellar.core import ttest, enrich, get_heatmap, get_violin_plot


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
        
        clusters_DE2 = [{'label':'rest','value':'rest'}] + clusters_DE
        value_DE2 =  'rest'
        
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
    def _func(n1, cluster_id, cluster_id2):
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
            if (isinstance(cluster_id, int)):
                logger.error("No labels found in adata. Cannot find DE "
                         f"genes for cluster {cluster_id}.")
                raise PreventUpdate
            if (isinstance(cluster_id2, int)):
                logger.error("No labels found in adata. Cannot find DE "
                         f"genes for cluster {cluster_id2}.")
                raise PreventUpdate
                
        title = f"DE genes for cluster {cluster_id} vs {cluster_id2}"

        logger.info("Running t-Test.")

        # Run tests
        try:
            test = ttest(dbroot.adatas[an]['adata'], cluster_id, cluster_id2)
        except Exception as e:
            logger.error(str(e))
            logger.error("An error occurred while running t-Test.")
            raise PreventUpdate

        test['id'] = test['gene'].copy()

        return [{
            "name": i.upper(),
            "id": i
        } for i in test.columns if i != 'id'],\
            test.to_dict('records'), 'csv', title

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-de-table", "columns"),
        Output(prefix + "-de-table", "data"),
        Output(prefix + "-de-table", "export_format"),
        Output(prefix + "-de-analysis-title", "children"),
        Input(prefix + "-find-de-genes-btn", "n_clicks"),
        State(prefix + "-de-cluster-select", "value"),
        State(prefix + "-de-cluster-select2", "value"),
        prevent_initial_call=True
    )(get_update_de_table_func(prefix, an))


# Enrichment
def get_enrichment_table_func(prefix, an):
    def _func(n1, gene_set, de_genes_no, de_data):
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

        if de_data is None or len(de_data) == 0:
            logger.info("Empty DE list encountered.")
            raise PreventUpdate

        # Make sure de_genes_no is an int
        try:
            de_genes_no = int(de_genes_no)
        except Exception as e:
            logger.error(str(e))
            logger.error("DE genes no. was not an int.")
            raise PreventUpdate

        if de_genes_no < 1:
            logger.error("No. genes to use is less than 1.")
            raise PreventUpdate

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
            logger.error("Error whle running enrichment analysis for gene "
                         f"set {gene_set}.")
            raise PreventUpdate

        return [{
            "name": i.upper(),
            "id": i
        } for i in enr.columns], enr.to_dict('records'), 'csv', title

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
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
    )(get_enrichment_table_func(prefix, an))


# Analysis plots
def get_plot_analysis_func(prefix, an):
    def _func(n1, n2, feature_list):
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

        if 'labels' not in dbroot.adatas[an]['adata'].obs:
            logger.info("No labels found. Cannot visualize features.")
            raise PreventUpdate

        if feature_list is None or len(feature_list) == 0:
            logger.info("No features selected.")
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if button_id == prefix + '-heatmap':
            try:
                fig = get_heatmap(dbroot.adatas[an]['adata'], feature_list)
            except Exception as e:
                logger.error(str(e))
                logger.error("Error when running heatmap.")
                raise PreventUpdate
        elif button_id == prefix + '-violin-plot':
            if len(feature_list) > 1:
                logger.warn("Cannot run violin plot with "
                            "more than one feature.")
                raise PreventUpdate
            try:
                fig = get_violin_plot(
                    dbroot.adatas[an]['adata'],
                    feature_list[0],
                    plot1=(an == 'a1'))
            except Exception as e:
                logger.error(str(e))
                logger.error("Error when running violin plot.")
                raise PreventUpdate

        return fig

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-analysis-plot", "figure"),

        Input(prefix + "-heatmap", "n_clicks"),
        Input(prefix + "-violin-plot", "n_clicks"),
        State(prefix + "-feature-list", "value"),
        prevent_initial_call=True
    )(get_plot_analysis_func(prefix, an))
