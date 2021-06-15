import dash
import numpy as np
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


def get_de_cluster_update_func(prefix, an):
    def _func(s1, s2, s3):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        clusters = []
        clusters_ann = []

        value = None

        if 'labels' in dbroot.adatas[an]['adata'].obs:
            unq_labels = np.unique(
                dbroot.adatas[an]['adata'].obs['labels'].to_numpy())
            value = prefix + "-cluster" + str(unq_labels[0])

            clusters = [{
                "label": "Cluster " + str(i),
                "value": prefix + "-cluster" + str(i)
            } for i in unq_labels]

            clusters_ann = clusters.copy()

        if 'subsets' in dbroot.adatas[an]['adata'].uns:
            subsets = list(dbroot.adatas[an]['adata'].uns['subsets'].keys())
            if len(clusters) == 0:
                value = prefix + "-subset-" + subsets[0]
            subsets = [{"label": s, "value": prefix + "-subset-" + s}
                       for s in subsets]

            clusters = clusters + subsets

        return clusters, value, clusters_ann, value, clusters, []
    return _func


# Loop over both plots
for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + '-de-cluster-select', "options"),
        Output(prefix + '-de-cluster-select', "value"),
        Output(prefix + '-annotation-select', "options"),
        Output(prefix + '-annotation-select', "value"),
        Output('ssclu-Leiden-' + prefix + '-checklist', "options"),
        Output('ssclu-Leiden-' + prefix + '-checklist', "value"),

        Input(prefix + '-cluster-list-signal', "data"),
        Input(prefix + '-subset-list-signal', "data"),
        Input(prefix + '-cluster-merge-signal', "data"),
        prevent_initial_call=True
    )(get_de_cluster_update_func(prefix, an))
