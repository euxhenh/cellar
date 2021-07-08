import dash
import numpy as np
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, dbroot, logger


def _fill_annotation(adata, cluster_id, value):
    if isinstance(cluster_id, list):
        cluster_id = cluster_id[0]
    if isinstance(value, list):
        value = value[0]
    value = str(value)

    if cluster_id.startswith('main-cluster'):
        cluster_id = cluster_id[len('main-cluster'):]
    elif cluster_id.startswith('side-cluster'):
        cluster_id = cluster_id[len('side-cluster'):]

    cluster_id = int(cluster_id)

    if 'annotations' not in adata.obs:
        adata.obs['annotations'] = np.array(
            [""] * adata.shape[0], dtype='U200')

    annotations_copy = adata.obs['annotations'].to_numpy().copy()
    annotations_copy[adata.obs['labels'].to_numpy() == cluster_id] = value
    adata.obs['annotations'] = annotations_copy


@app.callback(
    Output("annotation-signal", "data"),

    Input("annotation-store-btn", "n_clicks"),
    State("main-annotation-select", "value"),
    State("side-annotation-select", "value"),
    State("annotation-input", "value"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def signal_annotation_change(n1, id1, id2, value, actp):
    ctx = dash.callback_context
    if not ctx.triggered or n1 is None:
        raise PreventUpdate

    cluster_id = id1 if actp == 1 else id2
    an = 'a1' if actp == 1 else 'a2'

    _fill_annotation(dbroot.adatas[an]['adata'], cluster_id, value)

    return 1


def get_update_annotation_table(prefix, an):
    def _func(s1, s2, actp):
        data = [
            {
                "cluster_id": "N/A",
                "annotation": "N/A"
            }
        ]

        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            return data
        if 'adata' not in dbroot.adatas[an]:
            return data

        # Need labels and annotations keys to be populated
        if 'labels' not in dbroot.adatas[an]['adata'].obs:
            logger.warn("No labels found in adata.")
            return data
        if 'annotations' not in dbroot.adatas[an]['adata'].obs:
            logger.warn("No annotations found in adata.")
            return data

        # Get cluster ID's and the first index for each
        # so that we can also get annotations
        unq_labels, unq_indices = np.unique(
            dbroot.adatas[an]['adata'].obs['labels'].to_numpy(),
            return_index=True)
        unq_annotations = dbroot.adatas[an][
            'adata'][unq_indices].obs['annotations']

        data = [
            {'cluster_id': str(i), 'annotation': str(j)}
            for i, j in zip(unq_labels, unq_annotations)
            if j != ""
        ]

        return data
    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-annotation-table", "data"),

        Input("annotation-signal", "data"),
        Input("data-loaded-annotation-table-signal", "data"),
        State("active-plot", "data"),
        prevent_initial_call=True
    )(get_update_annotation_table(prefix, an))
