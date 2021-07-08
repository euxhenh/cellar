import numpy as np

import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot
from controller.cellar.utils.exceptions import InternalError


@app.callback(
    Output("main-subset-list-signal", "data"),
    Output("side-subset-list-signal", "data"),
    Output("main-subset-select", "options"),
    Output("side-subset-select", "options"),

    Input("subset-name-store-btn", "n_clicks"),
    State("subset-name-input", "value"),
    State("main-plot", "selectedData"),
    State("side-plot", "selectedData"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def store_subset(n1, subset_name, main_selected, side_selected, actp):
    if subset_name.startswith('Cluster'):
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    prefix = "main" if actp == 1 else "side"

    uns = dbroot.adatas[an]['adata'].uns
    selected = main_selected if actp == 1 else side_selected

    if selected is None:
        raise PreventUpdate
    if len(selected['points']) == 0:
        raise PreventUpdate

    if 'subsets' not in uns:
        uns['subsets'] = {}

    points = np.array([p['customdata'][0] for p in selected['points']])
    logger.info(f"Storing {len(points)} points into {subset_name}.")

    uns['subsets'][subset_name] = points

    subsets = sorted(list(uns['subsets'].keys()))
    subsets = [{
        "label": subset, "value": prefix + "-" + subset
    } for subset in subsets]

    if actp == 1:
        return 1, dash.no_update, subsets, dash.no_update
    return dash.no_update, 1, dash.no_update, subsets


@app.callback(
    Output("merge-plot-signal", "data"),
    Output("main-cluster-merge-signal", "data"),
    Output("side-cluster-merge-signal", "data"),

    Input("merge-subset-btn", "n_clicks"),
    State("main-subset-select", "value"),
    State("side-subset-select", "value"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def merge_subset(n1, main_select, side_select, actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    prefix = "main" if actp == 1 else "side"
    select = main_select if actp == 1 else side_select
    select = select[len(prefix + "-"):]

    if 'subsets' not in dbroot.adatas[an]['adata'].uns:
        raise InternalError("No subsets dict found.")
    if 'labels' not in dbroot.adatas[an]['adata'].obs:
        raise PreventUpdate

    max_cluster_id = np.max(
        dbroot.adatas[an]['adata'].obs['labels'].to_numpy())

    points = dbroot.adatas[an]['adata'].uns['subsets'][select]
    logger.info(f"Merging {len(points)} points from subset {select}.")

    labels = dbroot.adatas[an]['adata'].obs['labels'].to_numpy()
    labels[points] = max_cluster_id + 1
    dbroot.adatas[an]['adata'].obs['labels'] = labels

    if 'annotations' in dbroot.adatas[an]['adata'].obs:
        annotations = dbroot.adatas[an]['adata'].obs['annotations'].to_numpy()
        annotations[points] = ""
        dbroot.adatas[an]['adata'].obs['annotations'] = annotations

    if actp == 1:
        return 1, 1, dash.no_update
    return 1, dash.no_update, 1
