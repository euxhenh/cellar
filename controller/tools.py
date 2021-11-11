import numpy as np

import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot
from controller.cellar.utils.exceptions import InternalError
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification
from .operations import intg_filter
from .methods import intg_list


@app.callback(
    Output("main-subset-list-signal", "data"),
    Output("side-subset-list-signal", "data"),
    Output("main-subset-select", "options"),
    Output("side-subset-select", "options"),
    MultiplexerOutput("push-notification", "data"),

    Input("subset-name-store-btn", "n_clicks"),
    State("subset-name-input", "value"),
    State("main-plot", "selectedData"),
    State("side-plot", "selectedData"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def store_subset(n1, subset_name, main_selected, side_selected, actp):
    if subset_name.startswith('Cluster'):
        return [dash.no_update] * 4 + [_prep_notification(
            "Cannot use a name that starts with 'Cluster'.", "warning")]
    if len(subset_name) == 0:
        return [dash.no_update] * 4 + [_prep_notification(
            "Empty subset name.", "warning")]

    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

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

    uns['subsets'][subset_name] = points

    subsets = sorted(list(uns['subsets'].keys()))
    subsets = [{
        "label": subset, "value": prefix + "-" + subset
    } for subset in subsets]

    msg = f"Stored {len(points)} points into {subset_name}."
    logger.info(msg)

    if actp == 1:
        return 1, dash.no_update, subsets, dash.no_update, _prep_notification(
            msg, "info")
    return dash.no_update, 1, dash.no_update, subsets, _prep_notification(
        msg, "info")


@app.callback(
    Output("merge-plot-signal", "data"),
    Output("main-cluster-merge-signal", "data"),
    Output("side-cluster-merge-signal", "data"),
    MultiplexerOutput("push-notification", "data"),

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
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    prefix = "main" if actp == 1 else "side"
    select = main_select if actp == 1 else side_select

    if select is None:
        raise PreventUpdate

    select = select[len(prefix + "-"):]
    if len(select) == 0:
        raise PreventUpdate

    if 'subsets' not in dbroot.adatas[an]['adata'].uns:
        return [dash.no_update] * 3 + [_prep_notification(
            "No subsets dictionary found.", "warning")]
    if 'labels' not in dbroot.adatas[an]['adata'].obs:
        return [dash.no_update] * 3 + [_prep_notification(
            "Cannot merge without labels.", "warning")]

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
        return 1, 1, dash.no_update, _prep_notification(
            f"Merged subset {select}.", "info")
    return 1, dash.no_update, 1, _prep_notification(
        f"Merged subset {select}.", "info")


def get_integrate_func(an):
    def _func(n1, intg_method, actp, *settings):
        ctx = dash.callback_context
        if (actp == 1 and an != 'a1') or (actp == 2 and an != 'a2'):
            raise PreventUpdate
        if not ctx.triggered:
            raise PreventUpdate
        ref_an = 'a2' if an == 'a1' else 'a1'
        if an not in dbroot.adatas or ref_an not in dbroot.adatas:
            return [], dash.no_update, _prep_notification(
                "Must load a reference dataset.", "warn")
        if 'adata' not in dbroot.adatas[an] or 'adata' not in dbroot.adatas[ref_an]:
            return [], dash.no_update, _prep_notification(
                "Must load a reference dataset.", "warn")
        try:
            intg_filter(dbroot.adatas[an]['adata'], intg_method, settings,
                        extras={'ref': dbroot.adatas[ref_an]['adata']})
        except:
            raise InternalError("An error occurred during integration.")
        return [], 1, _prep_notification("Finished integrating.")
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-buf-integrate", "children"),
        MultiplexerOutput(prefix + "-other-features-change", "data"),
        MultiplexerOutput("push-notification", "data"),
        Input("run-integrate-btn", "n_clicks"),
        State("intg-methods-select", "value"),
        State("active-plot", "data"),
        [State(m['value'] + '-settings', 'children') for m in intg_list],
        prevent_initial_call=True
    )(get_integrate_func(an))
