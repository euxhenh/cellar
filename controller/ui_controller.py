import dash
import dash_bootstrap_components as dbc
from dash import html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot


# Manage collapsibles
def _get_toggler(collapsible_id, button_id):
    # Get function that toggles collapsible
    def _func(n, is_open):
        return not is_open
    return _func


for collapsible_id, button_id in [
    ("collapsible-dataset-bar", "change-dataset-btn"),
    ("collapsible-prep", "preprocessing-mode-btn"),
    ("dim-reduce-collapse", "dim-reduce-collapse-btn"),
    ("clustering-collapse", "clustering-collapse-btn"),
    ("annotation-collapse", "annotation-collapse-btn"),
    ("tools-collapse", "tools-collapse-btn"),
    ("session-collapse", "session-collapse-btn")
]:
    app.callback(
        Output(collapsible_id, "is_open"),
        Input(button_id, "n_clicks"),
        State(collapsible_id, "is_open"),
        prevent_initial_call=True
    )(_get_toggler(collapsible_id, button_id))


@app.callback(
    Output("main-plot-col", "width"),
    Output("side-plot-col", "width"),
    Output("main-plot-col", "className"),
    Output("side-plot-col", "className"),

    Output("dual-mode", "data"),
    Output("dual-mode-btn", "children"),

    Output("main-activate-btn", "disabled"),
    Output("side-activate-btn", "disabled"),

    Input("dual-mode-btn", "n_clicks"),

    State("dual-mode", "data"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def toggle_dual_mode(n1, is_dual, actp):
    w1, w2 = 12, 0
    d1, d2 = "block-display", "no-display"

    dm_label = [
        html.I(className="fas fa-clone me-2"),
        "Dual Mode"
    ]
    a1 = a2 = True

    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    if is_dual is None or not is_dual:
        w1, w2 = 6, 6
        d1, d2 = "block-display", "block-display"
        dm_label = [
            html.I(className="fas fa-square me-2"),
            "Single View"
        ]
        a1 = a2 = False
    elif is_dual:
        if actp == 2:
            w1, w2 = 0, 12
            d1, d2 = "no-display", "block-display"

    return w1, w2, d1, d2, not is_dual, dm_label, a1, a2


@app.callback(
    Output("main-cardheader", "className"),
    Output("side-cardheader", "className"),
    Output("active-plot", "data"),
    Output("main-analysis-col", "className"),
    Output("side-analysis-col", "className"),

    Output("main-annotation-addon", "className"),
    Output("side-annotation-addon", "className"),
    Output("main-subset-select", "className"),
    Output("side-subset-select", "className"),
    Output("shape-signal-atoggle", "data"),

    Output("ssclu-Leiden-main-checklist", "className"),
    Output("ssclu-Leiden-side-checklist", "className"),

    Output("main-annotation-table-row", "className"),
    Output("side-annotation-table-row", "className"),

    Input("main-activate-btn", "n_clicks"),
    Input("side-activate-btn", "n_clicks"),

    State("dual-mode", "data"),
    State("active-plot", "data")
)
def activate_plot(n1, n2, dual_mode, actp):
    mc1, mc2 = "active-plot-bg", "non-active-plot-bg"
    new_actp = 1
    mac1, mac2 = "block-display", "no-display"

    ctx = dash.callback_context
    if not ctx.triggered or not dual_mode:
        to_return = [dash.no_update] * 14
        to_return[2] = 1
        return to_return

    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    if button_id == "main-activate-btn" and actp == 1:
        raise PreventUpdate
    if button_id == "side-activate-btn" and actp == 2:
        raise PreventUpdate

    if button_id == "side-activate-btn":
        mc1, mc2 = mc2, mc1
        new_actp = 2
        mac1, mac2 = mac2, mac1

    return mc1, mc2, new_actp, mac1, \
        mac2, mac1, mac2, mac1, mac2, 1, \
        mac1, mac2, mac1 + " mb-2", mac2 + " mb-2"
