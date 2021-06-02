import dash
import os
import dash_core_components as dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot
from .cellar.core import read_adata


@app.callback(
    Output("export-session-d", "data"),
    Input("export-session-btn", "n_clicks"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def export_session(n1, actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    if an not in dbroot.adatas:
        raise PreventUpdate

    title = dbroot.adatas[an]['name'] + '_cellar.h5ad'
    path = os.path.join('data/tmp', title)
    dbroot.adatas[an]['adata'].write_h5ad(path, compression=9)

    return dcc.send_file(path)


@app.callback(
    Output("export-annotations-d", "data"),
    Input("export-annotations-btn", "n_clicks"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def export_session(n1, actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    if an not in dbroot.adatas:
        raise PreventUpdate

    title = dbroot.adatas[an]['name'] + '_obs_cellar.csv'
    return dcc.send_data_frame(dbroot.adatas[an]['adata'].obs.to_csv, title)
