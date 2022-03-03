import dash
import os
from dash import dcc
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from flask import request

from app import app, logger, dbroot


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
    path = os.path.join('tmp', title)
    logger.info("Compressing AnnData object.")

    if os.path.isfile(path):
        os.remove(path)

    dbroot.adatas[an]['adata'].write_h5ad(path, compression=9)

    return dcc.send_file(path)


@app.callback(
    Output("export-annotations-d", "data"),
    Input("export-annotations-btn", "n_clicks"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def export_annotations(n1, actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    if an not in dbroot.adatas:
        raise PreventUpdate

    title = dbroot.adatas[an]['name'] + '_obs_cellar.csv'
    return dcc.send_data_frame(dbroot.adatas[an]['adata'].obs.to_csv, title)


# @app.callback(
#     Output("shutdown-signal", "data"),
#     Input("shutdown-btn", "n_clicks"),
#     prevent_initial_call=True
# )
# def shutdown(n1):
#     ctx = dash.callback_context
#     if not ctx.triggered:
#         raise PreventUpdate

#     func = request.environ.get('werkzeug.server.shutdown')
#     if func is None:
#         raise RuntimeError('Not running with the Werkzeug Server')

#     logger.info("Shutting down...")
#     func()
