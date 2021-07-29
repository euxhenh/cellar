import os
import gc

import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot
from .cellar.core import read_adata, cl_add_gene_symbol
from .cellar.utils.misc import is_sparse
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification


@app.callback(
    Output("shape-signal-load", "data"),
    Output("feature-list-signal", "data"),
    Output("dataset-load", "style"),
    Output("data-loaded-plot-signal", "data"),
    Output("data-loaded-annotation-table-signal", "data"),
    MultiplexerOutput("push-notification", "data"),

    Input("load-dataset-btn", "n_clicks"),
    State("dataset-dropdown", "value"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def load_dataset(n1, dname, actp):
    """
    Loads dataset to the active plot given its path.
    """
    ctx = dash.callback_context
    if not ctx.triggered or n1 is None:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    # Initialize empty dataset
    dbroot.adatas[an] = {}
    gc.collect()

    try:
        dbroot.adatas[an]['adata'] = read_adata(dname)
    except Exception as e:
        logger.error(str(e))
        error_msg = "Error encountered while reading file. Maybe the file " +\
            "is not formatted properly, or invalid entries were found."
        logger.error(error_msg)
        return [dash.no_update] * 5 + [_prep_notification(error_msg, "danger")]

    dbroot.adatas[an]['name'] = os.path.splitext(
        os.path.basename(dname))[0]
    logger.info(f"Read {dname} info {an}.")
    if is_sparse(dbroot.adatas[an]['adata'].X):
        logger.info("Found Sparse Matrix.")

    return 1, 1, {}, 1, 1, dash.no_update


@app.callback(
    Output("no-samples", "children"),
    Output("no-features", "children"),

    Input("shape-signal-load", "data"),
    Input("shape-signal-load-prep", "data"),
    Input("shape-signal-load-prep-atac", "data"),
    Input("shape-signal-atoggle", "data"),

    State("active-plot", "data"),
    prevent_initial_call=True
)
def update_shape(s1, s2, s3, s4, actp):
    """
    When a new dataset is loaded, or the active plot is switched,
    return the shape of the data. Otherwise, return N/A.
    """
    nos, nof = 'N/A', 'N/A'

    ctx = dash.callback_context
    if not ctx.triggered or (s1 is None and s2 is None and s3 is None):
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    nos, nof = dbroot.adatas[an]['adata'].shape

    return nos, nof


@app.callback(
    Output("main-feature-list", "options"),
    Output("side-feature-list", "options"),
    MultiplexerOutput("push-notification", "data"),

    Input("feature-list-signal", "data"),
    Input("feature-list-signal-prep", "data"),
    Input("feature-list-signal-prep-atac", "data"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def update_feature_list(s1, s2, s3, actp):
    ctx = dash.callback_context
    if not ctx.triggered or (s1 is None and s2 is None):
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    if dbroot.adatas[an]['adata'].shape[1] > 100000:
        error_msg = "Too many features found. Skipping feature list."
        logger.warn(error_msg)
        to_return = [dash.no_update] * 2
        to_return[actp - 1] = []
        return to_return + [
            _prep_notification(error_msg, "warning")]

    try:
        if 'gene_symbols' not in dbroot.adatas[an]['adata'].var:
            cl_add_gene_symbol(dbroot.adatas[an]['adata'])
    except Exception as e:
        logger.error(str(e))
        error_msg = "Error occurred when adding gene symbols."
        logger.error(error_msg)
        return [dash.no_update] * 2 + [_prep_notification(error_msg, "danger")]

    features = dbroot.adatas[an][
        'adata'].var['gene_symbols'].to_numpy().astype(str)
    unique_index = \
        dbroot.adatas[an]['adata'].var_names.to_numpy().astype(str)

    to_return = [dash.no_update] * 2
    to_return[actp - 1] = [{'label': f, 'value': g}
                           for f, g in zip(features, unique_index)]

    return to_return + [dash.no_update]
