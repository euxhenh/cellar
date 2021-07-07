import os

import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger, dbroot
from .cellar.core import read_adata, cl_add_gene_symbol
from .cellar.utils.misc import is_sparse


@app.callback(
    Output("shape-signal-upload", "data"),
    Output("feature-list-signal", "data"),
    Output("dataset-load", "style"),
    Output("data-loaded-plot-signal", "data"),

    Input("load-dataset-btn", "n_clicks"),
    State("dataset-dropdown", "value"),
    State("active-plot", "data"),
    prevent_initial_call=True
)
def load_dataset(n1, dname, actp):
    ctx = dash.callback_context
    if not ctx.triggered or n1 is None:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    dbroot.adatas[an] = {}
    dbroot.adatas[an]['adata'] = read_adata(dname)
    dbroot.adatas[an]['name'] = os.path.splitext(
        os.path.basename(dname))[0]
    logger.info(f"Read {dname} info {an}.")
    if is_sparse(dbroot.adatas[an]['adata'].X):
        logger.info("Found Sparse Matrix.")

    return 1, 1, {}, 1


@app.callback(
    Output("no-samples", "children"),
    Output("no-features", "children"),

    Input("shape-signal-upload", "data"),
    Input("shape-signal-upload-prep", "data"),
    Input("shape-signal-upload-prep-atac", "data"),
    Input("shape-signal-atoggle", "data"),

    State("active-plot", "data"),
    prevent_initial_call=True
)
def update_shape(s1, s2, s3, s4, actp):
    nos, nof = 'N/A', 'N/A'

    ctx = dash.callback_context
    if not ctx.triggered or (s1 is None and s2 is None and s3 is None):
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'

    if an in dbroot.adatas:
        nos, nof = dbroot.adatas[an]['adata'].shape

    return nos, nof


@app.callback(
    Output("main-feature-list", "options"),
    Output("side-feature-list", "options"),

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

    if an in dbroot.adatas:
        if dbroot.adatas[an]['adata'].shape[1] > 100000:
            logger.warn("Too many features found. Skipping feature list.")
            raise PreventUpdate

        if 'gene_symbols' not in dbroot.adatas[an]['adata'].var:
            cl_add_gene_symbol(dbroot.adatas[an]['adata'])

        features = dbroot.adatas[an][
            'adata'].var['gene_symbols'].to_numpy().astype(str)
        unique_index = \
            dbroot.adatas[an]['adata'].var_names.to_numpy().astype(str)

        to_return = [dash.no_update] * 2
        to_return[actp - 1] = [{'label': f, 'value': g}
                               for f, g in zip(features, unique_index)]

        return to_return

    raise PreventUpdate
