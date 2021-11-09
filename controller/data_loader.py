import os
import gc

import numpy as np
import dash
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import dash_html_components as html
import dash_bootstrap_components as dbc

from app import app, logger, dbroot
from controller.cellar.utils.exceptions import UserError
from .cellar.core import read_adata, cl_add_gene_symbol
from .cellar.utils.misc import is_sparse, _check_proteins
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification


@app.callback(
    Output("shape-signal-load", "data"),
    Output("feature-list-signal", "data"),
    Output("dataset-load", "style"),
    Output("data-loaded-plot-signal", "data"),
    Output("data-loaded-annotation-table-signal", "data"),
    Output("main-data-load-clean", "data"),
    Output("side-data-load-clean", "data"),
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
        adata = read_adata(dname)
        if adata.shape[0] <= 2:
            raise UserError("Dataset contains less than 3 points.")
        dbroot.adatas[an]['adata'] = adata
    except Exception as e:
        logger.error(str(e))
        error_msg = "Error encountered while reading file. Maybe the file " +\
            "is not formatted properly, or invalid entries were found."
        logger.error(error_msg)
        return [dash.no_update] * 7 + [_prep_notification(error_msg, "danger")]

    dbroot.adatas[an]['name'] = os.path.splitext(
        os.path.basename(dname))[0]
    logger.info(f"Read {dname} info {an}.")
    if is_sparse(dbroot.adatas[an]['adata'].X):
        logger.info("Found Sparse Matrix.")

    mclean = 1 if an == 'a1' else dash.no_update
    sclean = 1 if an == 'a2' else dash.no_update

    return 1, 1, {}, 1, 1, mclean, sclean, dash.no_update


def _prettify_time(s):
    s = int(np.round(s))
    hrs, s = divmod(s, 3600)
    m, s = divmod(s, 60)
    t = f"{hrs}h " if hrs > 0 else ""
    t += f"{m}m " if m > 0 else ""
    t += f"{s}s"
    return t


def _get_singler_warning(nos, nof):
    nos, nof = nos / 1000, nof / 1000
    singler = 120 if nos < 5 else 0.18 * nos**2 + 17 * nos + 28
    msg = "Running SingleR on this dataset can " + \
        "take between " + _prettify_time(singler) + " to " + \
        _prettify_time(2 * singler) + ". Do you " + \
        "want to continue?"
    return msg


def _get_run_times(nos, nof):
    nos /= 1000
    nof /= 1000

    table_header = [
        html.Thead(html.Tr([html.Th("Method"), html.Th("Est. Runtime")]))
    ]

    # Assume quadratic model on data
    umap = 5 if nos < 5 else 0.002 * nos**2 + 0.21 * nos + 8
    if nof > 10:
        umap *= (nof / 5)
    umaprow = html.Tr([html.Td("PCA (40 PC) + UMAP"), html.Td(
        _prettify_time(umap))])

    # Assume quadratic model on data
    tsne = 15 if nos < 5 else 0.02 * nos**2 + 0.13 * nos + 20
    if nof > 10:
        tsne *= (nof / 5)
    tsnerow = html.Tr([html.Td("PCA (40 PC) + t-SNE"), html.Td(
        _prettify_time(tsne))])

    pca = 1 if nos < 5 else 0.00031 * nos ** 2 + 0.0465 * nos + 0.97
    if nof > 10:
        pca *= (nof / 5)
    pcarow = html.Tr([html.Td("PCA (40 PC) + PCA"), html.Td(
        _prettify_time(pca))])

    # Assume linear model on data (since using approx nn)
    leiden = 2 if nos < 5 else 0.7 * nos
    leidenrow = html.Tr([html.Td("Leiden Clustering (auto + max iter)"),
                         html.Td(_prettify_time(leiden))])

    # Assume linear model on data
    de = 1 if nos < 5 else 1.4 * nos
    if nof > 5:
        de *= (nof / 5)
    derow = html.Tr([html.Td("DE Analysis (1 vs rest)"), html.Td(
        _prettify_time(de))])

    # Need to run 2 pca's and 2 umaps so this should be a good approx
    scanpyrow = html.Tr([html.Td("Scanpy Ingest (Label Transfer)"), html.Td(
        _prettify_time(2 * umap) + " - " + _prettify_time(4 * umap))])

    # Assume quadratic model on data
    singler = 120 if nos < 5 else 0.18 * nos**2 + 17 * nos + 28
    singlerrow = html.Tr([html.Td("SingleR (Label Transfer)"), html.Td(
        _prettify_time(singler) + " - " + _prettify_time(2 * singler))])

    table_body = [html.Tbody([
        umaprow, tsnerow, pcarow, leidenrow, derow, scanpyrow, singlerrow])]
    table = dbc.Table(table_header + table_body, bordered=False)
    return table


@app.callback(
    Output("no-samples", "children"),
    Output("no-features", "children"),
    Output("run-time-popover-body", "children"),

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

    return nos, nof, [
        html.Span("Rough estimates based on similar datasets"),
        html.Br(), html.Br(),
        _get_run_times(nos, nof)
    ]


@app.callback(
    Output("main-feature-list", "options"),
    Output("side-feature-list", "options"),
    Output("main-feature-list-spatial", "options"),
    Output("side-feature-list-spatial", "options"),
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
    adata = dbroot.adatas[an]['adata']

    if adata.shape[1] > 40000:
        error_msg = "Too many features found. Skipping feature list."
        logger.warn(error_msg)
        to_return = [dash.no_update] * 4
        to_return[actp - 1] = to_return[actp + 1] = []
        return to_return + [
            _prep_notification(error_msg, "warning")]

    try:
        if 'gene_symbols' not in adata.var:
            cl_add_gene_symbol(adata)
    except Exception as e:
        logger.error(str(e))
        error_msg = "Error occurred when adding gene symbols."
        logger.error(error_msg)
        return [dash.no_update] * 4 + [_prep_notification(error_msg, "danger")]

    features = adata.var['gene_symbols'].to_numpy().astype(str)
    unique_index = adata.var_names.to_numpy().astype(str)
    # Check for CITE-seq data
    msg = None
    proteins = _check_proteins(adata)
    if proteins is not None:
        features = np.concatenate([features, proteins])
        unique_index = np.concatenate([unique_index, proteins])
        msg = "Found proteins in file. These were added as features " +\
            "prefixed by Protein_."

    to_return = [dash.no_update] * 4
    to_return[actp - 1] = to_return[actp + 1] = [
        {'label': f, 'value': g}
        for f, g in zip(features, unique_index)]

    return to_return + [
        dash.no_update if msg is None else _prep_notification(msg, "info")]
