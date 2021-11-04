from enum import Enum
import re

import dash
import dash_core_components as dcc
import plotly
from scipy.cluster.hierarchy import single
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from layout.misc import empty_figure
from .cellar.core import (clear_x_emb_dependends, get_clu_figure,
                          get_dim_figure, get_expression_figure,
                          get_reset_figure)
from .cellar.utils.exceptions import InternalError, InvalidArgument, UserError
from .methods import clu_list, dim_list, lbt_list, ssclu_list, vis_list
from .operations import (clu_filter, dim_reduce_filter, lbt_filter,
                         ssclu_filter, vis_filter)
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification
from scipy.sparse.linalg.eigen.arpack import ArpackNoConvergence
from .cellar.utils.colors import PALETTE
from .data_loader import _get_singler_warning


class Signal(int, Enum):
    DIM_REDUCE = 101
    CLUSTER = 201
    MERGE = 250
    SS_CLUSTER = 251
    PALETTE = 280
    LABEL_TRANSFER = 301
    FEATURE_EXP = 401
    ANNOTATION = 501
    DATA_LOAD = 801
    RESET = 901


@app.callback(
    Output("main-plot-signal-code", "data"),
    Output("side-plot-signal-code", "data"),
    Output("singler-runtime-modal", "is_open"),
    Output("singler-runtime-body", "children"),
    MultiplexerOutput("push-notification", "data"),

    Input("dim-run-btn", "n_clicks"),
    Input("label-run-btn", "n_clicks"),
    Input("main-expression", "n_clicks"),
    Input("side-expression", "n_clicks"),
    Input("main-clear-expression-btn", "n_clicks"),
    Input("side-clear-expression-btn", "n_clicks"),
    Input("annotation-signal", "data"),
    Input("merge-plot-signal", "data"),
    Input("data-loaded-plot-signal", "data"),
    Input("data-loaded-plot-signal-prep", "data"),
    Input("data-loaded-plot-signal-prep-atac", "data"),
    Input("main-apply-palette-signal", "data"),
    Input("side-apply-palette-signal", "data"),
    Input("singler-runtime-yes", "n_clicks"),
    Input("singler-runtime-no", "n_clicks"),

    State("active-plot", "data"),
    State("label-tabs", "active_tab"),
    prevent_initial_call=True
)
def signal_plot(
        n1, n2, mexp, sexp, c1, c2, ans, mps, dlps, dlpsp, dlpspa, mapp, sapp,
        singlery, singlerno, actp, actt):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = 'a1' if actp == 1 else 'a2'
    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    to_return = [dash.no_update] * 5
    index = 0 if actp == 1 else 1
    button_id = ctx.triggered[0]["prop_id"].split(".")[0]

    # logger.info(f"Running button with id {button_id}.")
    if button_id == "dim-run-btn":
        to_return[index] = Signal.DIM_REDUCE
    elif button_id == "label-run-btn":
        if actt == 'clu':
            to_return[index] = Signal.CLUSTER
        elif actt == 'ssclu':
            to_return[index] = Signal.SS_CLUSTER
        elif actt == 'lbt':
            if 'a1' not in dbroot.adatas or 'a2' not in dbroot.adatas:
                to_return[-1] = _prep_notification(
                    "Label Transfer can only be run if a reference dataset "
                    "has been loaded under Dual Mode.", "warning")
                return to_return
            if dbroot.adatas[an]['adata'].shape[0] >= 10_000:
                to_return[2] = True
                nos, nof = dbroot.adatas[an]['adata'].shape
                to_return[3] = _get_singler_warning(nos, nof)
                return to_return
            to_return[index] = Signal.LABEL_TRANSFER
        else:
            raise InternalError(f"No tab with name {actt} found.")
    elif button_id == "main-expression" or button_id == "side-expression":
        to_return[index] = Signal.FEATURE_EXP
    elif button_id == "main-clear-expression-btn" or \
            button_id == "side-clear-expression-btn":
        to_return[index] = Signal.RESET
    elif button_id == "annotation-signal":
        if ans is not None:
            to_return[index] = Signal.ANNOTATION
    elif button_id == "merge-plot-signal":
        if mps is not None:
            to_return[index] = Signal.MERGE
    elif button_id == "data-loaded-plot-signal"\
            or button_id == "data-loaded-plot-signal-prep"\
            or button_id == "data-loaded-plot-signal-prep-atac":
        if dlps is not None or dlpsp is not None or dlpspa is not None:
            to_return[index] = Signal.DATA_LOAD
    elif button_id == "main-apply-palette-signal" or \
            button_id == "side-apply-palette-signal":
        to_return[index] = Signal.PALETTE
    elif button_id == "singler-runtime-yes":
        to_return[2] = False
        to_return[index] = Signal.LABEL_TRANSFER
    elif button_id == "singler-runtime-no":
        to_return[2] = False

    return to_return


def get_update_plot_func(an, prefix):
    def update_plot(
            s_code,
            dim_method, vis_method, clu_method, ssclu_method, lbt_method,
            feature_list, feature_range,
            *settings):
        ctx = dash.callback_context
        if not ctx.triggered or s_code is None:
            return dash.no_update

        title = dbroot.adatas[an]['name']

        # logger.info(f"Code {s_code}.")

        try:
            if s_code == Signal.DIM_REDUCE:
                # Reduce dimensions and prepare figure
                dim_reduce_filter(
                    dbroot.adatas[an]['adata'], dim_method, settings)
                clear_x_emb_dependends(dbroot.adatas[an]['adata'])
                vis_filter(
                    dbroot.adatas[an]['adata'], vis_method, settings)
                return get_dim_figure(dbroot.adatas[an]['adata'], title),\
                    dash.no_update, dash.no_update
            elif s_code == Signal.CLUSTER:
                # Cluster and prepare figure
                clu_filter(dbroot.adatas[an]['adata'], clu_method, settings)
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]), 1, dash.no_update
            elif s_code == Signal.SS_CLUSTER:
                ssclu_filter(dbroot.adatas[an]['adata'], ssclu_method,
                             settings, extras={'actp': 1 if an == 'a1' else 2})
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]), 1, dash.no_update
            elif s_code == Signal.LABEL_TRANSFER:
                ref_an = 'a2' if an == 'a1' else 'a1'
                lbt_filter(dbroot.adatas[an]['adata'], lbt_method, settings,
                           extras={'ref': dbroot.adatas[ref_an]['adata']})
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]), 1, dash.no_update
            elif s_code == Signal.PALETTE:
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]),\
                    dash.no_update, dash.no_update
            elif s_code == Signal.FEATURE_EXP:
                # Show gene expression levels and prepare figure
                if feature_list is None:
                    raise PreventUpdate
                exp = get_expression_figure(
                    dbroot.adatas[an]['adata'], feature_list, feature_range)
                return exp, dash.no_update, dash.no_update
            elif s_code == Signal.RESET:
                return get_reset_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]),\
                    dash.no_update, dash.no_update
            elif s_code == Signal.ANNOTATION:
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]),\
                    dash.no_update, dash.no_update
            elif s_code == Signal.MERGE:
                return get_clu_figure(
                    dbroot.adatas[an]['adata'], title,
                    palette=dbroot.palettes[prefix]), 1, dash.no_update
            elif s_code == Signal.DATA_LOAD:
                if 'x_emb_2d' in dbroot.adatas[an]['adata'].obsm:
                    if 'labels' in dbroot.adatas[an]['adata'].obs:
                        notif = _prep_notification(
                            "Found embeddings and labels.", icon="info")

                        return get_clu_figure(
                            dbroot.adatas[an]['adata'], title,
                            palette=dbroot.palettes[prefix]), 1, notif
                    notif = _prep_notification(
                        "Found embeddings on file.", icon="info")

                    return get_dim_figure(dbroot.adatas[an]['adata'], title),\
                        1, notif
                notif = _prep_notification(
                    "No embeddings found. Please, run preprocessing first if " +
                    "needed and then reduce dimensions.",
                    icon="info")

                return empty_figure, 1, notif
            else:
                raise InternalError(f"No signal with id {s_code} found.")
        except ArpackNoConvergence as anc:
            logger.error(str(anc))
            error_msg = "No eigenvectors converged. Try a different method."
            logger.warn(error_msg)
            return [dash.no_update] * 2 + [_prep_notification(
                error_msg, "warning")]
        except (InternalError, InvalidArgument, UserError) as e:
            logger.error(str(e))
            error_msg = str(e)
            return [dash.no_update] * 2 + [_prep_notification(
                error_msg, "danger")]
        except Exception as e:
            logger.error(str(e))
            error_msg = "Could not process request."
            logger.error(error_msg)
            return [dash.no_update] * 2 + [_prep_notification(
                error_msg, "danger")]

    return update_plot


# Loop over both plots
for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-plot", "figure"),
        Output(prefix + "-cluster-list-signal", "data"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-plot-signal-code", "data"),

        State("dim-methods-select", "value"),
        State("vis-methods-select", "value"),
        State("clu-methods-select", "value"),
        State("ssclu-methods-select", "value"),
        State("lbt-methods-select", "value"),
        State(prefix + "-feature-list", "value"),
        State(prefix + "-feature-rangeslider", "value"),
        [State(m['value'] + '-settings', 'children') for m in dim_list],
        [State(m['value'] + '-settings', 'children') for m in clu_list],
        [State(m['value'] + '-settings', 'children') for m in vis_list],
        [State(m['value'] + '-settings', 'children') for m in ssclu_list],
        [State(m['value'] + '-settings', 'children') for m in lbt_list],
        prevent_initial_call=True
    )(get_update_plot_func(an, prefix))


def get_download_plot_func(an):
    def _func(n1, figure, im_format, width, height, scale):

        try:
            width = int(width)
            height = int(height)
            scale = float(scale)

            filepath = 'tmp/plot.' + im_format

            plotly.io.write_image(
                figure, filepath, format=im_format,
                width=width, height=height, scale=scale)
        except Exception as e:
            logger.error(str(e))
            error_msg = "Error occurred when exporting plot."
            logger.error(error_msg)
            return dash.no_update, _prep_notification(error_msg, "danger")

        return dcc.send_file(filepath), dash.no_update

    return _func


for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + "-download-plot-buf", "data"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-plot-download-btn", "n_clicks"),

        State(prefix + "-plot", "figure"),
        State(prefix + "-plot-download-format", "value"),
        State(prefix + "-plot-download-width", "value"),
        State(prefix + "-plot-download-height", "value"),
        State(prefix + "-plot-download-scale", "value"),
        prevent_initial_call=True
    )(get_download_plot_func(an))


def get_clear_palette_func():
    def _func(n1):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        return [1] + [""] * len(PALETTE)
    return _func


for prefix in ['main', 'side']:
    app.callback(
        Output(prefix + "-clear-palette-signal", "data"),
        [Output(prefix + f'-color-input-{i}', "value")
         for i in range(len(PALETTE))],
        Input(prefix + "-clear-palette-btn", "n_clicks"),
        prevent_initial_call=True
    )(get_clear_palette_func())


def get_palette_bg_switcher_func(prefix):
    def _func(clear, *colors):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]

        if button_id == prefix + "-clear-palette-signal":
            return [{'background-color': c} for c in PALETTE]

        color_no = int(button_id.split('-')[-1])
        match = re.search(r'^(?:[0-9a-fA-F]{6})$', colors[color_no])

        to_return = [dash.no_update] * len(colors)
        if colors[color_no] == "":
            # dbroot.palettes[prefix][color_no] = PALETTE[color_no]
            to_return[color_no] = {'background-color': PALETTE[color_no]}
        elif match:
            # dbroot.palettes[prefix][color_no] = '#' + colors[color_no]
            to_return[color_no] = {'background-color': '#' + colors[color_no]}
        return to_return

    return _func


for prefix in ['main', 'side']:
    app.callback(
        [Output(prefix + f'-color-bg-{i}', "style")
         for i in range(len(PALETTE))],
        Input(prefix + "-clear-palette-signal", "data"),
        [Input(prefix + f'-color-input-{i}', "value")
         for i in range(len(PALETTE))],
        prevent_initial_call=True
    )(get_palette_bg_switcher_func(prefix))


def get_apply_palette_func(prefix):
    def _func(n1, *colors):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        new_pal = PALETTE.copy()

        for i in range(len(PALETTE)):
            if colors[i] is None or colors[i] == "":
                continue

            match = re.search(r'^(?:[0-9a-fA-F]{6})$', colors[i])
            if match:
                new_pal[i] = '#' + colors[i]
            else:
                # Only apply if not errors
                return dash.no_update, dash.no_update

        dbroot.palettes[prefix] = new_pal
        return 1, False
    return _func


for prefix in ['main', 'side']:
    app.callback(
        Output(prefix + "-apply-palette-signal", "data"),
        Output(prefix + "-color-palette-popover", "is_open"),
        Input(prefix + "-apply-palette-btn", "n_clicks"),
        [State(prefix + f'-color-input-{i}', "value")
         for i in range(len(PALETTE))],
        prevent_initial_call=True
    )(get_apply_palette_func(prefix))
