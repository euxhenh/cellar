import io
import os
import shutil
import tarfile
from base64 import b64decode

import dash
import dash_core_components as dcc
import numpy as np
import matplotlib
import plotly.express as px
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .cellar.utils.tile_generator import (generate_10x_spatial,
                                          generate_tile)
from .cellar.core import adjScoreProteinsCODEX, adjScoreClustersCODEX
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification
from layout.misc import empty_spatial_figure, empty_colocalization_figure


def get_parse_tar_gz_func(an):
    def _func(contents, filename, data_type):
        content_type, content_string = contents.split(',')
        decoded = b64decode(content_string)

        if data_type == 'spatial-10x':
            extract_path = f'tmp/{an}/s10x'
        elif data_type == 'spatial-codex':
            extract_path = f'tmp/{an}/codex'
        else:
            raise PreventUpdate

        if filename.endswith('tar.gz'):
            logger.info(f"Extracting tar.gz file at {extract_path}.")
            try:
                tar = tarfile.open(fileobj=io.BytesIO(decoded))
                if os.path.isdir(extract_path):
                    shutil.rmtree(extract_path)
                tar.extractall(extract_path)
                tar.close()
            except Exception as e:
                logger.error(str(e))
                error_msg = "Couldn't extract tar.gz file."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        else:
            return dash.no_update, _prep_notification(
                f'{filename} format not recognized.', "danger")

        return {}, _prep_notification("Finished extracting file.", "info")

    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-buf-load", "style"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + '-upload-spatial', 'contents'),
        State(prefix + '-upload-spatial', 'filename'),
        State(prefix + "-spatial-type-dropdown", "value"),
        prevent_initial_call=True
    )(get_parse_tar_gz_func(an))


def get_generate_tile_func(an, prefix):
    def _func(n1, clean, data_type):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return empty_spatial_figure, dash.no_update

        owner = None
        if data_type == 'spatial-10x':
            # if not os.path.isdir(f'tmp/{an}/s10x'):
            #     raise PreventUpdate
            try:
                tile, owner = generate_10x_spatial(
                    f'tmp/{an}/s10x/spatial/detected_tissue_image.jpg',
                    f'tmp/{an}/s10x/spatial/tissue_positions_list.csv',
                    f'tmp/{an}/s10x/spatial/scalefactors_json.json',
                    adata=dbroot.adatas['a1']['adata'],
                    in_tissue=True,
                    palette=dbroot.palettes[prefix])
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred when generating 10x spatial tile."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        elif data_type == 'spatial-codex':
            tile_list = os.listdir('data/codex_tile')
            fname = dbroot.adatas[an]['name']

            try:
                if fname in tile_list:
                    logger.info("Found CODEX tile locally.")

                    tile, owner = generate_tile(
                        f'data/codex_tile/{fname}/images',
                        f'data/codex_tile/{fname}/data.csv',
                        adata=dbroot.adatas[an]['adata'],
                        palette=dbroot.palettes[prefix])
                else:
                    if not os.path.isdir(f'tmp/{an}/codex'):
                        raise PreventUpdate

                    tile, owner = generate_tile(
                        f'tmp/{an}/codex/images',
                        f'tmp/{an}/codex/data.csv',
                        adata=dbroot.adatas[an]['adata'],
                        palette=dbroot.palettes[prefix])
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred when generating CODEX tile. " +\
                    "Have the necessary supplementary files been uploaded?"
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        else:
            raise PreventUpdate

        if tile is None:
            error_msg = "No tile was generated."
            logger.error(error_msg)
            return dash.no_update, _prep_notification(error_msg, "danger")

        logger.info(f"Generated tile with shape {tile.shape}. Scaling...")

        matplotlib.image.imsave(
            f'tmp/spatial_{an}_tile.png', tile.astype(np.uint8))

        ho, wo = tile.shape[:2]
        scaler = 1000 / max(wo, ho)
        w, h = int(scaler * wo), int(scaler * ho)

        fig = px.imshow(tile, width=w, height=h)

        if owner is not None and 'labels' in dbroot.adatas[an]['adata'].obs:
            owner_cp = owner.copy()
            owner[owner < 0] = 0
            customdata = dbroot.adatas[an][
                'adata'].obs['labels'].to_numpy()[owner]
            customdata[owner_cp < 0] = -1
            fig.update(data=[{
                'customdata': customdata,
                'hovertemplate': 'Cluster ID: %{customdata}'}])
        fig.update_layout(coloraxis_showscale=False)
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig, _prep_notification(
            f"Generated tile with shape {tile.shape[:2]}", "info")
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-tile", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-generate-tile-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        State(prefix + "-spatial-type-dropdown", "value"),
        prevent_initial_call=True
    )(get_generate_tile_func(an, prefix))


def get_generate_cluster_scores_func(an, prefix):
    def _func(n1, clean):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return empty_colocalization_figure, dash.no_update

        tile_list = os.listdir('data/codex_tile')
        fname = dbroot.adatas[an]['name']

        if fname not in tile_list:
            msg = "Could not find spatial data."
            logger.warn()
            return dash.no_update, _prep_notification(msg, icon="warning")

        csv_path = f'data/codex_tile/{fname}/data.csv'
        res = adjScoreClustersCODEX(dbroot.adatas[an]['adata'], csv_path)

        x_cord, y_cord = res['f'].astype(int), res['g'].astype(int)
        scores = res['score'].astype(float)
        n = max(x_cord.max(), y_cord.max()) + 1
        heatmap = np.zeros((n, n))
        heatmap[x_cord, y_cord] = scores

        fig = px.imshow(
            heatmap,
            x=np.arange(n).astype(str),
            y=np.arange(n).astype(str),
            labels={
                'color': 'score'
            },
            color_continuous_scale='viridis'
        )
        # fig.update(data=[{
        #     'customdata': np.round(res['q'], 3),
        #     'hovertemplate': 'q-value: %{customdata}'}])
        return fig, dash.no_update
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-cluster-scores", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-generate-cluster-scores-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        prevent_initial_call=True
    )(get_generate_cluster_scores_func(an, prefix))


def get_generate_protein_scores_func(an, prefix):
    def _func(n1, clean):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return empty_colocalization_figure, dash.no_update

        tile_list = os.listdir('data/codex_tile')
        fname = dbroot.adatas[an]['name']

        if fname not in tile_list:
            msg = "Could not find spatial data."
            logger.warn()
            return dash.no_update, _prep_notification(msg, icon="warning")

        csv_path = f'data/codex_tile/{fname}/data.csv'
        res = adjScoreProteinsCODEX(dbroot.adatas[an]['adata'], csv_path)

        x_cord, y_cord = res['f'].astype(int), res['g'].astype(int)
        scores = res['score'].astype(float)
        n = max(x_cord.max(), y_cord.max()) + 1
        heatmap = np.zeros((n, n))
        heatmap[x_cord, y_cord] = scores

        fig = px.imshow(
            heatmap,
            x=np.arange(n).astype(str),
            y=np.arange(n).astype(str),
            labels={
                'color': 'score'
            },
            color_continuous_scale='viridis'
        )
        # fig.update(data=[{
        #     'customdata': np.round(res['q'], 3),
        #     'hovertemplate': 'q-value: %{customdata}'}])
        return fig, dash.no_update
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-protein-scores", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-generate-protein-scores-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        prevent_initial_call=True
    )(get_generate_protein_scores_func(an, prefix))


def get_download_tile_func(an):
    def _func(n1):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate

        filepath = f'tmp/spatial_{an}_tile.png'

        if not os.path.isfile(filepath):
            return dash.no_update, _prep_notification(
                "No tile image found.")

        return dcc.send_file(filepath), dash.no_update

    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-download-tile-buf", "data"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-download-tile-btn", "n_clicks"),
        prevent_initial_call=True
    )(get_download_tile_func(an))
