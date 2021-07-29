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
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification


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


def get_generate_tile_func(an):
    def _func(n1, data_type):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate

        if (data_type == 'spatial-10x'):
            # if not os.path.isdir(f'tmp/{an}/s10x'):
            #     raise PreventUpdate
            try:
                tile = generate_10x_spatial(
                    f'tmp/{an}/s10x/spatial/detected_tissue_image.jpg',
                    f'tmp/{an}/s10x/spatial/tissue_positions_list.csv',
                    f'tmp/{an}/s10x/spatial/scalefactors_json.json',
                    adata=dbroot.adatas['a1']['adata'],
                    in_tissue=True)
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

                    tile = generate_tile(
                        f'data/codex_tile/{fname}/images',
                        f'data/codex_tile/{fname}/data.csv',
                        adata=dbroot.adatas[an]['adata'])
                else:
                    if not os.path.isdir(f'tmp/{an}/codex'):
                        raise PreventUpdate

                    tile = generate_tile(
                        f'tmp/{an}/codex/images',
                        f'tmp/{an}/codex/data.csv',
                        adata=dbroot.adatas[an]['adata'])
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred when generating CODEX tile."
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
        State(prefix + "-spatial-type-dropdown", "value"),
        prevent_initial_call=True
    )(get_generate_tile_func(an))


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
