import io
import os
import shutil
import tarfile
from base64 import b64decode

import plotly.express as px
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from controller.cellar.utils.exceptions import IncorrectFileFormat
from controller.cellar.utils.tile_generator import (generate_10x_spatial,
                                                    generate_tile)


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
            logger.info("Extracting tar.gz file.")
            tar = tarfile.open(fileobj=io.BytesIO(decoded))
            if os.path.isdir(extract_path):
                shutil.rmtree(extract_path)
            tar.extractall(extract_path)
            tar.close()
        else:
            raise IncorrectFileFormat(f'{filename} format not recognized.')

        return {}

    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-buf-load", "style"),

        Input(prefix + '-upload-spatial', 'contents'),
        State(prefix + '-upload-spatial', 'filename'),
        State(prefix + "-spatial-type-dropdown", "value"),
        prevent_initial_call=True
    )(get_parse_tar_gz_func(an))


def get_generate_tile_func(an):
    def _func(n1, data_type):
        if an not in dbroot.adatas:
            raise PreventUpdate

        if (data_type == 'spatial-10x'):
            # if not os.path.isdir(f'tmp/{an}/s10x'):
            #     raise PreventUpdate

            tile = generate_10x_spatial(
                f'tmp/{an}/s10x/detected_tissue_image.jpg',
                f'tmp/{an}/s10x/tissue_positions_list.csv',
                f'tmp/{an}/s10x/scalefactors_json.json',
                adata=dbroot.adatas['a1']['adata'],
                in_tissue=True)
        elif data_type == 'spatial-codex':
            if not os.path.isdir(f'tmp/{an}/codex'):
                raise PreventUpdate

            tile = generate_tile(
                f'tmp/{an}/codex/images',
                f'tmp/{an}/codex/data.csv',
                adata=dbroot.adatas[an]['adata'])
        else:
            raise PreventUpdate

        logger.info(f"Generated tile with shape {tile.shape}. Scaling...")

        ho, wo = tile.shape[:2]
        scaler = 1000 / max(wo, ho)
        w, h = scaler * wo, scaler * ho

        fig = px.imshow(tile, width=w, height=h)
        fig.update_layout(coloraxis_showscale=False)
        fig.update_xaxes(showticklabels=False)
        fig.update_yaxes(showticklabels=False)

        return fig
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-tile", "figure"),

        Input(prefix + "-generate-tile-btn", "n_clicks"),
        State(prefix + "-spatial-type-dropdown", "value"),
        prevent_initial_call=True
    )(get_generate_tile_func(an))
