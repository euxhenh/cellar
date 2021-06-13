import os
import shutil
from base64 import b64encode, b64decode
import tarfile
import io

import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px

from app import app, logger, dbroot
from controller.cellar.utils.exceptions import IncorrectFileFormat
from controller.cellar.utils.tile_generator import generate_tile


def get_parse_tar_gz_func(an):
    def _func(contents, filename):
        content_type, content_string = contents.split(',')
        decoded = b64decode(content_string)

        if filename.endswith('tar.gz'):
            tar = tarfile.open(fileobj=io.BytesIO(decoded))
            if os.path.isdir(f'tmp/{an}'):
                shutil.rmtree(f'tmp/{an}')
            tar.extractall(f'tmp/{an}')
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
        prevent_initial_call=True
    )(get_parse_tar_gz_func(an))


def get_generate_tile_func(an):
    def _func(n1):
        if not os.path.isdir(f'tmp/{an}'):
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate

        tile = generate_tile(
            f'tmp/{an}/images',
            f'tmp/{an}/data.csv',
            adata=dbroot.adatas[an]['adata'])

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
        prevent_initial_call=True
    )(get_generate_tile_func(an))
