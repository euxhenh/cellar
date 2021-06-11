import os
from base64 import b64encode, b64decode
import tarfile
import io

import dash
import dash_html_components as html
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
import plotly.express as px

from app import app, logger, dbroot
from controller.cellar.utils.tile_generator import generate_tile


def parse_tar_gz(contents, filename):
    content_type, content_string = contents.split(',')
    decoded = b64decode(content_string)

    if filename.endswith('tar.gz'):
        tar = tarfile.open(fileobj=io.BytesIO(decoded))
        tar.extractall('tmp')
        tar.close()


@app.callback(
    Output("load-dataset-btn", "style"),

    Input('main-upload-spatial', 'contents'),
    State('main-upload-spatial', 'filename'),
    prevent_initial_call=True
)
def upload_spatial(contents, filename):
    parse_tar_gz(contents, filename)


@app.callback(
    Output("main-tile", "children"),

    Input("main-generate-tile-btn", "n_clicks"),
    prevent_initial_call=True
)
def gen(n1):
    tile = generate_tile('tmp/images', 'tmp/data.csv',
                         adata=dbroot.adatas['a1']['adata'])
    fig = px.imshow(tile, width=700)
    fig.update_layout(coloraxis_showscale=False)
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)

    img_bytes = fig.to_image(format="png")
    encoding = b64encode(img_bytes).decode()
    img_b64 = "data:image/png;base64," + encoding
    return html.Img(src=img_b64, style={'width': '100%'})

    return fig
