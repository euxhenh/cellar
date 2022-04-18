import io
import os
import shutil
import tarfile
from base64 import b64decode

import dash
from dash import dcc
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from .cellar.utils.tile_generator import (generate_10x_spatial,
                                          generate_tile)
from .cellar.utils.misc import get_title_from_feature_list
from .cellar.core import adjScoreProteinsCODEX, adjScoreClustersCODEX
from .cellar.core import adjScoreClusters10x
from .cellar.core import cl_get_expression
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification
from layout.misc import empty_spatial_figure, empty_colocalization_figure


def get_parse_tar_gz_func(an):
    def _func(contents, filename, data_type):
        print(data_type)
        if data_type == 'spatial-10x':
            extract_path = f'tmp/{an}/s10x'
        elif data_type == 'spatial-codex':
            extract_path = f'tmp/{an}/codex'
        else:
            return dash.no_update, _prep_notification("Please select a data type.", "info")

        if filename.endswith('tar.gz'):
            logger.info(f"Extracting tar.gz file at {extract_path}.")
            try:
                content_type, content_string = contents.split(',')
                decoded = b64decode(content_string)
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


def _determine_spatial_colors(adata, feature_list):
    """
    Determine whether to use labels or protein expression. If neither,
    will return None.
    """
    colors = None
    if feature_list is None or len(feature_list) == 0:
        if 'labels' in adata.obs:
            colors = adata.obs['labels'].to_numpy().astype(int)
    else:
        feature_list = np.array(feature_list, dtype='U200').flatten()
        colors = cl_get_expression(adata, feature_list).astype(float)
    return colors


def get_generate_tile_func(an, prefix):
    def _func(n1, clean, data_type, feature_list):
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

        adata = dbroot.adatas[an]['adata']
        colors = _determine_spatial_colors(adata, feature_list)
        savepath = f'tmp/spatial_{an}_tile.png'
        owner = None

        if data_type == 'spatial-10x':
            try:
                tile, owner = generate_10x_spatial(
                    f'tmp/{an}/s10x/spatial/detected_tissue_image.jpg',
                    f'tmp/{an}/s10x/spatial/tissue_positions_list.csv',
                    f'tmp/{an}/s10x/spatial/scalefactors_json.json',
                    adata=adata,
                    colors=colors,
                    in_tissue=False,
                    savepath=savepath,
                    palette=dbroot.palettes[prefix])
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred when generating 10x spatial tile."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        elif data_type == 'spatial-codex':
            tile_list = os.listdir('data/codex_tile')
            fname = dbroot.adatas[an]['name']
            # Selected a server CODEX dataset
            if fname in tile_list:
                logger.info("Found CODEX tile locally.")
                impath = f'data/codex_tile/{fname}/images'
                datapath = f'data/codex_tile/{fname}/data.csv'
            # In case no necessary spatial tile information has been uploaded
            elif not os.path.isdir(f'tmp/{an}/codex'):
                if 'x' not in adata.obs or 'y' not in adata.obs:
                    error_msg = "No spatial files have been uploaded."
                    logger.error(error_msg)
                    return dash.no_update, _prep_notification(error_msg, "warn")
                impath = None
                datapath = None
            else:
                impath = f'tmp/{an}/codex/images'
                datapath = f'tmp/{an}/codex/data.csv'

            try:
                tile, owner = generate_tile(
                    impath, datapath, adata=adata, colors=colors,
                    palette=dbroot.palettes[prefix], savepath=savepath)
            except Exception as e:
                logger.error(str(e))
                error_msg = "Error occurred when generating CODEX tile."
                logger.error(error_msg)
                return dash.no_update, _prep_notification(error_msg, "danger")
        else:
            msg = "Please select a data type."
            return dash.no_update, _prep_notification(msg, "info")

        logger.info(f"Generated tile with shape {tile.shape}. Scaling...")
        ho, wo = tile.shape[:2]
        scaler = 1000 / max(wo, ho)
        w, h = int(scaler * wo), int(scaler * ho)

        title = None
        if feature_list is not None and len(feature_list) >= 1:
            title = get_title_from_feature_list(adata, feature_list)

        fig = px.imshow(tile, width=w, height=h,
                        color_continuous_scale='magma',
                        binary_compression_level=9,
                        binary_format='jpg', title=title)

        if owner is not None and 'labels' in adata.obs:
            owner_cp = owner.copy()
            owner = owner.clip(min=0)
            customdata = adata.obs['labels'].to_numpy()[owner]
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
        State(prefix + "-feature-list-spatial", "value"),
        prevent_initial_call=True
    )(get_generate_tile_func(an, prefix))


def _remap_indices(x, old, new):
    sort_idx = old.argsort()
    mapped_idx = sort_idx[
        np.searchsorted(old, x, sorter=sort_idx)]
    remapped = new[mapped_idx]
    return remapped


def get_generate_cluster_scores_func(an, prefix):
    def _func(n1, clean, data_type, n_neighbors):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate
        adata = dbroot.adatas[an]['adata']

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return empty_colocalization_figure, dash.no_update

        if data_type == 'spatial-10x':
            if 'spatial_dict' in adata.uns:
                csv_path = None
            else:
                csv_path = f'tmp/{an}/s10x/spatial/tissue_positions_list.csv'

            res = adjScoreClusters10x(
                adata, csv_path, n_neighbors=n_neighbors)
        elif data_type == 'spatial-codex':
            if 'x' in adata.obs and 'y' in adata.obs:
                res = adjScoreClustersCODEX(
                    adata, None, n_neighbors=n_neighbors)
            else:
                tile_list = os.listdir('data/codex_tile')
                fname = dbroot.adatas[an]['name']

                if fname not in tile_list:
                    msg = "Could not find spatial data or data type " + \
                        "not supported."
                    logger.warn(msg)
                    return dash.no_update, _prep_notification(
                        msg, icon="warning")

                csv_path = f'data/codex_tile/{fname}/data.csv'
                res = adjScoreClustersCODEX(
                    adata, csv_path, n_neighbors=n_neighbors)
        else:
            return dash.no_update, _prep_notification(
                "Please select a data type.", icon="info")

        x_cord = res['f'].to_numpy().astype(int)
        y_cord = res['g'].to_numpy().astype(int)
        scores = res['score'].astype(float)
        q = np.round(res['q'].astype(float), 4)
        unq_labels = np.unique(dbroot.adatas[an]['adata'].obs['labels'])
        n = len(unq_labels)
        x_cord = _remap_indices(x_cord, unq_labels, np.arange(n))
        y_cord = _remap_indices(y_cord, unq_labels, np.arange(n))
        heatmap, qvals = np.zeros((n, n)), np.zeros((n, n))
        heatmap[x_cord, y_cord] = scores
        heatmap[y_cord, x_cord] = scores
        # heatmap = np.log10(heatmap + 1)
        qvals[x_cord, y_cord] = q
        qvals[y_cord, x_cord] = q
        heatmap = np.flip(heatmap, axis=1)
        qvals = np.flip(qvals, axis=1)

        fig = go.Figure(data=[go.Heatmap(
            x=np.repeat(unq_labels, n).astype(str),
            y=np.flip(np.tile(unq_labels, n).astype(str)),
            z=heatmap.flatten(),
            text=qvals.astype(str).flatten(),
            colorscale='magma',
            hovertemplate="Cluster x: %{x}<br>Cluster y: " +
            "%{y}<br>score: %{z}<br>q-val: %{text}"
        )])
        fig.update_layout(
            width=500,
            height=500,
            autosize=False,
            xaxis=dict(tickmode='linear'),
            yaxis=dict(tickmode='linear'),
            legend_title_text="score"
        )
        return fig, dash.no_update
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-cluster-scores", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-generate-cluster-scores-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        State(prefix + "-spatial-type-dropdown", "value"),
        State(prefix + "-spatial-cluster-scores-nneigh", "value"),
        prevent_initial_call=True
    )(get_generate_cluster_scores_func(an, prefix))


def get_generate_protein_scores_func(an, prefix):
    def _func(n1, clean, data_type, n_neighbors):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate
        if an not in dbroot.adatas:
            raise PreventUpdate
        if 'adata' not in dbroot.adatas[an]:
            raise PreventUpdate
        adata = dbroot.adatas[an]['adata']

        button_id = ctx.triggered[0]["prop_id"].split(".")[0]
        if button_id == prefix + "-data-load-clean":
            return empty_colocalization_figure, dash.no_update

        if 'x' in adata.obs and 'y' in adata.obs:
            res = adjScoreProteinsCODEX(
                adata, None, n_neighbors=n_neighbors)
        else:
            tile_list = os.listdir('data/codex_tile')
            fname = dbroot.adatas[an]['name']

            if fname not in tile_list:
                msg = "Could not find spatial data or data type not supported."
                logger.warn(msg)
                return dash.no_update, _prep_notification(msg, icon="warning")

            csv_path = f'data/codex_tile/{fname}/data.csv'
            res = adjScoreProteinsCODEX(
                dbroot.adatas[an]['adata'], csv_path, n_neighbors=n_neighbors)

        features = dbroot.adatas[an]['adata'].var['gene_symbols'].to_numpy()
        x_cord, y_cord = res['f'].astype(int), res['g'].astype(int)
        scores = res['score'].astype(float)

        n = features.shape[0]
        heatmap, qvals = np.zeros((n, n)), np.zeros((n, n))
        heatmap[x_cord, y_cord] = scores
        heatmap[y_cord, x_cord] = scores
        # heatmap = np.log10(heatmap + 1)
        q = np.round(res['q'].astype(float), 4)
        qvals[x_cord, y_cord] = q
        qvals[y_cord, x_cord] = q
        heatmap = np.flip(heatmap, axis=1)
        qvals = np.flip(qvals, axis=1)

        fig = go.Figure(data=[go.Heatmap(
            x=np.repeat(features, n).astype(str),
            y=np.flip(np.tile(features, n).astype(str)),
            z=heatmap.flatten(),
            text=qvals.astype(str).flatten(),
            colorscale='magma',
            hovertemplate="Cluster x: %{x}<br>Cluster y: " +
            "%{y}<br>score: %{z}<br>q-val: %{text}"
        )])
        fig.update_layout(
            width=500,
            height=500,
            autosize=False,
            xaxis=dict(tickmode='linear'),
            yaxis=dict(tickmode='linear'),
            legend_title_text="score"
        )
        return fig, dash.no_update
    return _func


for prefix, an in zip(["main", "side"], ["a1", "a2"]):
    app.callback(
        Output(prefix + "-protein-scores", "figure"),
        MultiplexerOutput("push-notification", "data"),

        Input(prefix + "-generate-protein-scores-btn", "n_clicks"),
        Input(prefix + "-data-load-clean", "data"),
        State(prefix + "-spatial-type-dropdown", "value"),
        State(prefix + "-spatial-protein-scores-nneigh", "value"),
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
