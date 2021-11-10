from enum import auto

import dash_bio as dashbio
import numpy as np
import plotly.colors as pc
import plotly.express as px
import plotly.graph_objects as go

from ..utils.colors import PALETTE, get_col_i
from ..utils.exceptions import InternalError, UserError
from ..utils.misc import (_clip_2_range, _filter_outliers,
                          _gene_value_2_symbol, get_title_from_feature_list,
                          get_title_from_other_list)
from ._tools import (_collect_x_from_other, _collect_x_from_vars,
                     cl_add_gene_symbol, cl_get_expression)


def get_dim_figure(adata, title):
    """
    Returns a plotly figure showing a scatter plot of the 2D embeddings.
    """
    if 'x_emb_2d' not in adata.obsm:
        raise InternalError("x_emb_2d not found in adata.")

    indices = np.arange(adata.shape[0]).astype('str')
    names = adata.obs.index.to_numpy().astype('str')
    hover_name = [i + ": " + j for i, j in zip(indices, names)]

    fig = px.scatter(
        x=adata.obsm['x_emb_2d'][:, 0],
        y=adata.obsm['x_emb_2d'][:, 1],
        opacity=0.8,
        custom_data=[np.arange(adata.shape[0])],  # indices
        hover_name=hover_name,
        render_mode='webgl')  # for speed

    fig.update_layout(
        title=title,
        showlegend=True,
        plot_bgcolor='#FFFFFF',
        dragmode='lasso',
        xaxis={'visible': False},
        yaxis={'visible': False}
    )

    return fig


def get_clu_figure(adata, title, palette=PALETTE):
    """
    Returns a plotly figure showing a scatter plot of the 2D embeddings along
    with the cluster assignments. Will also parse 'annotations'
    if they are present in adata.
    """
    if 'x_emb_2d' not in adata.obsm or 'labels' not in adata.obs:
        raise InternalError("x_emb_2d or labels not found in adata.")

    hover_data = {}

    unq_colors = np.unique(adata.obs['labels'].to_numpy())
    maxval = np.max(unq_colors)
    unq_colors = list(unq_colors.astype(str))
    pal = [px.colors.convert_colors_to_same_type(i)[0][0] for i in palette]

    if 'annotations' in adata.obs:
        hover_data['Annotation'] = adata.obs['annotations'].to_numpy()

    indices = np.arange(adata.shape[0]).astype('str')
    names = adata.obs.index.to_numpy().astype('str')
    hover_name = [i + ": " + j for i, j in zip(indices, names)]

    fig = px.scatter(
        x=adata.obsm['x_emb_2d'][:, 0],
        y=adata.obsm['x_emb_2d'][:, 1],
        hover_name=hover_name,
        hover_data=hover_data,
        custom_data=[np.arange(adata.shape[0])],  # indices
        color=adata.obs['labels'].astype('str'),
        category_orders={'color': unq_colors},
        opacity=0.8,
        labels={'color': 'Cluster ID'},
        color_discrete_map={str(i): get_col_i(pal, i)
                            for i in range(maxval + 1)},
        render_mode='webgl')

    fig.update_layout(
        title=title,
        showlegend=True,
        plot_bgcolor='#FFFFFF',
        dragmode='lasso',
        xaxis={'visible': False},
        yaxis={'visible': False}
    )

    return fig


def get_reset_figure(adata, title, palette=PALETTE):
    if 'labels' in adata.obs:
        return get_clu_figure(adata, title, palette=palette)

    return get_dim_figure(adata, title)


def _check_features(feature_values, adata):
    all_feats = adata.var_names.to_numpy().astype(str)
    if not set(list(feature_values)).issubset(set(list(all_feats))):
        raise InternalError(f"Some features not found in adata.")


def _check_other(other_values, adata):
    prefix = other_values[0].split(':')[0]
    other_values = [i.split(':')[1] for i in other_values]
    all_feats = np.array(adata.uns[prefix]).astype('U200')
    if not set(list(other_values)).issubset(set(list(all_feats))):
        raise InternalError(f"Some features not found in adata.")


def get_expression_figure(
        adata, feature_values=None, feature_range=None,
        other_values=None, auto_scale=True):
    """
    Returns a plotly figure with the expression level of a gene
    or list of genes. In case a list of genes is used, then the
    expression value is calculated as shown in
    https://jingtaowang22.github.io/cellar_docs/#/analysis
    under "Feature Visualization".

    Parameters
    __________
    adata: anndata.AnnData object
    feature_values: list or str
        List of gene names
    feature_range: (float, float)
        Contains upper and lower thresholds to use for filtering
        high and low expression values. Useful when filtering outliers.
    other_values: same as above. used for other features.
    """
    if 'x_emb_2d' not in adata.obsm:
        raise InternalError("No 2d embeddings found.")

    if feature_values is not None:
        feature_values = np.array(feature_values, dtype='U200').flatten()
        _check_features(feature_values, adata)
        title = get_title_from_feature_list(adata, feature_values)
        expression = cl_get_expression(adata, feature_values)
        single_feature = len(feature_values) == 1
    else:
        other_values = np.array(other_values, dtype='U200').flatten()
        _check_other(other_values, adata)
        title = get_title_from_other_list(adata, other_values)
        expression = cl_get_expression(adata, other_names=other_values)
        single_feature = len(other_values) == 1
    # Construct title with gene symbols
    # If more than 3 genes were used, we append dots
    hover_data = {}
    if 'annotations' in adata.obs:
        hover_data['Annotation'] = adata.obs['annotations'].to_numpy()
    if 'labels' in adata.obs:
        hover_data['Cluster ID'] = adata.obs['labels'].to_numpy()

    if single_feature and not auto_scale:
        expression = _clip_2_range(expression, feature_range)
        hover_data['Normalized Val.'] = expression.copy()
    elif auto_scale:
        expression = _filter_outliers(expression)

    exp_max = expression.max()

    fig = px.scatter(
        x=adata.obsm['x_emb_2d'][:, 0],
        y=adata.obsm['x_emb_2d'][:, 1],
        hover_name=adata.obs.index.to_numpy().astype('str'),
        hover_data=hover_data,
        custom_data=[np.arange(adata.shape[0])],  # indices
        color=expression,
        opacity=0.8,
        labels={
            'color': 'Clipped Val.' if single_feature else 'Min Co-Exp.'},
        range_color=[expression.min(), 1 if exp_max == 0 else exp_max],
        render_mode='webgl',
        color_continuous_scale=px.colors.sequential.Magma
    )

    fig.update_layout(
        title=title,
        showlegend=True,
        plot_bgcolor='#FFFFFF',
        dragmode='lasso',
        xaxis={'visible': False},
        yaxis={'visible': False}
    )

    return fig


def get_heatmap(adata, feature_list=None, other_list=None):
    """
    Returns a plotly figure with a heatmap of gene expression values
    in feature_list. Uses dashbio.Clustergram.
    """
    if 'labels' not in adata.obs:
        raise InternalError("No labels found in adata.")

    if feature_list is not None:
        feature_list = np.array(feature_list).flatten()
        _check_features(feature_list, adata)
        x = _collect_x_from_vars(feature_list, adata)
        symbols = _gene_value_2_symbol(feature_list, adata)
    else:
        other_list = np.array(other_list).flatten()
        _check_other(other_list, adata)
        x = _collect_x_from_other(other_list, adata)
        symbols = [i.split(':')[1] for i in other_list]

    unq_labels = np.unique(adata.obs['labels'])
    aves = []

    for label in unq_labels:
        row = x[np.where(adata.obs['labels'] == label)]
        row = row.mean(axis=0)
        aves.append(row)
    aves = np.array(aves)

    cluster = 'all'
    if aves.shape[1] == 1:
        cluster = 'row'
    if aves.shape[0] == 1:
        cluster = 'col'

    fig = dashbio.Clustergram(
        data=aves,
        column_labels=list(symbols),
        row_labels=['Cluster' + str(i) for i in unq_labels],
        color_map=pc.get_colorscale('Magma'),
        cluster=cluster,
        center_values=True,
        line_width=2,
        width=600
    )

    fig.update_xaxes(tickangle=-90)

    return fig


def get_violin_plot(
        adata, feature_values=None, feature_range=None,
        other_values=None, palette=PALETTE, auto_scale=True):
    """
    Returns a plotly figure with a heatmap of gene expression values
    in feature_values. If multiple genes are selected, then the
    expression value is calculated as shown in
    https://jingtaowang22.github.io/cellar_docs/#/analysis
    under "Feature Visualization".

    Parameters
    __________
    adata: anndata.AnnData object
    feature_values: list or str
        List of gene names
    feature_range: (float, float)
        Contains upper and lower thresholds to use for filtering
        high and low expression values. Useful when filtering outliers.
    other_values: same as above. used for other features.
    palette: list
        Color palette to use.
    """
    if feature_values is not None:
        feature_values = np.array(feature_values, dtype='U200').flatten()
        _check_features(feature_values, adata)
        title = get_title_from_feature_list(adata, feature_values)
        vect = cl_get_expression(adata, feature_values)
        single_feature = len(feature_values) == 1
    else:
        other_values = np.array(other_values, dtype='U200').flatten()
        _check_other(other_values, adata)
        title = get_title_from_other_list(adata, other_values)
        vect = cl_get_expression(adata, other_names=other_values)
        single_feature = len(other_values) == 1

    if single_feature and not auto_scale:
        vect = _clip_2_range(vect, feature_range)
    elif auto_scale:
        vect = _filter_outliers(vect)

    unq_labels = np.unique(adata.obs['labels'])
    pal = [px.colors.convert_colors_to_same_type(i)[0][0] for i in palette]

    if 'gene_symbols' not in adata.var:
        cl_add_gene_symbol(adata)

    fig = go.Figure()

    for i, label in enumerate(unq_labels):
        fig.add_trace(
            go.Violin(
                y=vect[adata.obs['labels'] == label],
                name='Cluster ' + str(label),
                box_visible=True,
                meanline_visible=True,
                marker=go.violin.Marker(color=get_col_i(pal, label))
            )
        )

    fig.update_layout(
        title={'text': title},
        showlegend=False,
        plot_bgcolor='#FFFFFF',
        xaxis={'zeroline': False, 'showgrid': False},
        yaxis={'zeroline': False, 'showgrid': False},
    )

    return fig
