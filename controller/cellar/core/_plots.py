import numpy as np
import dash_bio as dashbio
import plotly.graph_objects as go
import plotly.express as px
import plotly.colors as pc

from ..utils.exceptions import InternalError, UserError
from ._tools import cl_add_gene_symbol, cl_get_expression
from ..utils.misc import is_sparse, get_title_from_feature_list
from ..utils.colors import PALETTE, get_col_i


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


def get_expression_figure(adata, feature_values, feature_range):
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
    """
    feature_values = np.array(feature_values, dtype='U200').flatten()

    for feat in feature_values:
        if feat not in adata.var_names.to_numpy():
            raise InternalError(f"Feature {feat} not found in var_names.")

    if 'x_emb_2d' not in adata.obsm:
        raise InternalError("No 2d embeddings found.")

    # Construct title with gene symbols
    # If more than 3 genes were used, we append dots
    title = get_title_from_feature_list(adata, feature_values)
    hover_data = {}

    if 'annotations' in adata.obs:
        hover_data['Annotation'] = adata.obs['annotations'].to_numpy()
    if 'labels' in adata.obs:
        hover_data['Cluster ID'] = adata.obs['labels'].to_numpy()

    expression = cl_get_expression(adata, feature_values)

    single_feature = len(feature_values) == 1
    if single_feature:
        if feature_range[0] > feature_range[1]:
            raise InternalError("Incorrect feature range found.")

        hover_data['Normalized Val.'] = expression.copy()
        new_min = np.min(expression[expression >= feature_range[0]])
        new_max = np.max(expression[expression <= feature_range[1]])
        expression[expression < feature_range[0]] = new_min
        expression[expression > feature_range[1]] = new_max

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


def get_heatmap(adata, feature_list):
    """
    Returns a plotly figure with a heatmap of gene expression values
    in feature_list. Uses dashbio.Clustergram.
    """
    if 'labels' not in adata.obs:
        raise InternalError("No labels found in adata.")

    feature_list = np.array(feature_list).flatten()

    if not set(feature_list).issubset(set(adata.var_names.to_numpy())):
        raise UserError("Some features not found in data.")

    unq_labels = np.unique(adata.obs['labels'])
    aves = []

    for label in unq_labels:
        row = adata[:, feature_list].X[np.where(
            adata.obs['labels'] == label)]
        if is_sparse(row):
            row = np.asarray(row.todense())
        row = row.mean(axis=0)
        aves.append(row)
    aves = np.array(aves)

    if 'gene_symbols' not in adata.var:
        cl_add_gene_symbol(adata)

    cluster = 'all'
    if aves.shape[1] == 1:
        cluster = 'row'
    if aves.shape[0] == 1:
        cluster = 'col'

    fig = dashbio.Clustergram(
        data=aves,
        column_labels=list(adata[:, feature_list].var['gene_symbols']),
        row_labels=['Cluster' + str(i) for i in unq_labels],
        color_map=pc.get_colorscale('Magma'),
        cluster=cluster,
        center_values=True,
        line_width=2,
        width=600
    )

    fig.update_xaxes(tickangle=-90)

    return fig


def get_violin_plot(adata, feature_values, feature_range, palette=PALETTE):
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
    palette: list
        Color palette to use.
    """
    for feature in feature_values:
        if feature not in adata.var_names:
            raise UserError(f"Feature {feature} not found in adata.")

    title = get_title_from_feature_list(adata, feature_values)
    vect = cl_get_expression(adata, feature_values)
    single_feature = len(feature_values) == 1
    if single_feature:
        if feature_range[0] > feature_range[1]:
            raise InternalError("Incorrect feature range found.")

        # hover_data['Normalized Val.'] = vect.copy()
        new_min = np.min(vect[vect >= feature_range[0]])
        new_max = np.max(vect[vect <= feature_range[1]])
        vect[vect < feature_range[0]] = new_min
        vect[vect > feature_range[1]] = new_max

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
