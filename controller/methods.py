from .cellar.core import cl_PCA, cl_TruncatedSVD
from .cellar.core import cl_Leiden

dim_list = [
    {'label': 'PCA', 'value': 'dim-PCA', 'func': cl_PCA},
    {'label': 'Truncated SVD', 'value': 'dim-Truncated-SVD',
     'func': cl_TruncatedSVD},
    # {'label': 'UMAP', 'value': 'dim-UMAP'},
    # {'label': 'Diffusion Map', 'value': 'dim-Diffusion-Map'},
    # {'label': 'MDS', 'value': 'dim-MDS'},
    # {'label': 'Isomap', 'value': 'dim-Isomap'},
]


vis_list = [
    {'label': 'UMAP', 'value': 'vis-UMAP'},
    {'label': 'TSNE', 'value': 'vis-TSNE'},
    {'label': 'PCA', 'value': 'vis-PCA'},
    {'label': 'Truncated SVD', 'value': 'vis-Truncated-SVD'},
    {'label': 'Diffusion Map', 'value': 'vis-Diffusion-Map'},
    {'label': 'MDS', 'value': 'vis-MDS'},
    {'label': 'Isomap', 'value': 'vis-Isomap'},
]

clu_list = [
    {'label': 'Leiden', 'value': 'clu-Leiden', 'func': cl_Leiden},
    # {'label': 'KMeans', 'value': 'clu-KMeans'},
    # {'label': 'KMedoids', 'value': 'clu-KMedoids'},
    # {'label': 'Spectral Clustering', 'value': 'clu-Spectral-Clustering'},
    # {'label': 'Agglomerative Clustering',
    #  'value': 'clu-Agglomerative-Clustering'},
    # {'label': 'Cluster Ensemble', 'value': 'clu-Cluster-Ensemble'},
]


ssclu_list = [
    {'label': 'Constrained Leiden', 'value': 'ssclu-Constrained-Leiden'},
    {'label': 'Constrained KMeans', 'value': 'ssclu-Constrained-KMeans'},
    {'label': 'Seeded KMeans', 'value': 'ssclu-Seeded-KMeans'},
    {'label': 'KNN Filter', 'value': 'ssclu-KNN-Filter'}
]


lbt_list = [
    {'label': 'SingleR', 'value': 'lbt-SingleR'},
    {'label': 'Scanpy Ingest', 'value': 'lbt-Scanpy-Ingest'}
]


def find_method(cat_list, method_value):
    for c in cat_list:
        if c['value'] == method_value:
            return c
