from .cellar.core import (cl_Agglomerative, cl_Ingest, cl_KMeans, cl_KMedoids,
                          cl_kPCA, cl_Leiden, cl_MDS, cl_PCA, cl_TSNE,
                          cl_Diffmap, cl_SingleR, cl_cisTopic, cl_ExactLT,
                          cl_SpectralClustering, cl_ssLeiden, cl_TruncatedSVD,
                          cl_UMAP, cl_uncertainty)

dim_list = [
    {'label': 'PCA', 'value': 'dim-PCA', 'func': cl_PCA},
    {'label': 'Truncated SVD', 'value': 'dim-Truncated-SVD',
     'func': cl_TruncatedSVD},
    {'label': 'Kernel PCA', 'value': 'dim-Kernel-PCA', 'func': cl_kPCA},
    {'label': 'Diffusion Map', 'value': 'dim-Diffmap', 'func': cl_Diffmap},
    {'label': 'MDS', 'value': 'dim-MDS', 'func': cl_MDS},
    {'label': 'UMAP', 'value': 'dim-UMAP', 'func': cl_UMAP},
    {'label': 'cisTopic', 'value': 'dim-cisTopic', 'func': cl_cisTopic}
]


vis_list = [
    {'label': 'UMAP', 'value': 'vis-UMAP', 'func': cl_UMAP},
    {'label': 't-SNE', 'value': 'vis-TSNE', 'func': cl_TSNE},
    {'label': 'PCA', 'value': 'vis-PCA', 'func': cl_PCA},
    {'label': 'Diffusion Map', 'value': 'vis-Diffmap', 'func': cl_Diffmap},
    {'label': 'MDS', 'value': 'vis-MDS', 'func': cl_MDS},
    {'label': 'Truncated SVD', 'value': 'vis-Truncated-SVD',
     'func': cl_TruncatedSVD},
    {'label': 'Kernel PCA', 'value': 'vis-Kernel-PCA', 'func': cl_kPCA}
]

clu_list = [
    {'label': 'Leiden', 'value': 'clu-Leiden', 'func': cl_Leiden},
    {'label': 'KMeans', 'value': 'clu-KMeans', 'func': cl_KMeans},
    {'label': 'KMedoids', 'value': 'clu-KMedoids', 'func': cl_KMedoids},
    {'label': 'Spectral Clustering', 'value': 'clu-Spectral-Clustering',
     'func': cl_SpectralClustering},
    {'label': 'Agglomerative Clustering',
     'value': 'clu-Agglomerative', 'func': cl_Agglomerative},
    {'label': 'Uncertainty Clustering',
        'value': 'clu-Uncertainty-Clustering', 'func': cl_uncertainty}
    # {'label': 'Cluster Ensemble', 'value': 'clu-Cluster-Ensemble','func':},
]


ssclu_list = [
    {'label': 'Constrained Leiden',
        'value': 'ssclu-Constrained-Leiden', 'func': cl_ssLeiden},

    # {'label': 'Constrained KMeans', 'value': 'ssclu-Constrained-KMeans'},
    # {'label': 'Seeded KMeans', 'value': 'ssclu-Seeded-KMeans'},
    # {'label': 'KNN Filter', 'value': 'ssclu-KNN-Filter'}
]


lbt_list = [
    {'label': 'Scanpy Ingest',
        'value': 'lbt-Scanpy-Ingest', 'func': cl_Ingest},
    {'label': 'SingleR', 'value': 'lbt-SingleR', 'func': cl_SingleR},
    {'label': 'Cell ID Based', 'value': 'lbt-exact', 'func': cl_ExactLT}
]


def find_method(cat_list, method_value):
    for c in cat_list:
        if c['value'] == method_value:
            return c
