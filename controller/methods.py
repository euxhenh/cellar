from .cellar.core import cl_PCA, cl_TruncatedSVD, cl_IncrementalPCA, cl_KernelPCA 
from .cellar.core import cl_TSNE, cl_Isomap, cl_MDS, cl_SpectralEmbedding, cl_FeatureAgglomeration
from .cellar.core import cl_dm, cl_umap
from .cellar.core import cl_Leiden

dim_list = [
    {'label': 'PCA', 'value': 'dim-PCA', 'func': cl_PCA},
    {'label': 'Truncated SVD', 'value': 'dim-Truncated-SVD',
     'func': cl_TruncatedSVD},
    {'label': 'Increamental PCA', 'value': 'dim-IPCA', 'func': cl_IncrementalPCA},
    {'label': 'Kernel PCA', 'value': 'dim-KPCA', 'func': cl_KernelPCA},
    {'label': 'TSNE', 'value': 'dim-TSNE', 'func': cl_TSNE},
    {'label': 'UMAP', 'value': 'dim-UMAP', 'func': cl_umap},
    {'label': 'Diffusion Map', 'value': 'dim-Diffusion-Map', 'func':cl_dm},
    {'label': 'MDS', 'value': 'dim-MDS', 'func': cl_MDS},
    {'label': 'Isomap', 'value': 'dim-Isomap', 'func': cl_Isomap},
    {'label': 'Spectral Embedding', 'value': 'dim-Spectral-Embedding', 'func': cl_SpectralEmbedding},
    {'label': 'Feature Agglomeration', 'value': 'dim-Feature-Agglomeration', 'func': cl_FeatureAgglomeration}
    
    
]


vis_list = [
    {'label': 'UMAP', 'value': 'vis-UMAP','func':cl_umap},
    {'label': 'TSNE', 'value': 'vis-TSNE', 'func':cl_TSNE},
    {'label': 'PCA', 'value': 'vis-PCA', 'func':cl_PCA},
    {'label': 'Truncated SVD', 'value': 'vis-Truncated-SVD','func':cl_TruncatedSVD},
    {'label': 'Diffusion Map', 'value': 'vis-Diffusion-Map' ,'func':cl_dm},
    {'label': 'MDS', 'value': 'vis-MDS' ,'func':cl_MDS},
    {'label': 'Isomap', 'value': 'vis-Isomap', 'func':cl_Isomap},
    {'label': 'Increamental PCA', 'value': 'vis-IPCA', 'func': cl_IncrementalPCA},
    {'label': 'Kernel PCA', 'value': 'vis-KPCA', 'func': cl_KernelPCA},
    {'label': 'Spectral Embedding', 'value': 'vis-Spectral-Embedding', 'func': cl_SpectralEmbedding},
    {'label': 'Feature Agglomeration', 'value': 'vis-Feature-Agglomeration', 'func': cl_FeatureAgglomeration}
    
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
