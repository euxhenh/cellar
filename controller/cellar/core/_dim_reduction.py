from sklearn.decomposition import PCA, TruncatedSVD, KernelPCA
from sklearn.manifold import MDS, TSNE
from pydiffmap import diffusion_map as dm
from umap import UMAP
from anndata._core.sparse_dataset import SparseDataset
from scipy.sparse import issparse


func_map = {
    'cl_PCA': PCA,
    'cl_TruncatedSVD': TruncatedSVD,
    'cl_kPCA': KernelPCA,
    'cl_MDS': MDS,
    'cl_UMAP': UMAP,
    'cl_TSNE': TSNE,
    'cl_Diffmap': dm.DiffusionMap.from_sklearn
}


def get_func(func_name):
    def _func(adata, key, x_to_use, **kwargs):
        for k in kwargs:
            if kwargs[k] == '':
                kwargs[k] = None

        if x_to_use == 'x':
            x_to_use = adata.X
            # Load sparse matrix to memory since cannot work with
            # HDF5 in backed mode
            if isinstance(adata.X, SparseDataset) or issparse(adata.X):
                if adata.isbacked:
                    x_to_use = x_to_use.to_memory()
        else:
            x_to_use = adata.obsm['x_emb']

        comp_key = 'n_evecs' if func_name == 'cl_Diffmap' else 'n_components'
        if comp_key not in kwargs:
            kwargs[comp_key] = 2

        fitter = func_map[func_name](**kwargs)
        adata.obsm[key] = fitter.fit_transform(x_to_use)

    return _func


for func_name in func_map.keys():
    globals()[func_name] = get_func(func_name)


def clear_x_emb_dependends(adata):
    if 'x_emb_2d' in adata.obsm:
        adata.obsm.pop('x_emb_2d')
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')
