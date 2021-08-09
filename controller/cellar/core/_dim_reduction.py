import anndata2ri
import numpy as np
import rpy2.robjects as ro
from rpy2.robjects import r, numpy2ri
from rpy2.robjects.packages import importr
from anndata._core.sparse_dataset import SparseDataset
from pydiffmap import diffusion_map as dm
from scipy.sparse import issparse, csr_matrix
from sklearn.decomposition import PCA, KernelPCA, TruncatedSVD
from sklearn.manifold import MDS, TSNE
from umap import UMAP

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


def _correct_bin_names(bin_names):
    for i in range(len(bin_names)):
        bin_names[i] = bin_names[i].replace(
            ':', '_', bin_names[i].count(':')-1)
    return bin_names


def cl_cisTopic(
        adata, key, x_to_use, topics=40, iterations=150, **kwargs):
    topics = int(topics)
    topic_list = [topics, topics + 5, topics - 5]
    mat = adata.to_memory().X.T.copy()

    if issparse(mat):
        mat = mat.tocoo()
        r_Matrix = importr("Matrix")
        mat = r_Matrix.sparseMatrix(
            i=ro.IntVector(mat.row + 1),
            j=ro.IntVector(mat.col + 1),
            x=ro.FloatVector(mat.data),
            dims=ro.IntVector(mat.shape))
    else:
        mat = numpy2ri.py2rpy(mat)

    var_names = _correct_bin_names(adata.var_names.to_numpy())
    mat = r("`rownames<-`")(mat, ro.vectors.StrVector(var_names))
    mat = r("`colnames<-`")(mat, ro.vectors.StrVector(adata.obs.index))

    cisTopic = importr('cisTopic')

    print('Creating cisTopic object.')
    cc = cisTopic.createcisTopicObject(mat)
    print('Created cisTopic object.')

    print('Running LDA Models...')
    cc = cisTopic.runWarpLDAModels(
        cc,
        topic=numpy2ri.py2rpy(topic_list),
        nCores=min(4, len(topic_list)),
        iterations=iterations,
        addModels=False,
        returnType='selectedModel')

    cellassign = cisTopic.modelMatSelection(cc, 'cell', 'Probability')
    cellassign = np.array(cellassign).T.copy()  # Transpose back

    adata.obsm[key] = cellassign


def clear_x_emb_dependends(adata):
    if 'x_emb_2d' in adata.obsm:
        adata.obsm.pop('x_emb_2d')
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')
