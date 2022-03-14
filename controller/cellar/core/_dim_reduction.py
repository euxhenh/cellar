import gc

import anndata2ri
import numpy as np
import rpy2.robjects as ro
from anndata._core.sparse_dataset import SparseDataset
from app import logger
from controller.cellar.utils.exceptions import InvalidArgument
from pydiffmap import diffusion_map as dm
from rpy2.robjects import numpy2ri, r
from rpy2.robjects.packages import importr
from scipy.sparse import csr_matrix, issparse
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
    """
    Use a functional to generate a function for each dimensionality reduction
    method since they all have the same interface.
    """
    method = func_name[3:]  # Parse the method name

    def _func(adata, key, x_to_use, **kwargs):
        """
        Reduces the dimensionality of the data using the 'func_name' method.

        Parameters
        __________

        adata: anndata.AnnData object
        key: str
            Key to store the reduced data under adata.obsm
        x_to_use: str
            Can be 'x' or 'x_emb'. If set to 'x', will use adata.X
            to reduce the data. Otherwise will use adata.obsm['x_emb'].
            We need the latter when this function is called to find 2D
            embeddings.
        kwargs: dict
            Any additional arguments passed to the constructor of func_name.
        """
        # Empty input boxes are parsed as empty strings
        for k in kwargs:
            if kwargs[k] == '':
                kwargs[k] = None

        if x_to_use == 'x':
            x_to_use = adata.X
            # Load sparse matrix to memory since cannot work with
            # HDF5 in backed mode
            if isinstance(adata.X, SparseDataset) or issparse(adata.X):
                if func_name not in ['cl_TruncatedSVD', 'cl_UMAP']:
                    raise InvalidArgument(
                        "Sparse data is not supported using the selected "
                        "reduction method. "
                        "Please choose TruncatedSVD or UMAP.")
                if adata.isbacked:
                    x_to_use = x_to_use.to_memory()
        else:
            x_to_use = adata.obsm['x_emb']

        # Diffusion maps use a different parameter name for the number of comp
        comp_key = 'n_evecs' if func_name == 'cl_Diffmap' else 'n_components'
        # If no number of components was found in kwargs, assume this
        # method was run for visualizing the data and set n_components to 2.
        if comp_key not in kwargs:
            kwargs[comp_key] = 2

        mins = min(x_to_use.shape[0], x_to_use.shape[1])
        if kwargs[comp_key] >= mins:
            raise InvalidArgument(
                "Number of components is higher than or equal to " +
                f"min(samples, features) = {mins}.")

        fitter = func_map[func_name](**kwargs)
        adata.obsm[key] = fitter.fit_transform(x_to_use)
        adata.uns[key] = kwargs.copy()
        adata.uns[key]['method'] = method

    return _func


# Assing every method to a global function
for func_name in func_map.keys():
    globals()[func_name] = get_func(func_name)


def _correct_bin_names(bin_names):
    for i in range(len(bin_names)):
        bin_names[i] = bin_names[i].replace(
            ':', '_', bin_names[i].count(':')-1)
    return bin_names


def cl_cisTopic(adata, key, x_to_use, topics=40, iterations=150, **kwargs):
    """
    In Cellar, cisTopic is meant to be used with scATAC-seq data.
    https://www.nature.com/articles/s41592-019-0367-1
    This method uses LDA to infer cis regulatory topics. We use it here
    solely as a "dimensionality reduction" method where the topics
    found can serve as components. Since cisTopic is only available for R,
    we rely on rpy2 to call R functions.

    Parameters
    __________

    adata: anndata.AnnData object
    key: str
        Key to store the reduced data under adata.obsm
    x_to_use: str
        Ignored. Present only for consistency.
    topics: int
        Number of topics to consider. Will run cisTopic for
        topics - 5, topics, and topics + 5 and select the best one.
    iterations: int
        Number of iterations.
    kwargs: dict
        Ignored. Present only for consistency.
    """
    topics = int(topics)
    topic_list = [topics, topics + 5, topics - 5]
    # Unfortunately, with most R functions we cannot use backed mode
    # so we have to load adata into memory. This can potentially lead to
    # memory issues if multiple adatas are found in memory at the same time.
    # Also transpose matrix since cisTopic accepts data in (bin, cell) format.
    mat = adata.to_memory().X.T.copy()

    # If mat is sparse, then we convert mat to an R friendly sparse format
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
    # Set row and column names
    var_names = _correct_bin_names(adata.var_names.to_numpy())
    mat = r("`rownames<-`")(mat, ro.vectors.StrVector(var_names))
    mat = r("`colnames<-`")(mat, ro.vectors.StrVector(adata.obs.index))
    cisTopic = importr('cisTopic')
    logger.info('Creating cisTopic object.')
    cc = cisTopic.createcisTopicObject(mat)

    logger.info('Running LDA Models. This may take a while...')
    cc = cisTopic.runWarpLDAModels(
        cc,
        topic=numpy2ri.py2rpy(topic_list),
        nCores=2,  # Careful with this, since each run duplicates the data
        iterations=iterations,
        addModels=False,
        returnType='selectedModel')

    cellassign = cisTopic.modelMatSelection(cc, 'cell', 'Probability')
    # Transpose back and convert to float32
    cellassign = np.array(cellassign).T.copy().astype(np.float32)

    adata.obsm[key] = cellassign
    adata.uns[key] = {
        'method': 'cisTopic',
        'topics': topics,
        'iterations': iterations
    }
    # Clean mat
    del mat
    gc.collect()


def clear_x_emb_dependends(adata):
    if 'x_emb_2d' in adata.obsm:
        adata.obsm.pop('x_emb_2d')
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')
