from itertools import combinations
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import scipy.sparse as sp

import anndata2ri
from anndata._core.sparse_dataset import SparseDataset
from controller.cellar.utils.exceptions import UserError

from ._neighbors import get_spatial_knn_graph


def adjScoreProteinsCODEX(
        adata, path_to_df, n_neighbors=5, key='spatial_nneigh'):
    """
    Computes a neighbors graph using the spatial tile
    and use it to compute a colocalization score of the proteins.
    See https://github.com/CamaraLab/STvEA

    Parameters
    __________
    adata: anndata.AnnData object
    path_to_df: str
        The path to data.csv file containing spatial information.
        In particular, it must contain the columns 'rx', 'ry', 'rid' (optional)
        where 'rx' and 'ry' correspond to the x and y coordinates,
        while 'rid' contains the order of the cells in adata.
    n_neighbors: int
        Number of neighbors to compute.
    key: str
        If adata is not None, will use this key to store the adjacency
        matrix in adata.obsp
    """
    if key in adata.obsp and key in adata.uns and\
            adata.uns[key]['n_neighbors'] == n_neighbors:
        adj = adata.obsp[key]
    else:
        adj = get_spatial_knn_graph(
            path_to_df, n_neighbors=n_neighbors, adata=adata, key=key)

    try:
        STvEA = importr('STvEA')
        Matrix = importr('Matrix')
    except:
        raise ImportError("Could not import R Libraries.")

    x_cords, y_cords = adj.nonzero()
    # Add one to coordinates since R is retarded
    adj_mat = Matrix.sparseMatrix(
        i=ro.IntVector(x_cords + 1),
        j=ro.IntVector(y_cords + 1),
        x=1,
        dims=ro.IntVector(np.array([*adj.shape])))

    if isinstance(adata.X, SparseDataset) or sp.issparse(adata.X):
        protein_mat = anndata2ri.scipy2ri.py2rpy(sp.csr_matrix(adata.X))
    else:
        protein_mat = ro.numpy2ri.py2rpy(np.array(adata.X))

    protein_names = list(adata.var.index.to_numpy())
    protein_mat = ro.r("`colnames<-`")(
        protein_mat, ro.IntVector(list(range(adata.shape[1]))))

    protein_pairs = list(combinations(list(range(adata.shape[1])), 2))
    protein_pairs = ro.r.matrix(
        [ro.IntVector([i[0], i[1]]) for i in protein_pairs])

    res = STvEA.AdjScoreProteins_internal(
        adj_mat, protein_mat, protein_pairs, num_cores=2)
    res = ro.pandas2ri.rpy2py_dataframe(res)
    return res


def adjScoreClustersCODEX(
        adata, path_to_df, n_neighbors=3, key='spatial_nneigh',
        labels_key='labels'):
    """
    Computes a neighbors graph using the spatial tile
    and use it to compute a colocalization score of the clusters.
    See https://github.com/CamaraLab/STvEA

    Parameters
    __________
    adata: anndata.AnnData object
    path_to_df: str
        The path to data.csv file containing spatial information.
        In particular, it must contain the columns 'rx', 'ry', 'rid' (optional)
        where 'rx' and 'ry' correspond to the x and y coordinates,
        while 'rid' contains the order of the cells in adata.
    n_neighbors: int
        Number of neighbors to compute.
    key: str
        If adata is not None, will use this key to store the adjacency
        matrix in adata.obsp
    """
    if labels_key not in adata.obs:
        raise UserError("No labels found in adata. Cannot compute " +
                        "colocalization score.")

    if key in adata.obsp and key in adata.uns and\
            adata.uns[key]['n_neighbors'] == n_neighbors:
        adj = adata.obsp[key]
    else:
        adj = get_spatial_knn_graph(
            path_to_df, n_neighbors=n_neighbors, adata=adata, key=key)
    try:
        STvEA = importr('STvEA')
        Matrix = importr('Matrix')
    except:
        raise ImportError("Could not import R Libraries.")

    x_cords, y_cords = adj.nonzero()
    # Add one to coordinates since R is retarded
    r_mat = Matrix.sparseMatrix(
        i=ro.IntVector(x_cords + 1),
        j=ro.IntVector(y_cords + 1),
        x=1,
        dims=ro.IntVector(np.array([*adj.shape])))

    labels = ro.IntVector(adata.obs[labels_key])

    res = STvEA.AdjScoreClustersCODEX_internal(r_mat, labels, num_cores=2)
    res = ro.pandas2ri.rpy2py_dataframe(res)
    return res
