from itertools import combinations_with_replacement as cwr
import numpy as np
from rpy2.robjects.packages import importr
import rpy2.robjects as ro
import scipy.sparse as sp

import anndata2ri
from anndata._core.sparse_dataset import SparseDataset
from controller.cellar.utils.exceptions import UserError

from ._neighbors import get_spatial_knn_graph
from app import logger


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
    if adata.shape[0] > 5_000:
        logger.warn("Too many samples were found. Subsampling...")
        adj = get_spatial_knn_graph(
            path_to_df, n_neighbors=n_neighbors, adata=adata, key=key,
            subsample=True, subsample_n=5000)
    elif key in adata.obsp and key in adata.uns and\
            adata.uns[key]['n_neighbors'] == n_neighbors:
        adj = adata.obsp[key]
    else:
        adj = get_spatial_knn_graph(
            path_to_df, n_neighbors=n_neighbors, adata=adata, key=key)

    # We will only consider non-zero rows/cols.
    adata = adata[adj.getnnz(1) > 0]
    # We use the fact that adj is symmetric
    adj = adj[adj.getnnz(1) > 0][:, adj.getnnz(0) > 0]
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
        xx_cords, yy_cords = adata.X.nonzero()
        protein_mat = Matrix.sparseMatrix(
            i=ro.IntVector(xx_cords + 1),
            j=ro.IntVector(yy_cords + 1),
            x=adata.X[xx_cords, yy_cords],
            dims=ro.IntVector(np.array([*adata.shape])))
    else:
        protein_mat = ro.numpy2ri.py2rpy(np.array(adata.X))

    colnames = list(range(1, adata.shape[1] + 1))
    protein_mat = ro.r("`colnames<-`")(protein_mat, ro.IntVector(colnames))

    protein_pairs = np.array(list(cwr(colnames, 2)))
    protein_pairs = ro.numpy2ri.py2rpy(np.array(protein_pairs))

    res = STvEA.AdjScoreProteins_internal(adj_mat, protein_mat, protein_pairs)
    res = ro.pandas2ri.rpy2py_dataframe(res)
    # Careful: Subtract 1
    res['f'] -= 1
    res['g'] -= 1
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
    adj_mat = Matrix.sparseMatrix(
        i=ro.IntVector(x_cords + 1),
        j=ro.IntVector(y_cords + 1),
        x=1,
        dims=ro.IntVector(np.array([*adj.shape])))

    labels = ro.IntVector(adata.obs[labels_key])

    res = STvEA.AdjScoreClustersCODEX_internal(adj_mat, labels, num_cores=2)
    res = ro.pandas2ri.rpy2py_dataframe(res)
    return res
