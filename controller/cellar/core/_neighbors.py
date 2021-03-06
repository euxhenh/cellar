from logging import log
import os

import faiss
import numpy as np
import pandas as pd
from app import logger
from numpy.core.fromnumeric import shape
from scipy.sparse import csr_matrix
from sklearn.neighbors import KDTree, kneighbors_graph
from ..utils.tile_generator import _read_verify_10x_df
from ..utils.exceptions import InternalError, UserError


def full_knn(x, n_neighbors=15, mode='connectivity'):
    """
    Compute the exact nearest neighbors graph using sklearn.
    """
    kg = kneighbors_graph(x, n_neighbors=n_neighbors, mode=mode)
    sources, targets = kg.nonzero()
    weights = np.array(kg[sources, targets])[0]

    return sources, targets, weights


def faiss_knn(x, n_neighbors=15):
    """
    Compute approximate neighbors using the faiss library. We are using
    the CPU based release of faiss.
    """
    n_samples = x.shape[0]
    n_features = x.shape[1]
    x = np.ascontiguousarray(x)  # important

    index = faiss.IndexHNSWFlat(n_features, min(15, n_samples))
    index.add(x)

    weights, targets = index.search(x, n_neighbors)

    sources = np.repeat(np.arange(n_samples), n_neighbors)
    targets = targets.flatten()
    # weights = weights.flatten()

    if -1 in targets:
        raise InternalError("Not enough neighbors were found. Please consider "
                            "reducing the number of neighbors.")
    # TODO fix weights
    return sources, targets, np.full(len(sources), 1)


def knn_auto(x, n_neighbors=15, mode='connectivity', method='auto'):
    if method == 'auto':
        if x.shape[0] > 5000:
            print("Dataset is too large. Finding approximate neighbors.")
            return faiss_knn(x, n_neighbors=n_neighbors)
        return full_knn(x, n_neighbors=n_neighbors, mode=mode)
    elif method == 'full':
        return full_knn(x, n_neighbors=n_neighbors, mode=mode)

    return faiss_knn(x, n_neighbors=n_neighbors)


def _get_coordinates(rid, n_samples, x, y):
    """
    Matches the IDs in rid and adata.obs.index and returns the appropriate
    x and y coordinates. It will raise an error if adata contains
    any IDs not present in rid.
    """
    if rid.shape[0] < n_samples:
        logger.warn("Found more samples than coordinates. Ignoring samples...")
    elif rid.shape[0] > n_samples:
        logger.warn("Found more coordinates than samples. Truncating...")

    indices = np.arange(n_samples)
    overlap, ind1, ind2 = np.intersect1d(rid, indices, return_indices=True)
    if overlap.size == 0:
        raise UserError("No rid values correspond to samples.")
    logger.info(f"Using {overlap.shape[0]} coordinates.")

    return rid[ind1], x[ind1], y[ind1], ind2


def _subsample_from_overlap(rid, n_samples, x, y, subsample_n=5000):
    """
    Looks at the overlap between arange(len(adata)) and rid and
    samples subsample_n points.
    """
    # NOTE: rid maps to the ordinal index in adata and NOT adata.obs.index
    common, ind1, ind2 = np.intersect1d(
        rid, np.arange(n_samples), return_indices=True)
    if common.shape[0] <= subsample_n:
        ind_of_ind = np.arange(len(common))
    else:
        ind_of_ind = np.random.choice(
            np.arange(len(common)), subsample_n, replace=False)
    df_idx = ind1[ind_of_ind]
    rid, x, y = rid[df_idx], x[df_idx], y[df_idx]
    sample_idx = ind2[ind_of_ind]

    return rid, x, y, sample_idx


def get_spatial_knn_graph(
        path_to_df, adata, n_neighbors=3, key='spatial_nneigh',
        subsample=False, subsample_n=5000):
    """
    Constructs a nearest neighbors graph using CODEX spatial tiles
    and return a sparse adjacency matrix.

    Parameters
    __________
    path_to_df: string
        Path to data.csv as generated by cytokit.
        It must contain the columns 'rx', 'ry', 'rid'
        where 'rx' and 'ry' correspond to the x and y coordinates,
        while 'rid' contains the order of the cells in adata.
    n_neighbors: int
        Number of neighbors to compute.
    adata: anndata.AnnData object
    key: str
        If adata is not None, will use this key to store the adjacency
        matrix in adata.obsp
    subsample: bool
        If set to True, will sample points instead of using the full data.

    Returns
    _______
    sparse.csr_matrix symmetric, binary adjacency matrix.
    If adata is not None, will also add the adjacency matrix to adata.obsp,
    indices if subsample is set to True
    """
    if 'x' in adata.obs and 'y' in adata.obs:
        logger.info("Reading x and y coordinates from adata.")
        x = adata.obs['x'].to_numpy().astype(float)
        y = adata.obs['y'].to_numpy().astype(float)
        rid = np.arange(adata.shape[0])
    elif not os.path.exists(path_to_df):
        raise UserError("No data.csv file containing spatial info was found.")
    else:
        # Read dataframe
        data = pd.read_csv(path_to_df)
        # We will only be using x and y coordinates
        if 'rx' not in data or 'ry' not in data:
            raise UserError("data.csv contains incorrect columns. " +
                            "One of 'rx' or 'ry' not found.")
        x = data['rx'].to_numpy().astype(float)
        y = data['ry'].to_numpy().astype(float)
        rid = data['rid'].to_numpy().astype(int)

    if subsample:
        rid, x, y, sample_idx = _subsample_from_overlap(
            rid, adata.shape[0], x, y, subsample_n)
    else:
        # Check shapes and get overlapping IDs
        rid, x, y, sample_idx = _get_coordinates(rid, adata.shape[0], x, y)

    ######### Neighbors ##########
    # Stack into a single data array
    coords = np.array([x, y]).T
    # Use a KD-Tree for fast neighbors computation
    kdt = KDTree(coords, leaf_size=5)
    if n_neighbors + 1 >= coords.shape[0]:
        logger.warn(f"Not enough samples found for f{n_neighbors}." +
                    f"Switching to f{coords.shape[0] // 2 + 1} neighbors.")
        n_neighbors = coords.shape[0] // 2 + 1
    # Add one neighbor since we do not want to include self
    nn_indices = kdt.query(coords, return_distance=False, k=n_neighbors + 1)
    nn_indices = nn_indices[:, 1:]  # Remove self
    ######### Neighbors ##########

    # Construct sparse matrix from nn_indices
    n, d = nn_indices.shape
    x_cord = np.repeat(sample_idx, d)
    adj = csr_matrix(
        (np.full(n * d, 1), (x_cord, sample_idx[nn_indices.flat])),
        shape=(adata.shape[0], adata.shape[0]))
    adj = ((adj + adj.transpose()) > 0).astype(float)

    if not subsample:
        adata.obsp[key] = adj
        adata.uns[key] = {
            'n_neighbors': n_neighbors
        }

    return adj


def get_spatial_knn_graph_10x(
        path_to_df, adata, n_neighbors=3, key='spatial_nneigh',
        subsample=False, subsample_n=5000):
    """
    Constructs a nearest neighbors graph using 10x spatial tiles
    and return a sparse adjacency matrix.

    Parameters
    __________
    path_to_df: string
        Path to data.csv, the file with spatial coordinates. If None, check
        whether 'sptial_dict' is already stored in adata.uns
    n_neighbors: int
        Number of neighbors to compute.
    adata: anndata.AnnData object
    key: str
        If adata is not None, will use this key to store the adjacency
        matrix in adata.obsp
    subsample: bool
        If set to True, will sample points instead of using the full data.

    Returns
    _______
    sparse.csr_matrix symmetric, binary adjacency matrix.
    If adata is not None, will also add the adjacency matrix to adata.obsp,
    indices if subsample is set to True
    """

    if 'spatial_dict' in adata.uns:
        spatial_dict = adata.uns['spatial_dict']
    else:
        if not os.path.exists(path_to_df):
            raise UserError("No data.csv file containing spatial info was found.")
        spatial_dict = _read_verify_10x_df(path_to_df,in_tissue=False)


    ######### filter dic##########
    # some points are not in adata
    topop=[]
    for k in spatial_dict.keys():
        if k not in adata.obs.index:
            topop.append(k)
    for k in topop:
        spatial_dict.pop(k)

    ######### Neighbors ##########

    coords = np.array(list(spatial_dict.values())).astype('float')
    # Use a KD-Tree for fast neighbors computation
    kdt = KDTree(coords, leaf_size=5)
    if n_neighbors + 1 >= coords.shape[0]:
        logger.warn(f"Not enough samples found for f{n_neighbors}." +
                    f"Switching to f{coords.shape[0] // 2 + 1} neighbors.")
        n_neighbors = coords.shape[0] // 2 + 1
    # Add one neighbor since we do not want to include self
    nn_indices = kdt.query(coords, return_distance=False, k=n_neighbors + 1)
    nn_indices = nn_indices[:, 1:]  # Remove self
    ######### Neighbors ##########

    # map index in spatial_dict to adata
    dic_to_a = []
    for k in spatial_dict.keys():
        aidx = list(adata.obs.index).index(k)
        dic_to_a.append(aidx)

    # Construct sparse matrix from nn_indices
    n, d = nn_indices.shape

    dic_idx = list(spatial_dict.keys())
    #idx=[]
    idx=np.arange(n)
    nidx=[]
    for aidx in adata.obs.index:
        if aidx in dic_idx:
            #idx.append(dic_idx.index(aidx))
            didx = dic_idx.index(aidx)
            nn_idx_in_dic = nn_indices[didx]
            for i in nn_idx_in_dic:
                nidx_in_adata = dic_to_a[i]
                nidx.append(nidx_in_adata)

    nidx=np.array(nidx)
    idx=np.repeat(idx,d)

    adj = csr_matrix(
        (np.full(n * d, 1), (idx, nidx.flatten())),
        shape=(adata.shape[0], adata.shape[0]))
    adj = ((adj + adj.transpose()) > 0).astype(float)

    if not subsample:
        adata.obsp[key] = adj
        adata.uns[key] = {
            'n_neighbors': n_neighbors
        }

    return adj
