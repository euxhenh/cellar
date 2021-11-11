from math import floor

import numpy as np
from sklearn.preprocessing import normalize
from scipy.sparse import csc_matrix
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects.vectors import IntVector

from app import logger
from controller.cellar.utils.colors import interpolate_grayimage
from ..utils.exceptions import InternalError
from ..utils.misc import is_sparse

try:
    STvEA = importr('STvEA')
except:
    raise ImportError("Could not import R packages.")


def cl_STvEA(
        adata, key=None, x_to_use=None,
        clean_codex=True, codex_clean_model="gaussian",
        codex_normalize=False, clean_cite=True,
        cite_clean_model="gaussian", cite_normalize=True,
        num_chunks=8, k_anchor=20, k_filter=100,
        k_score=80, k_weight=100, extras={}):
    """
    Integrates CITE-seq and CODEX data. Uses an anchor correction method
    to map features between the two. Must be used on a CODEX dataset and
    will expand that dataset by adding expresion levels from the genes
    that were mapped from CITE-seq.

    Parameters
    __________
    adata: anndata.AnnData object
        Corresponds to the CODEX dataset.
    key, x_to_use: Ignored. Present for consistency.
    normalize_cite, normalize_codex: bool
        If True, will normalize each data type to the range (0, 1)
        to improve the performance of STvEA.
    extras: dict
        Must have a 'ref' key that corresponds to the reference adata.
        Corresponds to the CITE-seq data.
        The reference adata must have the adata.obsm['proteins'],
        adata.uns['proteins'], and adata.obsm['x_emb'] keys populated.

    Will add gene expression data under adata.obsm['genes'] and
    gene names under adata.uns['genes'].
    """
    if 'ref' not in extras:
        raise InternalError("No reference dataset found.")
    adata_ref = extras['ref']
    if 'proteins' not in adata_ref.obsm:
        raise InternalError("No protein data found in reference data.")
    if 'proteins' not in adata_ref.uns:
        raise InternalError("No protein names found in reference data.")
    if 'x_emb' not in adata_ref.obsm:
        raise InternalError("No embeddings found in reference data.")
    # First find overlapping proteins
    cite_protein_names = np.char.upper(np.array(
        adata_ref.uns['proteins']).astype(str).flatten())
    if len(cite_protein_names) != adata_ref.obsm['proteins'].shape[1]:
        raise InternalError("Number of protein names does not match " +
                            "protein expression matrix width.")
    codex_protein_names = np.char.upper(adata.var_names.to_numpy().astype(str))
    protein_overlap, codex_ind, cite_ind = np.intersect1d(
        codex_protein_names, cite_protein_names, return_indices=True)
    if len(protein_overlap) < 2:
        raise InternalError(
            "Found little overlap between CODEX and CITE-seq proteins. " +
            "Please make sure the two datasets are using the same " +
            "naming convention.")
    # Construct R matrices for codex-protein, cite-protein, and cite-latent
    codex_mat = adata.X[:, codex_ind]
    if is_sparse(codex_mat):
        codex_mat = codex_mat.todense()
    codex_mat = np.array(codex_mat)
    if clean_codex:
        if codex_clean_model == "gaussian":
            # Following https://github.com/CamaraLab/STvEA/blob/master/R/data_processing.R#L76
            codex_mat = codex_mat - np.min(codex_mat)
            avg_cell_total = np.mean(np.sum(codex_mat, axis=1))
            codex_mat = normalize(codex_mat, norm='l1') * avg_cell_total
    codex_mat_r = ro.numpy2ri.py2rpy(codex_mat)
    codex_mat_r = ro.r("`colnames<-`")(
        codex_mat_r, ro.StrVector(protein_overlap))
    codex_mat_r = ro.r("`rownames<-`")(
        codex_mat_r, ro.StrVector(adata.obs_names.to_numpy()))
    if clean_codex:
        logger.info("Cleaning CODEX data.")
        if codex_clean_model == "gaussian":
            codex_mat_r = STvEA.CleanCODEX_gaussian_internal(codex_mat_r)
        else:
            codex_mat_r = STvEA.CleanCODEX_nb_internal(
                codex_mat_r, num_cores=1, normalize=codex_normalize)

    cite_mat = adata_ref.obsm['proteins'][:, cite_ind]
    if is_sparse(cite_mat):
        cite_mat = cite_mat.todense()
    cite_mat = np.array(cite_mat)
    if clean_cite:
        if cite_clean_model == "gaussian":
            # Following https://github.com/CamaraLab/STvEA/blob/master/R/data_processing.R#L140
            cite_mat = normalize(cite_mat, norm='l1')
            cite_mat = np.log(cite_mat * 1e4 + 1)
    cite_mat_r = ro.numpy2ri.py2rpy(cite_mat)
    cite_mat_r = ro.r("`colnames<-`")(
        cite_mat_r, ro.StrVector(protein_overlap))
    cite_mat_r = ro.r("`rownames<-`")(
        cite_mat_r, ro.StrVector(adata_ref.obs_names.to_numpy()))
    if clean_cite:
        logger.info("Cleaning CITE-seq data.")
        if cite_clean_model == "gaussian":
            cite_mat_r = STvEA.CleanCITE_gaussian_internal(
                cite_mat_r, num_cores=1)
        else:
            cite_mat_r = STvEA.CleanCITE_nb_internal(
                cite_mat_r, num_cores=1, factr=1e-9, normalize=cite_normalize)

    cite_emb = np.array(adata_ref.obsm['x_emb'])
    cite_emb_r = ro.numpy2ri.py2rpy(cite_emb)

    num_cc = len(protein_overlap) - 1
    logger.info("Running STvEA mapping.")
    corrected_codex_r = STvEA.MapCODEXtoCITE_internal(
        cite_mat_r, codex_mat_r, cite_emb_r,
        num_chunks=num_chunks, num_cores=1, num_cc=num_cc,
        k_anchor=k_anchor, k_filter=k_filter, k_score=k_score,
        k_weight=k_weight)

    logger.info("Getting transfer matrix.")
    transfer_mat_r = STvEA.GetTransferMatrix_internal(
        cite_mat_r,
        corrected_codex_r,
        max(1, floor(cite_mat.shape[0] * 0.002)),
        0.1
    )
    # No need to subtract 1 since these indices start from 0 for dgCMatrix
    data = np.array(transfer_mat_r.slots['x']).astype(float)
    indices = np.array(transfer_mat_r.slots['i'])
    indptr = np.array(transfer_mat_r.slots['p'])
    transfer_mat = csc_matrix((data, indices, indptr),
                              (codex_mat.shape[0], cite_mat.shape[0]))
    cite_exp_mat = adata_ref.X / np.sum(adata_ref.X, axis=1, keepdims=True)
    cite_exp_mat[np.isnan(cite_exp_mat)] = 0
    codex_exp_mat = transfer_mat.dot(cite_exp_mat)
    adata.obsm['genes'] = codex_exp_mat
    if 'gene_symbols' in adata_ref.var:
        genes = adata_ref.var['gene_symbols']
    else:
        genes = adata_ref.var_names
    adata.uns['genes'] = genes.to_numpy().astype(str)
