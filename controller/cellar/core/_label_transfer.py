import gc
import numpy as np
import pandas as pd
import scanpy as sc
from app import logger
from controller.cellar.utils.exceptions import InternalError

# SingleR stuff
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, r
import rpy2.robjects as ro
import anndata2ri
import scipy.sparse as sp


def _transfer_annotations(ref_labels, ref_anns, labels):
    """
    Transfers annotations by matching labels.
    """
    main_anns = np.array([""] * len(labels), dtype='U200')
    unq_labels, unq_indices = np.unique(ref_labels, return_index=True)

    for unq_label, unq_index in zip(unq_labels, unq_indices):
        main_anns[labels == unq_label] = ref_anns[unq_index]
    return main_anns


def cl_ExactLT(
        adata, key='labels', x_to_use=None, clear_annotations=True, extras={}):
    """
    Exact label transfer based on Cell IDs. Cell IDs between the reference
    and main datasets should match.
    """
    if 'ref' not in extras:
        raise InternalError("No reference adata found.")
    if 'labels' not in extras['ref'].obs:
        raise InternalError("No labels found in reference dataset.")

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    main_df = pd.DataFrame(index=adata.obs_names.to_numpy().astype(str))

    if 'annotations' in extras['ref'].obs:
        annotations = extras['ref'].obs['annotations'].to_numpy()
    else:
        annotations = np.array(
            [""] * extras['ref'].shape[0], dtype='U200')

    ref_df = pd.DataFrame(
        index=extras['ref'].obs_names.to_numpy().astype(str),
        columns=['labels', 'annotations'])
    ref_df['labels'] = extras['ref'].obs['labels'].to_numpy()
    ref_df['annotations'] = annotations

    df = pd.merge(
        left=main_df, right=ref_df, how='left',
        left_index=True, right_index=True)

    df['labels'].fillna(
        np.max(extras['ref'].obs['labels'].to_numpy()) + 1, inplace=True)
    df['annotations'].fillna("No Matching ID", inplace=True)

    adata.obs['labels'] = df['labels'].to_numpy().astype(int)
    adata.obs['annotations'] = df['annotations'].to_numpy().astype('U200')


def cl_Ingest(
        adata, key='labels', x_to_use=None, clear_annotations=True, extras={}):
    """
    Runs Scanpy Ingest to map labels from a reference adata in
    extras['ref'] to adata. Scanpy Ingest is based on projecting adata
    to a PCA space that has been fit on the reference adata. It then
    uses a knn classifier to map embeddings with labels.

    Parameters
    __________

    adata: anndata.AnnData object
    key: str
        Key to use for storing the labels under adata.obs
    x_to_use: str
        Ignored.
    clear_annotations: bool
        Set to True to clear adata.obs['annotations] if present.
    extras: dict
        Must contain a reference adata store in 'ref'.
    """
    if 'ref' not in extras:
        raise InternalError("No reference adata found.")
    if key not in extras['ref'].obs:
        raise InternalError(f"No {key} key was found in the reference data.")

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    # We can only use this method by restricting the datasets to their
    # overlapping genes
    vars_main = adata.var_names.to_numpy().astype('str')
    vars_ref = extras['ref'].var_names.to_numpy().astype('str')
    vars_common = np.intersect1d(vars_main, vars_ref)

    logger.info(f"Found {len(vars_common)} common variables.")
    if len(vars_common) <= 5:
        raise InternalError(
            "Not enough common genes found to run label transfer.")
    # Need to load to memory since we cannot make copies otherwise.
    # This can potentially lead to memory issues.
    adata_main = adata[:, vars_common].to_memory().copy()
    adata_ref = extras['ref'][:, vars_common].to_memory().copy()

    logger.info("Running Scanpy Ingest.")
    # Run Scanpy preprocessing. Ingest won't work if we use Cellar's
    # embeddings.
    sc.pp.pca(adata_main, min(40, len(vars_common) - 1))
    sc.pp.neighbors(adata_main)
    sc.pp.pca(adata_ref, min(40, len(vars_common) - 1))
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)
    sc.tl.ingest(adata_main, adata_ref, obs=key)
    labels = adata_main.obs[key].to_numpy().copy().astype(int)
    # Free memory
    del adata_main, adata_ref
    gc.collect()

    # Also transfer annotations if they are present
    if 'annotations' in extras['ref'].obs:
        adata.obs['annotations'] = _transfer_annotations(
            extras['ref'].obs[key].to_numpy(),
            extras['ref'].obs['annotations'].to_numpy(), labels)
    adata.obs[key] = labels


def cl_SingleR(
        adata, key='labels', x_to_use=None, clear_annotations=True, extras={}):
    """
    Runs SingelR to map labels from a reference adata in
    extras['ref'] to adata. SingleR is a correlation-based method that
    iteratively reduces the number of cell type candidates until one
    cell type is left. Since this runs for every cell separately, this method
    can be very slow.

    Parameters
    __________

    adata: anndata.AnnData object
    key: str
        Key to use for storing the labels under adata.obs
    x_to_use: str
        Ignored.
    clear_annotations: bool
        Set to True to clear adata.obs['annotations] if present.
    extras: dict
        Must contain a reference adata store in 'ref'.
    """
    if 'ref' not in extras:
        raise InternalError("No reference adata found.")
    if key not in extras['ref'].obs:
        raise InternalError("No labels found in reference dataset.")

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    # Restrict adatas to overlapping genes
    vars_main = adata.var_names.to_numpy().astype('str')
    vars_ref = extras['ref'].var_names.to_numpy().astype('str')
    vars_common = np.intersect1d(vars_main, vars_ref)

    logger.info(f"Found {len(vars_common)} common variables.")
    if len(vars_common) <= 5:
        raise InternalError(
            "Not enough common genes found to run label transfer.")
    # Need to load to memory since we cannot make copies otherwise.
    # This can potentially lead to memory issues.
    adata_main = adata[:, vars_common].to_memory().copy()
    adata_ref = extras['ref'][:, vars_common].to_memory().copy()

    def _prep_mat(ad):
        """
        Prepare an R friendly matrix.
        """
        mat = ad.X.T
        if sp.issparse(mat):
            mat = anndata2ri.scipy2ri.py2rpy(sp.csr_matrix(mat))
        else:
            mat = numpy2ri.py2rpy(mat)

        mat = r("`rownames<-`")(mat, ro.vectors.StrVector(ad.var_names))
        mat = r("`colnames<-`")(mat, ro.vectors.StrVector(ad.obs.index))
        return mat

    mat = _prep_mat(adata_main)
    mat2 = _prep_mat(adata_ref)

    labels = ro.vectors.StrVector(
        adata_ref.obs[key].values.astype('str').tolist())
    par = r('BiocParallel::SerialParam')()
    r('BiocParallel::register(BiocParallel::SerialParam())')
    try:
        s = importr('SingleR')
    except:
        raise ImportError("Could not import SingleR.")

    logger.info("Running SingleR.")
    labels = s.SingleR(
        test=mat,
        ref=mat2,
        labels=labels,
        BPPARAM=par)

    labels = list(labels.slots['listData'].rx2("labels"))
    labels = np.asarray(labels).astype(int)
    # Free memory
    del adata_main, adata_ref, mat, mat2
    gc.collect()

    # Also transfer annotations if they are present
    if 'annotations' in extras['ref'].obs:
        adata.obs['annotations'] = _transfer_annotations(
            extras['ref'].obs[key].to_numpy(),
            extras['ref'].obs['annotations'].to_numpy(), labels)
    adata.obs[key] = labels
