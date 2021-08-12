import gc
import numpy as np
import scanpy as sc
from app import logger
from controller.cellar.utils.exceptions import InternalError

# SingleR stuff
from rpy2.robjects.packages import importr
from rpy2.robjects import numpy2ri, r
import rpy2.robjects as ro
import anndata2ri
import scipy.sparse as sp


def cl_Ingest(
        adata, key='labels', x_to_use='x', clear_annotations=True, extras={}):
    if 'ref' not in extras:
        raise InternalError("No reference adata found.")
    if 'labels' not in extras['ref'].obs:
        raise InternalError("No labels found in reference dataset.")

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    vars_main = adata.var_names.to_numpy().astype('str')
    vars_ref = extras['ref'].var_names.to_numpy().astype('str')
    vars_common = np.intersect1d(vars_main, vars_ref)

    logger.info(f"Found {len(vars_common)} common variables.")
    if len(vars_common) <= 1:
        raise InternalError(
            "Not enough common genes found to run label transfer.")

    adata_main = adata[:, vars_common].to_memory().copy()
    adata_ref = extras['ref'][:, vars_common].to_memory().copy()

    logger.info("Running Scanpy Ingest.")
    sc.pp.pca(adata_main, min(40, len(vars_common) - 1))
    sc.pp.pca(adata_ref, min(40, len(vars_common) - 1))
    sc.pp.neighbors(adata_main)
    sc.pp.neighbors(adata_ref)
    sc.tl.umap(adata_ref)

    sc.tl.ingest(adata_main, adata_ref, obs='labels')

    labels = adata_main.obs['labels'].to_numpy().copy().astype(int)

    del adata_main
    del adata_ref
    gc.collect()

    # Also transfer annotations if they are present
    if 'annotations' in extras['ref'].obs:
        ref_labels = extras['ref'].obs['labels'].to_numpy()
        ref_anns = extras['ref'].obs['annotations'].to_numpy()
        main_anns = np.array([""] * adata.shape[0], dtype='U200')
        unq_labels, unq_indices = np.unique(ref_labels, return_index=True)

        for unq_label, unq_index in zip(unq_labels, unq_indices):
            main_anns[labels == unq_label] = ref_anns[unq_index]

        adata.obs['annotations'] = main_anns

    adata.obs['labels'] = labels


def cl_SingleR(
        adata, key='labels', x_to_use='x', clear_annotations=True, extras={}):
    if 'ref' not in extras:
        raise InternalError("No reference adata found.")
    if 'labels' not in extras['ref'].obs:
        raise InternalError("No labels found in reference dataset.")

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    vars_main = adata.var_names.to_numpy().astype('str')
    vars_ref = extras['ref'].var_names.to_numpy().astype('str')
    vars_common = np.intersect1d(vars_main, vars_ref)

    logger.info(f"Found {len(vars_common)} common variables.")
    if len(vars_common) <= 1:
        raise InternalError(
            "Not enough common genes found to run label transfer.")

    adata_main = adata[:, vars_common].to_memory().copy()
    adata_ref = extras['ref'][:, vars_common].to_memory().copy()

    mat = adata_main.X.T.copy()

    if sp.issparse(mat):
        mat = anndata2ri.scipy2ri.py2rpy(sp.csr_matrix(mat))
    else:
        mat = numpy2ri.py2rpy(mat)

    mat = r("`rownames<-`")(
        mat, ro.vectors.StrVector(adata_main.var_names))
    mat = r("`colnames<-`")(mat, ro.vectors.StrVector(adata_main.obs.index))

    mat2 = adata_ref.X.T.copy()
    if sp.issparse(mat2):
        mat2 = anndata2ri.scipy2ri.py2rpy(sp.csr_matrix(mat2))
    else:
        mat2 = numpy2ri.py2rpy(mat2)
    mat2 = r("`rownames<-`")(
        mat2, ro.vectors.StrVector(adata_ref.var_names))
    mat2 = r("`colnames<-`")(mat2, ro.vectors.StrVector(adata_ref.obs.index))

    labels = ro.vectors.StrVector(
        adata_ref.obs['labels'].values.astype('str').tolist())

    par = r('BiocParallel::SerialParam')()
    r('BiocParallel::register(BiocParallel::SerialParam())')
    s = importr('SingleR')

    logger.info("Running SingleR.")

    labels = s.SingleR(
        test=mat,
        ref=mat2,
        labels=labels,
        BPPARAM=par)

    labels = list(labels.slots['listData'].rx2("labels"))
    labels = np.asarray(labels).astype(int)

    del adata_main
    del adata_ref
    del mat
    del mat2
    gc.collect()

    # Also transfer annotations if they are present
    if 'annotations' in extras['ref'].obs:
        ref_labels = extras['ref'].obs['labels'].to_numpy()
        ref_anns = extras['ref'].obs['annotations'].to_numpy()
        main_anns = np.array([""] * adata.shape[0], dtype='U200')
        unq_labels, unq_indices = np.unique(ref_labels, return_index=True)

        for unq_label, unq_index in zip(unq_labels, unq_indices):
            main_anns[labels == unq_label] = ref_anns[unq_index]

        adata.obs['annotations'] = main_anns

    adata.obs['labels'] = labels
