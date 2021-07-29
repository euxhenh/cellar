import gc
import numpy as np
import scanpy as sc
from app import logger
from controller.cellar.utils.exceptions import InternalError


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

    logger.info("Running Scanpy Ingest.")

    vars_main = adata.var_names.to_numpy().astype('str')
    vars_ref = extras['ref'].var_names.to_numpy().astype('str')
    vars_common = np.intersect1d(vars_main, vars_ref)

    logger.info(f"Found {len(vars_common)} common variables.")
    if len(vars_common) <= 1:
        raise InternalError(
            "Not enough common genes found to run label transfer.")

    adata_main = adata[:, vars_common].to_memory().copy()
    adata_ref = extras['ref'][:, vars_common].to_memory().copy()

    # Add neighbors
    sc.pp.pca(adata_main)
    sc.pp.pca(adata_ref)
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
