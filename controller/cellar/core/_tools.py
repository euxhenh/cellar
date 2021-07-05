import anndata
import numpy as np
import pandas as pd
from app import logger
from controller.cellar.utils.exceptions import InternalError
from controller.cellar.utils.exceptions import IncorrectFileFormat
from controller.cellar.utils.misc import is_sparse


def read_adata(path, mode='r'):
    try:
        return anndata.read_h5ad(path, mode)
    except Exception as e:
        logger.info(str(e))
        raise IncorrectFileFormat


def cl_get_expression(adata, var_names, op='min'):
    """
    Given an AnnData object and a list of var_names, return the
    expression (if single var) or co-expression (if more than one var)
    value of those features.

    Co-expression is calculated by first scaling every feature column
    to the (0, 1) range and then applying the 'op' operation across
    all participating features for every sample point.
    """
    var_names = np.array(var_names, dtype=str).flatten()

    if len(var_names) == 1:  # Don't normalize if single feature
        x = adata[:, var_names[0]].X

        if is_sparse(x):
            x = np.asarray(x.todense())

        return x.flatten()

    x = adata[:, var_names].X

    if is_sparse(x):
        x = np.asarray(x.todense())

    rang = np.ptp(x, axis=0)  # range of values

    x = (x - x.min(axis=0)) / (rang + np.finfo(float).eps)
    x = np.clip(x, 0, 1)  # just in case there is rounding error

    if op == 'min':
        return np.min(x, axis=1)
    elif op == 'sum':
        return np.sum(x, axis=1)
    else:
        raise InternalError(f"No operation {op} has been implemented.")


hgnc_mart = pd.read_csv('data/HGNC_Mart.csv')


def cl_add_gene_symbol(adata, spliton='.'):
    adata.var_names_make_unique(join='-')

    var_names_trimmed = adata.var_names.to_numpy().astype(str)
    var_names_trimmed = [v.split(spliton)[0] for v in var_names_trimmed]

    # determine format
    if var_names_trimmed[0][:4] == 'ENSG':
        ens_sym = pd.Series(
            hgnc_mart['Approved symbol'].values,
            index=hgnc_mart['Ensembl gene ID']).to_dict()

        gene_symbols = [ens_sym.get(i, i) for i in var_names_trimmed]
    else:
        gene_symbols = var_names_trimmed

    adata.var['gene_symbols'] = np.array(gene_symbols, dtype='U200')
