import anndata
import numpy as np
from numpy.core.fromnumeric import var
import pandas as pd
from app import logger
from controller.cellar.utils.exceptions import InternalError, InvalidArgument
from controller.cellar.utils.exceptions import IncorrectFileFormat
from controller.cellar.utils.misc import is_sparse, _check_proteins


def read_adata(path, mode='r'):
    try:
        return anndata.read_h5ad(path, mode)
    except Exception as e:
        logger.info(str(e))
        raise IncorrectFileFormat


def _split_var_protein(var_names, adata):
    in_idx = np.isin(var_names, adata.var_names)
    vars = var_names[in_idx]
    proteins = var_names[~in_idx]
    return vars, proteins


def _collect_var_protein(var_names, adata):
    vars, proteins = _split_var_protein(var_names, adata)
    all_proteins = _check_proteins(adata)
    xvar = adata[:, vars].X if len(vars) > 0 else None
    xprot = adata.obsm['protein.X'][:, np.isin(
        all_proteins, proteins)] if len(proteins) > 0 else None
    if xvar is not None and is_sparse(xvar):
        xvar = np.array(xvar.todense())
    if xprot is not None and is_sparse(xprot):
        xprot = np.array(xprot.todense())
    if xvar is not None and xprot is not None:
        x = np.hstack([xvar, xprot])
    elif xvar is None:
        x = xprot
    else:
        x = xvar
    return x


def cl_get_expression(adata, var_names, op='min'):
    """
    Given an AnnData object and a list of var_names, return the
    expression (if single var) or co-expression (if more than one var)
    value of those features. Make sure to check if any of the features
    corresponds to a protein name in case of CITE-seq data.

    Co-expression is calculated by first scaling every feature column
    to the (0, 1) range and then applying the 'op' operation across
    all participating features for every sample point.
    """
    var_names = np.array(var_names, dtype=str).flatten()
    x = _collect_var_protein(var_names, adata)
    if x.shape[1] == 1:  # single feature
        return x.flatten()
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
    var_names_trimmed = [i.upper() for i in var_names_trimmed]

    # determine format
    if var_names_trimmed[0][:4] == 'ENSG':
        ens_sym = pd.Series(
            hgnc_mart['Approved symbol'].values,
            index=hgnc_mart['Ensembl gene ID']).to_dict()

        gene_symbols = [ens_sym.get(i, i) for i in var_names_trimmed]
    else:
        gene_symbols = var_names_trimmed

    adata.var['gene_symbols'] = np.array(gene_symbols, dtype='U200')
