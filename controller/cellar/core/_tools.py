import anndata
import numpy as np
from numpy.core.fromnumeric import var
import pandas as pd
from app import logger
from pyensembl import EnsemblRelease
from controller.cellar.utils.exceptions import InternalError, InvalidArgument
from controller.cellar.utils.exceptions import IncorrectFileFormat
from controller.cellar.utils.misc import is_sparse


def read_adata(path, mode='r'):
    try:
        return anndata.read_h5ad(path, mode)
    except Exception as e:
        logger.info(str(e))
        raise IncorrectFileFormat


def _collect_x_from_vars(var_names, adata):
    xvar = adata[:, var_names].X if len(var_names) > 0 else None
    if is_sparse(xvar):
        xvar = np.array(xvar.todense())
    return xvar


def _collect_x_from_other(other_names, adata):
    prefix = other_names[0].split(':')[0]
    other_names = [i.split(':')[1] for i in other_names]
    if prefix not in adata.obsm:
        raise InternalError(f"No data with prefix {prefix} found in obsm.")

    all_other_names = np.array(adata.uns[prefix]).astype('U200')
    indices = []
    for oname in other_names:
        indices.append(np.where(all_other_names == oname)[0][0])
    indices = np.array(indices)

    xvar = adata.obsm[prefix][:, indices]
    if is_sparse(xvar):
        xvar = np.array(xvar.todense())
    return xvar


def cl_get_expression(adata, var_names=None, other_names=None, op='min'):
    """
    Given an AnnData object and a list of var_names, return the
    expression (if single var) or co-expression (if more than one var)
    value of those features. Make sure to check if any of the features
    corresponds to a protein name in case of CITE-seq data.

    Co-expression is calculated by first scaling every feature column
    to the (0, 1) range and then applying the 'op' operation across
    all participating features for every sample point.
    """
    if var_names is not None:
        var_names = np.array(var_names, dtype=str).flatten()
        x = _collect_x_from_vars(var_names, adata)
    else:
        other_names = np.array(other_names, dtype=str).flatten()
        x = _collect_x_from_other(other_names, adata)
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


def cl_add_gene_symbol(adata, spliton='.'):
    adata.var_names_make_unique(join='-')
    var_names_trimmed = adata.var_names.to_numpy().astype(str)
    var_names_trimmed = [v.split(spliton)[0] for v in var_names_trimmed]
    var_names_trimmed = [i.upper() for i in var_names_trimmed]
    # determine format
    if var_names_trimmed[0][:4] == 'ENSG':
        data = EnsemblRelease(104)
        gene_symbols = []
        for i in var_names_trimmed:
            try:
                gene_name = data.gene_name_of_gene_id(i)
                gene_symbols.append(gene_name)
            except:
                gene_symbols.append(i)
    elif var_names_trimmed[0][:4] == 'ENSM':
        data = EnsemblRelease(104, species='mouse')
        gene_symbols = []
        for i in var_names_trimmed:
            try:
                gene_name = data.gene_name_of_gene_id(i)
                gene_symbols.append(gene_name)
            except ValueError:
                gene_symbols.append(i)
    else:
        gene_symbols = var_names_trimmed

    adata.var['gene_symbols'] = np.char.upper(gene_symbols)
