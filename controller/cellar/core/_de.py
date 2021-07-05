import diffxpy.api as de
import numpy as np
import gseapy as gp

from controller.cellar.utils.exceptions import InternalError


def ttest(adata, cluster_id, alpha=0.05):
    if 'labels' not in adata.obs and isinstance(cluster_id, int):
        raise InternalError("No labels found in adata.")

    grouping = np.zeros(adata.shape[0])
    if isinstance(cluster_id, str):
        if 'subsets' not in adata.uns:
            raise InternalError("'subsets' key not found in adata.")
        if cluster_id not in adata.uns['subsets']:
            raise InternalError("subset name not found in adata.")

        grouping[adata.uns['subsets'][cluster_id]] = 1
    else:
        grouping[adata.obs['labels'].to_numpy() == cluster_id] = 1

    if 'gene_symbols' in adata.var:
        gene_names = adata.var['gene_symbols']
    else:
        gene_names = adata.var.index.to_numpy()

    test = de.test.t_test(
        data=adata.X,
        grouping=grouping,
        gene_names=gene_names,
        is_logged=True
    )

    test = test.summary(qval_thres=alpha, fc_upper_thres=1)
    test = test.drop(['zero_mean', 'zero_variance'], axis=1)
    test = test.sort_values(by='log2fc', ascending=False)
    test = test.round(3)

    adata.uns['de_genes'] = test['gene'].copy()
    return test


def enrich(adata, gene_set, de_genes_list):
    background = None

    if gene_set == 'Cell Type':
        gene_set = 'data/cell_type_markers.gmt'
        background = 'hsapiens_gene_ensembl'
    else:
        gene_set = [gene_set]

    enr = gp.enrichr(
        gene_list=de_genes_list,
        gene_sets=gene_set,
        background=background,
        outdir=None,
        no_plot=True,
        cutoff=0.05
    )

    res = enr.results

    res.drop(['Gene_set', 'Old P-value', 'Old Adjusted P-value'],
             axis=1, inplace=True, errors='ignore')
    res = res.round(3)
    res.rename(columns={'P-value': 'pval',
                        'Adjusted P-value': 'qval'}, inplace=True)

    if gene_set == 'data/cell_type_markers.gmt':
        res.sort_values(by=['qval'], inplace=True)

    return res
