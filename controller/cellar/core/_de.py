import diffxpy.api as de
import numpy as np
import gseapy as gp


from controller.cellar.utils.exceptions import InternalError, UserError


def ttest(adata, cluster_id, cluster_id2,
          alpha=0.05, fc_thresh=1, is_logged=True):
    if cluster_id == cluster_id2:
        raise UserError("Selected subsets cannot be the same.")

    if 'labels' not in adata.obs and isinstance(cluster_id, int):
        raise InternalError("No labels found in adata.")

    try:
        alpha = float(alpha)
    except:
        raise UserError("Incorrect alpha value specified.")

    if isinstance(cluster_id, str):
        if 'subsets' not in adata.uns:
            raise InternalError("'subsets' key not found in adata.")
        if cluster_id not in adata.uns['subsets']:
            raise InternalError("Subset name not found in adata.")

        indices1 = adata.uns['subsets'][cluster_id]
    else:
        indices1 = np.where(adata.obs['labels'].to_numpy() == cluster_id)[0]

    if isinstance(cluster_id2, str) and cluster_id2 != 'rest':
        if 'subsets' not in adata.uns:
            raise InternalError("'subsets' key not found in adata.")
        if cluster_id2 not in adata.uns['subsets']:
            raise InternalError("Subset2 name not found in adata.")

        indices2 = adata.uns['subsets'][cluster_id2]
    elif cluster_id2 == 'rest':
        indices2 = None
    else:
        indices2 = np.where(adata.obs['labels'].to_numpy() == cluster_id2)[0]

    if indices2 is None:
        grouping = np.zeros(adata.shape[0])
        grouping[indices1] = 1
        data = adata.X
    else:
        indices = np.concatenate([indices1, indices2])
        grouping = np.zeros(len(indices))
        grouping[:len(indices1)] = 1  # make 1 indices of subset 1, others 0
        # this will sort the data, subset 1 coming first
        data = adata[indices].X

    if 'gene_symbols' in adata.var:
        gene_names = adata.var['gene_symbols']
    else:
        gene_names = adata.var.index.to_numpy()

    test = de.test.t_test(
        data=data,
        grouping=grouping,
        gene_names=gene_names,
        is_logged=is_logged
    )
    test = test.summary(qval_thres=alpha, fc_upper_thres=fc_thresh)
    means = np.array(adata.X[indices1].mean(axis=0)).flatten()

    if indices2 is None:
        indices = np.arange(adata.X.shape[0])
        other_indices = np.delete(indices, indices1)
        other_means = np.array(adata.X[other_indices].mean(axis=0)).flatten()
    else:
        other_means = np.array(adata.X[indices2].mean(axis=0)).flatten()

    test['mean_set1'] = means[test.index.to_numpy()].astype(float)
    test['mean_set2'] = other_means[test.index.to_numpy()].astype(float)
    test = test.drop(['zero_mean', 'zero_variance', 'mean'], axis=1)
    test = test.sort_values(by='log2fc', ascending=False)
    test = test.round(3)

    # adata.uns['de_genes'] = test['gene'].copy()
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

    if res.empty:
        return res

    if gene_set == 'data/cell_type_markers.gmt':
        res.sort_values(by=['qval'], inplace=True)

    return res
