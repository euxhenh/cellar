import diffxpy.api as de
import numpy as np
import gseapy as gp
import json

from controller.cellar.utils.exceptions import InternalError


def ttest(adata, cluster_id, alpha=0.05):
    if 'labels' not in adata.obs:
        raise InternalError("No labels found in adata.")

    grouping = np.zeros(adata.shape[0])
    grouping[adata.obs['labels'].to_numpy() == cluster_id] = 1

    test = de.test.t_test(
        data=adata.X,
        grouping=grouping,
        gene_names=adata.var['gene_symbols'],
        is_logged=True
    )

    test = test.summary(qval_thres=alpha, fc_upper_thres=1)
    test = test.drop(['zero_mean', 'zero_variance'], axis=1)
    test = test.sort_values(by='log2fc', ascending=False)
    test = test.round(3)

    adata.uns['de_genes'] = test['gene'].copy()
    return test


def get_cell_type_dic(marker_file):
    '''could be removed if use .gmt file'''
    f=open(marker_file,'r')
    old_markers = json.loads(f.read())
    f.close()
    tissues = old_markers.keys()
    
    dic={}
    
    for tissue in tissues:
        for cell_type in old_markers[tissue].keys():
            new_type = tissue +' - ' + cell_type
            genes = old_markers[tissue][cell_type]
            GENES = []
            for g in genes:
                GENES.append(g.upper())
            dic[new_type] = GENES
            
    return dic
    


def enrich(adata, gene_set, de_genes_list):
    
    if gene_set=='Cell Types':
        
        '''choose one out of the 2 implementations: '''
        gene_set = 'cell_type_markers.gmt'
        #gene_set = get_cell_type_dic('cell_type_marker.json')
        '''choose one out of the 2 implementations'''
        
        background = 'hsapiens_gene_ensembl'
        
        enr = gp.enrichr(
            gene_list=de_genes_list,
            gene_sets=gene_set,
            background = background,
            outdir=None,
            no_plot=True,
            cutoff=0.05
        )
        
        res = enr.results

        res = res.drop(['Gene_set'], axis=1)
        res = res.round(3)
        res.rename(columns={'P-value': 'pval',
                        'Adjusted P-value': 'qval'}, inplace=True)

        return res

        
    else:
        gene_set = [gene_set]

        enr = gp.enrichr(
            gene_list=de_genes_list,
            gene_sets=gene_set,
            outdir=None,
            no_plot=True,
            cutoff=0.05
        )

        res = enr.results

        res = res.drop(['Gene_set', 'Old P-value', 'Old Adjusted P-value'], axis=1)
        res = res.round(3)
        res.rename(columns={'P-value': 'pval',
                        'Adjusted P-value': 'qval'}, inplace=True)

        return res
