import glob
import os
import numpy as np
import scipy.sparse
from anndata._core.sparse_dataset import SparseDataset


def is_sparse(x):
    if isinstance(x, SparseDataset):
        return True
    return scipy.sparse.issparse(x)


def _gene_value_protein_2_symbol(feature_values, adata):
    if 'gene_symbols' in adata.var:
        symbols = feature_values.copy()
        feats_in_idx = np.isin(feature_values, adata.var_names)
        symbols[feats_in_idx] = adata[
            :, feature_values[feats_in_idx]].var['gene_symbols']
    else:
        symbols = feature_values
    return symbols


def get_title_from_feature_list(adata, feature_values):
    symbols = _gene_value_protein_2_symbol(feature_values, adata)
    title = symbols[0]
    for i in range(1, min(3, len(symbols))):
        title += ", " + symbols[i]
    if len(symbols) > 3:
        title += "..."

    return title


def _check_proteins(adata):
    if 'protein.X' not in adata.obsm:
        return None
    if 'protein.var_names' in adata.uns:
        protein_names = np.array(adata.uns['protein.var_names']).astype(str)
    else:
        protein_names = np.arange(adata.obsm['protein.X'].shape[1])
        protein_names = protein_names.astype(str)
        adata.uns['protein.var_names'] = protein_names
    protein_names = np.char.add('Protein_', protein_names)
    return protein_names


def get_server_dataset_dict(root='data'):
    """
    Reads the datasets in the above folder and returns a dict
    {
        'path': 'dataset name', ...
    }
    """
    dataset_dict = {}
    # Split datasets so we can put uploaded first
    server_datasets = sorted(glob.glob(
        os.path.join(root, 'server', '**/*.h5ad'), recursive=True))
    uploaded_datasets = sorted(glob.glob(
        os.path.join(root, 'uploaded', '**/*.h5ad'), recursive=True))

    datasets = uploaded_datasets + server_datasets

    for dataset in datasets:
        paths = dataset.split("/")
        fn = paths[-1][:-5]  # filename, remove extension
        value = ""
        if paths[1] == 'uploaded':
            value += "uploaded / "
        for subdir in paths[2:-1]:
            value += subdir + " / "
        value += fn
        dataset_dict[dataset] = value

    return dataset_dict


if __name__ == "__main__":
    print(get_server_dataset_dict('data'))
