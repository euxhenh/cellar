import glob
import os
import scipy.sparse
from anndata._core.sparse_dataset import SparseDataset


def is_sparse(x):
    if isinstance(x, SparseDataset):
        return True
    return scipy.sparse.issparse(x)


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
