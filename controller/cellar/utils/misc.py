import glob
import os
import scipy.sparse
from anndata._core.sparse_dataset import SparseDataset


def is_sparse(x):
    if isinstance(x, SparseDataset):
        return True
    return scipy.sparse.issparse(x)


def get_server_dataset_dict(root='data/server'):
    """
    Reads the datasets in the above folder and returns a dict
    {
        'path': 'dataset name'
    }

    root dir is assumed to be split into centers/tissues/datasets.h5ad.
    """

    dataset_dict = {}
    datasets = glob.glob(os.path.join(root, '**/*.h5ad'), recursive=True)

    for dataset in datasets:
        paths = dataset.split("/")
        fn = paths[-1][:-5]
        value = ""
        for subdir in paths[2:-1]:
            value += subdir + " / "
        value += fn
        # value = paths[-3] + " / " + paths[-2] + " / " + paths[-1][:-5]
        dataset_dict[dataset] = value

    return dataset_dict


if __name__ == "__main__":
    print(get_server_dataset_dict(
        'data/server'))
