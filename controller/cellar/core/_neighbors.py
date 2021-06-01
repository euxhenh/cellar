import numpy as np
import faiss
from sklearn.neighbors import kneighbors_graph

from ..utils.exceptions import InternalError


def full_knn(x, n_neighbors=15, mode='connectivity'):
    kg = kneighbors_graph(x, n_neighbors=n_neighbors, mode=mode)
    sources, targets = kg.nonzero()
    weights = np.array(kg[sources, targets])[0]

    return sources, targets, weights


def faiss_knn(x, n_neighbors=15):
    n_samples = x.shape[0]
    n_features = x.shape[1]
    x = np.ascontiguousarray(x)

    index = faiss.IndexHNSWFlat(n_features, 15)
    index.add(x)

    D, I = index.search(x, n_neighbors)

    sources = np.repeat(np.arange(n_samples), n_neighbors)
    targets = I.flatten()
    weights = D.flatten()

    if -1 in targets:
        raise InternalError("Not enough neighbors were found. Please consider "
                            "reducing the number of neighbors.")
    return sources, targets, weights


def knn_auto(x, n_neighbors=15, mode='distance'):
    if x.shape[0] > 5000:
        print("Dataset is too large. Finding approximate neighbors.")
        return faiss_knn(x, n_neighbors=n_neighbors)
    return full_knn(x, n_neighbors=n_neighbors, mode=mode)
