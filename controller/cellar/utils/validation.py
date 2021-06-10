from ast import literal_eval

import numpy as np

from .exceptions import InvalidArgument
from .exceptions import MethodNotImplementedError

from anndata import AnnData
#from kneed import KneeLocator




def _validate_n_clusters(n_clusters, h):
    if n_clusters >= h:
        raise InappropriateArgument("Number of clusters needs to be less than "
                                    "the number of samples.")
    if n_clusters < 2:
        raise InappropriateArgument("Number of clusters needs to be greater "
                                    "than 1.")

def _validate_clu_n_clusters(clu_n_clusters, h):
    if isinstance(clu_n_clusters, str):
        try:
            clu_n_clusters = literal_eval(clu_n_clusters)
        except:
            raise InvalidArgument(
                "Incorrect format for the number of clusters.")

    if isinstance(clu_n_clusters, int):
        _validate_n_clusters(clu_n_clusters, h)
        return clu_n_clusters
    elif isinstance(clu_n_clusters, tuple):
        try:
            clu_n_clusters = list(range(*clu_n_clusters))
        except:
            raise InvalidArgument("Incorrect tuple specified for the number of "
                                  "clusters.")
        if len(clu_n_clusters) < 1:
            raise InappropriateArgument(
                "Empty tuple encountered for the number of "
                "clusters.")
    elif isinstance(clu_n_clusters, (list, np.ndarray)):
        if len(clu_n_clusters) < 1:
            raise InappropriateArgument(
                "Empty list encountered for the number of "
                "clusters.")
    else:
        raise InvalidArgument("Incorrect format for the number of clusters.")

    _validate_n_clusters(clu_n_clusters[0], h)
    _validate_n_clusters(clu_n_clusters[-1], h)

    return np.sort(clu_n_clusters)


def _validate_n_jobs(n_jobs):
    if isinstance(n_jobs, str):
        try:
            n_jobs = literal_eval(n_jobs)
        except:
            raise InvalidArgument("Incorrect number of jobs specified.")

    if isinstance(n_jobs, float):
        n_jobs = int(n_jobs)

    if isinstance(n_jobs, int):
        if n_jobs < -1:
            raise InvalidArgument("Incorrect number of jobs specified.")
        elif n_jobs > 8:
            raise InappropriateArgument("Number of jobs is too high.")
        elif n_jobs == 0:
            raise InappropriateArgument("Number of jobs is 0.")
    elif n_jobs is None:
        n_jobs = 1
    elif n_jobs is not None:
        raise InvalidArgument("Incorrect number of jobs specified.")

    return n_jobs

def _validate_ensemble_methods(ensemble_methods):
    methods = [
        "All",
        "KMeans",
        "KMedoids",
        "Spectral",
        "Agglomerative",
        # "DBSCAN",
        # "Birch",
        #"GaussianMixture",
        "Leiden"
    ]

    if isinstance(ensemble_methods, (list, np.ndarray)):
        for method in ensemble_methods:
            if method not in methods:
                raise MethodNotImplementedError(
                    "Incorrect method encountered in ensemble clustering.")
            if method == 'All':
                return methods[1:]
    elif ensemble_methods is None:
        return "default"
    elif isinstance(ensemble_methods, str):
        if ensemble_methods not in methods:
            raise MethodNotImplementedError(
                "Incorrect method encountered in ensemble clusering.")
        if ensemble_methods == 'All':
            return methods[1:]
        else:
            return [ensemble_methods]
    else:
        raise InvalidArgument(
            "Incorrect list provided for ensemble clustering.")
    return ensemble_methods



