from sklearn.metrics import (calinski_harabasz_score, davies_bouldin_score,
                             silhouette_score)

from ._unit import Unit


class Eval_Silhouette(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def get(self, x, labels):
        return silhouette_score(x, labels, **self.kwargs)


class Eval_DaviesBouldin(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def get(self, x, labels):
        # Return negative the result, because the db score
        # assumes a better clustering if the score is lower
        return -davies_bouldin_score(x, labels, **self.kwargs)


class Eval_CalinskiHarabasz(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    def get(self, x, labels):
        return calinski_harabasz_score(x, labels)
