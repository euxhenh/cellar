from sklearn.metrics import (calinski_harabasz_score, davies_bouldin_score,
                             silhouette_score)

from controller.cellar.utils.exceptions import InvalidArgument
from ._unit import Unit


class Eval_Silhouette(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    @staticmethod
    def get(x, labels):
        return silhouette_score(x, labels)


class Eval_DaviesBouldin(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    @staticmethod
    def get(x, labels):
        # Return negative the result, because the db score
        # assumes a better clustering if the score is lower
        return -davies_bouldin_score(x, labels)


class Eval_CalinskiHarabasz(Unit):
    def __init__(self, **kwargs):
        self.kwargs = kwargs

    @staticmethod
    def get(x, labels):
        return calinski_harabasz_score(x, labels)


def get_eval_obj(eval_obj='Silhouette'):
    if eval_obj == 'Silhouette':
        return Eval_Silhouette
    elif eval_obj == 'DaviesBouldin':
        return Eval_DaviesBouldin
    elif eval_obj == 'CalinskiHarabasz':
        return Eval_CalinskiHarabasz
    else:
        raise InvalidArgument("Evaluation method not found.")
