import logging

import numpy as np
from joblib import Parallel, delayed
from app import  logger

from ._evaluation import Eval_Silhouette
from ..utils.validation import _validate_n_jobs


def cluster_multiple(x, obj_def, k_list=np.array([2, 4, 8, 16]),
                     attribute_name='n_clusters', eval_obj=None, x_eval=None,
                     method_name='fit_predict', n_jobs=None, **kwargs):
    """
    Runs clustering for multiple n_clusters.

    Parameters
    __________

    x: array, shape (n_samples, n_features)
        The data array.

    obj_def: object name
        Object to be instantiated in this function.

    k_list: array, shape (n_trials,), dtype int, default [2, 4, 8, 16]
        Array containing the different values of clusters to try.

    attribute_name: string, default 'n_clusters'
        Name of the obj.attribute_name that corresponds to n_clusters.

    eval_obj: Eval or None, default None
        Evaluation object to compare performance of different trials.
        If set to None, it will be initialized to SilhouetteScore by default.

    method_name: string, default 'fit_predict'
        Name of the method obj.method_name(x) that returns labels given x.

    n_jobs: int or None, default None
        Number of jobs to use if multithreading. See
        https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.
        If None or 1, will not run multithreading.

    **kwargs: dictionary
        Dictionary of parameters that will get passed to obj_def
        when instantiating it.

    Returns
    _______

    top_y: array, dtype=obj.method_name(x).dtype, shape (n_samples,)
        List of labels that correspond to the best clustering k, as
        evaluated by eval_obj.

    """

    # Default evaluation object
    if eval_obj is None:
        eval_obj = Eval_Silhouette()

    # If n_jobs = -1, run all threads
    n_jobs = _validate_n_jobs(n_jobs)

    if n_jobs == 1:
        top_y, top_score, top_k = None, -np.Inf, -1
        score_list = []

        for i, k in enumerate(k_list):
            kwargs[attribute_name] = k
            # Cluster
            y = getattr(obj_def(**kwargs), method_name)(x)
            # Evaluate
            score = eval_obj.get(x if x_eval is None else x_eval, y)
            score_list.append(score)
            # Update if better (higher) score (saves memory)
            if score > top_score:
                top_y, top_score, top_k = y, score, k

            logger.info(
                "Finished clustering with k={0}. Score={1:.2f}.".format(k, score))
    else:
        logger.info("Running multiple threads with n_jobs={0}.".format(n_jobs))

        kwargs_list = []
        for i, k in enumerate(k_list):
            kcopy = kwargs.copy()
            kcopy[attribute_name] = k
            kwargs_list.append(kcopy)

        # Run clustering in parallel
        y_list = Parallel(n_jobs=n_jobs)(
            delayed(getattr(obj_def(**kwargs_list[i]), method_name))(x)
            for i in range(len(k_list)))
        # Run evaluation in parallel
        score_list = Parallel(n_jobs=n_jobs)(
            delayed(eval_obj.get)(x if x_eval is None else x_eval, y)
            for y in y_list)

        # Find the best score
        top_index = np.argmax(score_list)
        top_y, top_k = y_list[top_index], k_list[top_index]

        # Log scores
        for k, score in zip(k_list, score_list):
            logger.info(
                "Finished clustering with k={0}. Score={1:.2f}.".format(k, score))

    logger.info(
        "Finished clustering. Best score achieved for k={0}.".format(top_k))

    return top_y, score_list
