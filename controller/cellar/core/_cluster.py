import igraph as ig
import leidenalg as la
import numpy as np
from app import logger
from scipy.sparse import csr_matrix, issparse
from sklearn.cluster import AgglomerativeClustering, KMeans, SpectralClustering
from sklearn_extra.cluster import KMedoids

from ..utils.exceptions import UserError
from ..utils.validation import _validate_clu_n_clusters
from ._cluster_multiple import cluster_multiple
from ._evaluation import get_eval_obj
from ._neighbors import faiss_knn, knn_auto


def _has_neighbors(adata, params, key='neighs'):
    """
    Checks if there are any neighbors in adata and return those if new params
    equal existing ones in adata.uns[key].

    Parameters
    __________
    adata: anndata.AnnData object
    params: dict
        dictionary of neighbor parameters
    key: str
        Key to look for neighbors params under adata.uns
    """
    if key in adata.uns:
        old_params = adata.uns[key]
        if isinstance(old_params, dict) and isinstance(params, dict):
            if params == old_params and key in adata.obsp:
                if issparse(adata.obsp[key]):
                    return True
    return False


def cl_Leiden(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        partition_type='RBConfigurationVertexPartition', directed=True,
        graph_method='auto', n_neighbors=15, use_cached_neigh=True,
        use_weights=False, neigh_key='neighs', **kwargs):
    """
    Runs the leiden clustering algorithm:
    https://www.nature.com/articles/s41598-019-41695-z

    Parameters
    __________

    adata: anndata.AnnData object
    key: str
        Name of the key to use for labels under adata.obs.
    x_to_use: str
        Key name in adata.obsm to use as feature vectors. If set to 'x'
        then will use adata.X
    clear_annotations: bool
        Set to True to clear adata.obs['annotations'] before clustering.
    partition_type: str
        Can be one of [
            'RBConfigurationVertexPartition',
            'ModularityVertexPartition',
            'RBERVertexPartition',
            'CPMVertexPartition',
            'SurpriseVertexPartition'
        ]
        See https://leidenalg.readthedocs.io/en/stable/reference.html#mutablevertexpartition
    directed: bool
        Set to True to construct a directed neighbors graph and undirected
        otherwise.
    graph_method: str
        Method to use to construct the connectivity graph.
        Can be one of ['auto', 'full', 'approx']. If set to 'auto',
        will automatically choose between 'full' and 'approx' depending
        on the size of the dataset. 'full' uses an exact KNN computation
        to construct the graph, while 'approx' uses the faiss library
        to compute approximate nearest neighbors.
    n_neighbors: int
        Number of nearest neighbors to find using the method specified above.
    use_cached_neigh: bool
        If True, then will check adata if it already contains neighbors
        and if the parameters used to compute these neighbors match
        the current ones.
    use_weights: bool
        Set to True to construct a weighted graph (with weights being the
        distances between points).
    neigh_key: str
        Key to store neighbors in adata.obsp.
    **kwargs: dict
        Any additional arguments will be passed to la.find_partition
    """
    config = {}
    neigh_config = {}

    if clear_annotations and 'annotations' in adata.obs:
        _ = adata.obs.pop('annotations')
    # Using adata.X to perform clustering can be very slow
    neigh_config['x_to_use'] = x_to_use
    if x_to_use == 'x':
        x_to_use = adata.X
        if x_to_use.shape[1] > 500:
            logger.warn(
                "Data is very high-dimensional. This process can " +
                "be slow and may yield suboptimal results. Running " +
                "dimensionality reduction first is highly recommended.")
    else:
        if x_to_use not in adata.obsm:
            raise UserError("No embeddings found. Please " +
                            "run dimensionality reduction first.")
        x_to_use = adata.obsm[x_to_use]

    neigh_config['n_neighbors'] = n_neighbors
    neigh_config['mode'] = 'connectivity'
    neigh_config['graph_method'] = graph_method
    # Check if adata contains neighbors and the configuration is the same
    if use_cached_neigh and _has_neighbors(adata, neigh_config, neigh_key):
        logger.info("Using cached neighbors.")
        sources, targets = adata.obsp[neigh_key].nonzero()
    else:
        # (Re)Compute nearest neighbors
        if n_neighbors >= x_to_use.shape[0]:
            n_neighbors = x_to_use.shape[0] // 2 + 1
            logger.warn(
                "Number of neighbors is greater than " +
                f"the dataset size. Using {n_neighbors} neighbors instead.")
        sources, targets, weights = knn_auto(
            x_to_use, n_neighbors=n_neighbors,
            mode='connectivity', method=graph_method)
        adata.uns[neigh_key] = neigh_config
        adata.obsp[neigh_key] = csr_matrix(
            (weights, (sources, targets)), shape=(x_to_use.shape[0],) * 2)

    # weights = 1 / weights
    # TODO: construct weighted graph
    # kwargs['weights'] = weights if use_weights else None
    # Construct graph
    gg = ig.Graph(directed=directed)
    config['directed'] = directed
    gg.add_vertices(x_to_use.shape[0])
    # Add edges as list of pairs (unweighted)
    gg.add_edges(list(zip(list(sources), list(targets))))

    # Empty input boxes are parsed as empty strings
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None
    # These two partitions do not use resolution
    if partition_type in [
            'ModularityVertexPartition', 'SurpriseVertexPartition']:
        if 'resolution_parameter' in kwargs:
            _ = kwargs.pop('resolution_parameter')
    if not hasattr(la, partition_type):
        raise UserError(f"No partition type '{partition_type}' found.")
    # Perform the clustering
    part = la.find_partition(gg, getattr(la, partition_type), **kwargs)
    config['partition_type'] = partition_type
    # Merge the two dictionaries to also include kwargs passed to leidenalg
    config = {**config, **kwargs}
    config['method'] = 'Leiden'
    adata.obs[key] = np.array(part.membership, dtype=int)
    adata.uns[key] = config
    logger.info(f"Graph modularity: {gg.modularity(part.membership)}")


def _get_wrapper(x, obj_def, n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj='Silhouette', x_eval=None,
                 n_jobs=None, attribute_name='n_clusters', **kwargs):
    """
    Wrapper function for those classes which specify the number of clusters
    in advance and also have fit_predict implemented. Classes include:
    KMedoids, KMeans, SpectralClustering, AgglomerativeClustering, Birch.

    Parameters
    __________
    x: array, shape (n_samples, n_features)
        The data array.
    obj_def: object name
        Object to be instantiated in this function.
    n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
        Array containing the different values of clusters to try,
        or single int specifying the number of clusters,
        or tuple of the form (a, b, c) which specifies a range
        for (x=a; x<b; x+=c)
    eval_obj: Eval or None, default None
        Evaluation object to compare performance of different trials.
    x_eval: np.ndarray or None
        Array to be used for computing the clustering score. If set to None,
        then will use x.
    n_jobs: int or None, default None
        Number of jobs to use if multithreading. See
        https://joblib.readthedocs.io/en/latest/generated/joblib.Parallel.html.
    attribute_name: string, default 'n_clusters'
        Name of the obj.attribute_name that corresponds to n_clusters.
    **kwargs: dictionary
        Dictionary of parameters that will get passed to obj_def
        when instantiating it.
    """
    # Determine type of n_clusters passed
    k = _validate_clu_n_clusters(n_clusters, x.shape[0])
    eval_obj = get_eval_obj(eval_obj)

    # If n_clusters determined to be single integer
    if isinstance(k, int):
        kwargs[attribute_name] = k
        y = obj_def(**kwargs).fit_predict(x)
        score = eval_obj.get(x if x_eval is None else x_eval, y)
        logger.info(
            "Finished clustering with k={0}. Score={1:.2f}.".format(k,
                                                                    score))
        return y, score
    # If n_clusters determined to be a list of integers
    elif isinstance(k, (list, np.ndarray)):
        return cluster_multiple(
            x, obj_def=obj_def, k_list=k, attribute_name=attribute_name,
            eval_obj=eval_obj, x_eval=x_eval, method_name='fit_predict',
            n_jobs=n_jobs, **kwargs)


def cl_KMeans(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        eval_obj='Silhouette', **kwargs):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    if x_to_use not in adata.obsm:
        raise UserError("No embeddings found. Please " +
                        "run dimensionality reduction first.")

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    x_eval = None
    if 'x_emb_2d' in adata.obsm:
        x_eval = adata.obsm['x_emb_2d']

    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    logger.info("Running KMeans.")

    n_clusters = kwargs.pop('n_clusters')

    y, score = _get_wrapper(x_to_use, obj_def=KMeans, n_clusters=n_clusters,
                            eval_obj=eval_obj, x_eval=x_eval,
                            n_jobs=None, attribute_name='n_clusters', **kwargs)

    adata.obs[key] = y


def cl_KMedoids(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        eval_obj='Silhouette', **kwargs):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    if x_to_use not in adata.obsm:
        raise UserError("No embeddings found. Please " +
                        "run dimensionality reduction first.")

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    x_eval = None
    if 'x_emb_2d' in adata.obsm:
        x_eval = adata.obsm['x_emb_2d']

    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    logger.info("Running KMedoids.")

    n_clusters = kwargs.pop('n_clusters')

    y, score = _get_wrapper(x_to_use, obj_def=KMedoids, n_clusters=n_clusters,
                            eval_obj=eval_obj, x_eval=x_eval,
                            n_jobs=None, attribute_name='n_clusters', **kwargs)

    adata.obs[key] = y


def cl_SpectralClustering(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        eval_obj='Silhouette', **kwargs):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    if x_to_use not in adata.obsm:
        raise UserError("No embeddings found. Please " +
                        "run dimensionality reduction first.")

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    x_eval = None
    if 'x_emb_2d' in adata.obsm:
        x_eval = adata.obsm['x_emb_2d']

    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None
    if 'n_components' in kwargs:
        if kwargs['n_components'] == 0:
            kwargs.pop('n_components')

    sources, targets, weights = knn_auto(
        x_to_use, n_neighbors=kwargs['n_neighbors'], mode='connectivity')
    kwargs.pop('n_neighbors')

    n_samples = x_to_use.shape[0]
    nneighs = csr_matrix((weights, (sources, targets)),
                         shape=(n_samples, n_samples))

    logger.info("Running Spectral Clustering.")

    n_clusters = kwargs.pop('n_clusters')

    kwargs['affinity'] = 'precomputed_nearest_neighbors'

    y, score = _get_wrapper(nneighs, obj_def=SpectralClustering,
                            n_clusters=n_clusters, x_eval=x_eval,
                            eval_obj=eval_obj, n_jobs=None,
                            attribute_name='n_clusters', **kwargs)

    adata.obs[key] = y


def cl_Agglomerative(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        eval_obj='Silhouette', **kwargs):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    if x_to_use not in adata.obsm:
        raise UserError("No embeddings found. Please " +
                        "run dimensionality reduction first.")

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    x_eval = None
    if 'x_emb_2d' in adata.obsm:
        x_eval = adata.obsm['x_emb_2d']

    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    logger.info("Running Agglomerative Clustering.")

    n_clusters = kwargs.pop('n_clusters')

    y, score = _get_wrapper(x_to_use, obj_def=AgglomerativeClustering,
                            n_clusters=n_clusters,
                            eval_obj=eval_obj, x_eval=x_eval,
                            n_jobs=None, attribute_name='n_clusters', **kwargs)

    adata.obs[key] = y


def cl_uncertainty(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        n_neighbors=15, seed=None, method='centers', extras=None):

    #logger.info('uncertainty, uncertainty, uncertainty')

    fm = 1
    fstd = 1
    # an='a1'
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')

    if x_to_use not in adata.obsm:
        raise UserError("No embeddings found. Please " +
                        "run dimensionality reduction first.")

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    sources, targets, weights = faiss_knn(x_to_use, n_neighbors=n_neighbors)
    eps = 1
    if method == 'knn':
        source_labels = np.array(adata.obs['labels'][sources])
        target_labels = np.array(adata.obs['labels'][targets])
        same = np.array(source_labels == target_labels).reshape(
            (-1, n_neighbors)).astype(int)
        numerator = (same*weights.reshape((-1, n_neighbors))).sum(axis=1)
        denominator = ((1-same)*weights.reshape((-1, n_neighbors))).sum(axis=1)
        nonconformity = denominator/(numerator+eps)
        m = nonconformity.mean()
        std = nonconformity.std()
        uncertain = (nonconformity > (m*fm+std*fstd)).astype(int)
        print('total uncertain cells:', uncertain.sum())
        # subset 999 = uncertain
        subset_labels = (1-uncertain)*adata.obs['labels']+uncertain*999
        subsets = np.unique(subset_labels)
        adata.uns['subsets'] = {}
        c = 0
        for i in subsets:
            adata.uns['subsets'][str(i)] = []
        for i in subset_labels:
            adata.uns['subsets'][str(i)].append(c)
            c += 1
        for i in subsets:
            adata.uns['subsets'][str(i)] = np.array(
                adata.uns['subsets'][str(i)])

    elif method == 'centers':
        centers = []
        for i in range(np.max(adata.obs['labels'])+1):
            mask = np.array(adata.obs['labels']) == i
            center = (adata.obsm['x_emb'] *
                      mask.reshape((-1, 1))).sum(axis=0)/mask.sum()
            centers.append(center)
        np.array(centers)
        adata.uns['centers'] = centers
        centers = adata.uns['centers']
        x = adata.obsm['x_emb']
        x = x.reshape((x.shape[0], 1, x.shape[1]))
        dist = (x-centers)**2
        dist = np.sqrt(dist.sum(axis=2))  # dist to centers
        kml = np.argmin(dist, axis=1)
        min_d = dist.min(axis=1)
        dist.partition(1, axis=1)
        min2_d = dist[:, 1]
        uncertainty = min_d
        #margin = min_d - min2_d
        #margin /= dist.sum(axis=1)
        #margin = min_d
        m = uncertainty.mean()
        std = uncertainty.std()

        uncertain = (uncertainty > (m*fm+std*fstd)).astype(int)
        uncertain = np.logical_or(uncertain, kml != adata.obs['labels'])
        uncertain = np.array(uncertain)
        logger.info('total uncertain cells:'+str(uncertain.sum()))
        # subset 999 = uncertain
        subset_labels = (1-uncertain)*adata.obs['labels']+uncertain*-1
        adata.obs['labels'] = subset_labels

    return 901


''' ensemble clustering'''
'''
cluster_dict = {
    "KMeans": cl_KMeans,
    "KMedoids": cl_KMedoids,
    "Spectral": cl_SpectralClustering,
    "Agglomerative": cl_Agglomerative,
    "Leiden": cl_Leiden
}

def clu_wrap(method):
    """
    Args:
        method (string): Method to use in the given step.
    Returns:
        object (Unit): Object of the right type.
    """
    if method not in cluster_dict:
        raise NotImplementedError("{0} method not implemented.".format(method))
    return cluster_dict[method]



class Ens_HyperGraph(Unit):
    """
    Ensemble clustering method based on hyper-graph partitioning.
    See https://github.com/GGiecold/Cluster_Ensembles
    """

    def __init__(self, ensemble_methods=["KMedoids", "KMeans", "Spectral"],
                 n_clusters=np.array([2, 4, 8, 16]),
                 eval_obj=None, n_jobs=None, **kwargs):
        """
        Parameters
        __________

        ensemble_methods: list of clustering methods to use.
                            Should be a list of strings.

        n_clusters: array or int or tuple, dtype int, default [2, 4, 8, 16]
            Array containing the different values of clusters to try,
            or single int specifying the number of clusters,
            or tuple of the form (a, b, c) which specifies a range
            for (x=a; x<b; x+=c)

        eval_obj: Eval or None, default None
            Evaluation object to compare performance of different trials.

        n_jobs: int or None, default None
            Number of jobs to use if multithreading. See
            https://joblib.readthedocs.io/en/latest/generated/
            joblib.Parallel.html.

        **kwargs: dictionary
            Dictionary of parameters that will get passed to obj_def
            when instantiating it.

        """
        methods = _validate_ensemble_methods(ensemble_methods)

        if methods == "default":
            self.methods = ["KMedoids", "KMeans", "Spectral"]
        else:
            self.methods = methods
        self.n_clusters = n_clusters
        self.eval_obj = eval_obj
        self.n_jobs = n_jobs
        self.kwargs = kwargs

    def get(self, x):
        """
        Parameters
        __________
        x: array, shape (n_samples, n_features)
            The data array.

        Returns
        _______
        y: array, shape (n_samples,)
            List of labels that correspond to the best clustering k, as
            evaluated by eval_obj.

        """
        try:
            import Cluster_Ensembles as CE
        except ImportError:
            raise ImportError(
                "Manually install Cluster_Ensembles package if you "
                "wish to run ensemble clustering. For more information ",
                "see here: https://pypi.org/project/Cluster_Ensembles/")

        self.logger.info("Initializing Ensemble Clustering.")
        self.logger.info("Using the following methods:")
        self.logger.info(", ".join(self.methods))

        if len(self.methods) == 1:
            # No need to do ensemble if only one method
            return clu_wrap(self.methods[0])(
                eval_obj=self.eval_obj,
                n_clusters=self.n_clusters,
                n_jobs=self.n_jobs, **self.kwargs).get(x)
        elif len(self.methods) < 1:
            raise ValueError("No methods specified for ensemble clustering.")

        # initialize empty partition matrix
        partitions = np.zeros((len(self.methods), x.shape[0]))
        scores = []

        for i, method in enumerate(self.methods):
            clu_obj = clu_wrap(method)(
                eval_obj=self.eval_obj, n_clusters=self.n_clusters,
                n_jobs=self.n_jobs, **self.kwargs
            )
            partitions[i, :], score = clu_obj.get(x)
            scores.append(np.max(score))

        ensemble_labels = CE.cluster_ensembles(partitions.astype(np.int))

        return ensemble_labels, scores


def cl_ensemble(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        **kwargs):

    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]


    logger.info("Initializing Agglomerative.")

    # get floats
    n_clusters=kwargs.pop('n_clusters')
    if kwargs['distance_threshold']!=None:
        kwargs['distance_threshold']=float(kwargs['distance_threshold'])


    y,score = _get_wrapper(x_to_use, obj_def=Ens_HyperGraph,
                        n_clusters=n_clusters,
                        eval_obj= default_eval_obj, n_jobs=None,
                        attribute_name='n_clusters', **kwargs)

    adata.obs[key] = y
'''
