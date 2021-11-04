import igraph as ig
import leidenalg as la
import numpy as np
from scipy.sparse import csr_matrix
from app import logger
from controller.cellar.utils.exceptions import InternalError

from ..utils.exceptions import UserError
from ._neighbors import faiss_knn, knn_auto
from ._cluster import _has_neighbors


def cl_ssLeiden(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        partition_type='RBConfigurationVertexPartition', directed=False,
        graph_method='auto', n_neighbors=15, use_cached_neigh=True,
        use_weights=False, neigh_key='neighs', resolution_parameter=1,
        n_iterations=-1, fix_membership=True,
        main_constraints=None, side_constraints=None, extras=None):
    """
    Runs the leiden clustering algorithm with constraints:
    https://www.nature.com/articles/s41598-019-41695-z

    Parameters
    __________

    adata: anndata.AnnData object
    key: str
        Name of the key to add labels under adata.obs[key].
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
    resolution_parameter: float
        Resolution parameter for leidenalg.
    n_iterations: int
        Maximum number of iterations for leidenalg. If set to -1,
        will run until convergence.
    fix_membership: bool
        Set to True if the membership of the constraints should remain
        fixed during the iterations of the leidenalg.
    main_constraints, side_constraints: lists of str
        Lists of constrained clusters. The elements of the lists can have
        two possible formats. One being prefix + '-cluster' + i
        which will constrain the cluster with ID i, and the other format
        is prefix + '-subset-' sbt, which will constrain the subset
        with name sbt under adata.uns['subsets'].
    extras: dict
        Currently only being used to determine if we should use
        main_constraints (if PLOT 1 is active) or side_constraints.
    """
    config = {}
    neigh_config = {}

    annts = False
    if clear_annotations:
        if 'annotations' in adata.obs:
            annts = True
    if annts:
        annos_cp = adata.obs['annotations'].to_numpy().copy()
        adata.obs['annotations'] = np.array(
            [""] * adata.shape[0], dtype="U200")

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

    constraints = main_constraints if extras['actp'] == 1 else side_constraints
    prefix = 'main' if extras['actp'] == 1 else 'side'

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

    weights = None

    gg = ig.Graph(directed=directed)
    config['directed'] = directed
    gg.add_vertices(x_to_use.shape[0])
    gg.add_edges(list(zip(list(sources), list(targets))))

    initial_membership = np.arange(gg.vcount())
    is_membership_fixed = np.zeros(gg.vcount())

    labels = adata.obs[key].to_numpy()

    if annts:
        ant_dict = {}

    if constraints is not None:
        for constraint in constraints:
            if constraint.startswith(prefix + '-cluster'):
                cluster_id = int(constraint[len(prefix + '-cluster'):])
                cell_indices = np.where(labels == cluster_id)[0]
                if annts:
                    ant_dict[annos_cp[cell_indices][0]] = cell_indices
            else:
                subset_name = constraint[len(prefix + '-subset-'):]
                if subset_name not in adata.uns['subsets']:
                    raise InternalError(
                        f"Subset name {subset_name} not found.")

                cell_indices = adata.uns['subsets'][subset_name]

            first_label = initial_membership[cell_indices[0]]
            initial_membership[cell_indices] = first_label
            if fix_membership:
                is_membership_fixed[cell_indices] = 1

    is_membership_fixed = is_membership_fixed.astype(bool)

    optimiser = la.Optimiser()
    if not hasattr(la, partition_type):
        raise UserError(f"No partition type '{partition_type}' found.")
    if partition_type in [
            'ModularityVertexPartition', 'SurpriseVertexPartition']:
        part = getattr(la, partition_type)(
            gg, initial_membership=initial_membership.tolist(),
            weights=weights)
    else:
        part = getattr(la, partition_type)(
            gg, initial_membership=initial_membership.tolist(),
            weights=weights,
            resolution_parameter=resolution_parameter)
    config['partition_type'] = partition_type
    config['resolution_parameter'] = resolution_parameter
    config['method'] = 'Constrained Leiden'

    diff = optimiser.optimise_partition(
        part, n_iterations, is_membership_fixed.tolist())

    # Start from 0
    part.renumber_communities()

    adata.obs[key] = np.array(part.membership, dtype=int)
    adata.uns[key] = config

    if annts:
        labels = adata.obs[key].to_numpy()
        for key in ant_dict:
            ci = ant_dict[key][0]
            ci_label = labels[ci]
            adata.obs['annotations'][labels == ci_label] = key
