import igraph as ig
import leidenalg as la
import numpy as np

from ._neighbors import knn_auto, faiss_knn
from controller.cellar.utils.exceptions import InternalError
from app import logger


def cl_ssLeiden(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        partition_type='RBConfigurationVertexPartition', directed=False,
        graph_method='auto', n_neighbors=15, use_weights=False,
        resolution_parameter=1, n_iterations=-1, max_comm_size=0, seed=None,
        fix_membership=True, main_constraints=None, side_constraints=None,
        extras=None):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    constraints = main_constraints if extras['actp'] == 1 else side_constraints
    prefix = 'main' if extras['actp'] == 1 else 'side'

    sources, targets, weights = knn_auto(
        x_to_use, n_neighbors=n_neighbors,
        mode='connectivity', method=graph_method)

    if not use_weights:
        weights = None

    gg = ig.Graph(directed=directed)
    gg.add_vertices(x_to_use.shape[0])
    gg.add_edges(list(zip(list(sources), list(targets))))

    initial_membership = np.arange(gg.vcount())
    is_membership_fixed = np.zeros(gg.vcount())

    labels = adata.obs[key].to_numpy()

    if constraints is not None:
        for constraint in constraints:
            if constraint.startswith(prefix + '-cluster'):
                cluster_id = int(constraint[len(prefix + '-cluster'):])
                cell_indices = np.where(labels == cluster_id)[0]
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

    diff = optimiser.optimise_partition(
        part, n_iterations, is_membership_fixed.tolist())

    part.renumber_communities()

    adata.obs[key] = np.array(part.membership, dtype=int)
