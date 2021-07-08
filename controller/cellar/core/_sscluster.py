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
        main_constraints=None, side_constraints=None, extras=None):
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


def cl_uncertainty(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        n_neighbors=15, seed=None, method='centers', extras=None):

    logger.info('uncertainty, uncertainty, uncertainty')

    fm = 1
    fstd = 1
    # an='a1'
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
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
        margin = min_d - min2_d
        #margin /= dist.sum(axis=1)
        #margin = min_d
        m = margin.mean()
        std = margin.std()

        uncertain = (margin > (m*fm+std*fstd)).astype(int)
        uncertain = np.logical_or(uncertain, kml != adata.obs['labels'])
        uncertain = np.array(uncertain)
        logger.info('total uncertain cells:'+str(uncertain.sum()))
        # subset 999 = uncertain
        subset_labels = (1-uncertain)*adata.obs['labels']+uncertain*-1
        adata.obs['labels'] = subset_labels

    return 901
