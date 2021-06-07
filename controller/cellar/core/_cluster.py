import igraph as ig
import leidenalg as la
import numpy as np

from ._neighbors import knn_auto


def cl_Leiden(
        adata, key='labels', x_to_use='x_emb', clear_annotations=True,
        **kwargs):
    if clear_annotations:
        if 'annotations' in adata.obs:
            adata.obs.pop('annotations')
    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm[x_to_use]

    sources, targets, weights = knn_auto(
        x_to_use, n_neighbors=15, mode='distance')

    gg = ig.Graph(directed=True)
    gg.add_vertices(x_to_use.shape[0])
    gg.add_edges(list(zip(list(sources), list(targets))))

    part = la.find_partition(
        gg, la.RBConfigurationVertexPartition,
        weights=weights,
        n_iterations=-1,
        resolution_parameter=1)

    adata.obs[key] = np.array(part.membership, dtype=int)
