from sklearn.decomposition import PCA, TruncatedSVD, KernelPCA
from sklearn.manifold import MDS
from umap import UMAP


def cl_PCA(adata, key, x_to_use, **kwargs):
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm['x_emb']

    if 'n_components' not in kwargs:
        kwargs['n_components'] = 2

    pca = PCA(**kwargs)

    adata.obsm[key] = pca.fit_transform(x_to_use)


def cl_TruncatedSVD(adata, key, x_to_use, **kwargs):
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm['x_emb']

    if 'n_components' not in kwargs:
        kwargs['n_components'] = 2

    tsvd = TruncatedSVD(**kwargs)

    adata.obsm[key] = tsvd.fit_transform(x_to_use)


def cl_kPCA(adata, key, x_to_use, **kwargs):
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm['x_emb']

    if 'n_components' not in kwargs:
        kwargs['n_components'] = 2

    kpca = KernelPCA(**kwargs)

    adata.obsm[key] = kpca.fit_transform(x_to_use)


def cl_MDS(adata, key, x_to_use, **kwargs):
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm['x_emb']

    if 'n_components' not in kwargs:
        kwargs['n_components'] = 2

    mds = MDS(**kwargs)

    adata.obsm[key] = mds.fit_transform(x_to_use)


def cl_UMAP(adata, key, x_to_use, **kwargs):
    for k in kwargs:
        if kwargs[k] == '':
            kwargs[k] = None

    if x_to_use == 'x':
        x_to_use = adata.X
    else:
        x_to_use = adata.obsm['x_emb']

    if 'n_components' not in kwargs:
        kwargs['n_components'] = 2

    umap = UMAP(**kwargs)

    adata.obsm[key] = umap.fit_transform(x_to_use)


def clear_x_emb_dependends(adata):
    if 'x_emb_2d' in adata.obsm:
        adata.obsm.pop('x_emb_2d')
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')
