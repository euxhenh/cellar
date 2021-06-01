from sklearn.decomposition import PCA, TruncatedSVD


def cl_PCA(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    pca = PCA(**kwargs)

    adata.obsm[key] = pca.fit_transform(adata.X)


def cl_TruncatedSVD(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    tsvd = TruncatedSVD(**kwargs)

    adata.obsm[key] = tsvd.fit_transform(adata.X)


def clear_x_emb_dependends(adata):
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')
