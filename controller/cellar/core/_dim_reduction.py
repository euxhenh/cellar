from sklearn.decomposition import PCA, TruncatedSVD, IncrementalPCA, KernelPCA
from sklearn.manifold import TSNE
from sklearn.manifold import MDS
from sklearn.manifold import Isomap
from sklearn.manifold import SpectralEmbedding
from sklearn.cluster import FeatureAgglomeration
from pydiffmap.diffusion_map import DiffusionMap as dm
from umap import UMAP

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



def cl_IncrementalPCA(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    ipca = IncrementalPCA(**kwargs)

    adata.obsm[key] = ipca.fit_transform(adata.X)


def cl_KernelPCA(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    kpca = KernelPCA(**kwargs)

    adata.obsm[key] = kpca.fit_transform(adata.X)


def cl_TSNE(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    tsne = TSNE(**kwargs)

    adata.obsm[key] = tsne.fit_transform(adata.X)



def cl_MDS(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    mds = MDS(**kwargs)

    adata.obsm[key] = mds.fit_transform(adata.X)



def cl_Isomap(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    isomap = Isomap(**kwargs)

    adata.obsm[key] = isomap.fit_transform(adata.X)



def cl_SpectralEmbedding(adata, key='x_emb', **kwargs):
    if 'random_state' in kwargs:
        if kwargs['random_state'] == '':
            kwargs['random_state'] = None

    s_e_ = SpectralEmbedding(**kwargs)

    adata.obsm[key] = s_e_.fit_transform(adata.X)





def cl_FeatureAgglomeration(adata, key='x_emb', **kwargs):

    f_a_ = FeatureAgglomeration(**kwargs)
    adata.obsm[key] = f_a_.fit_transform(adata.X)


def cl_dm(adata, key='x_emb', **kwargs):

    adata.obsm[key] = dm.from_sklearn(**kwargs).fit_transform(adata.X)
    


def cl_umap(adata, key='x_emb', **kwargs):

    u = UMAP(**kwargs)
    adata.obsm[key] = u.fit_transform(adata.X)



def clear_x_emb_dependends(adata):
    if 'labels' in adata.obs:
        adata.obs.pop('labels')
    if 'annotations' in adata.obs:
        adata.obs.pop('annotations')



