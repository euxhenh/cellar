import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

'''Catalogue'''
# PCA
# Truncated SVD
# Incremental PCA
# Kernel PCA
# TSNE
# MDS
# Isomap
# Spectral Embedding

vis_pca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-PCA-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Whiten"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=False,
                            id="vis-PCA-whiten",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("SVD Solver"),
                        dbc.Select(
                            id="vis-PCA-solver",
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "full", "value": "full"},
                                {"label": "arpack", "value": "arpack"},
                                {"label": "randomized", "value": "randomized"}
                            ],
                            value="auto"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-PCA-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-PCA-settings",
    target="vis-PCA-btn",
    trigger="click"
)

vis_pca_settings_keys = {
    'vis-PCA-n-components': 'n_components',
    'vis-PCA-whiten': 'whiten',
    'vis-PCA-solver': 'svd_solver',
    'vis-PCA-random-state': 'random_state'
}


vis_tsvd_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Truncated SVD Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-Truncated-SVD-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Algorithm"),
                        dbc.Select(
                            id="vis-Truncated-SVD-solver",
                            options=[
                                {"label": "arpack", "value": "arpack"},
                                {"label": "randomized", "value": "randomized"}
                            ],
                            value="randomized"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Iterations"),
                        dbc.Input(
                            id="vis-Truncated-SVD-n-iter",
                            type="number",
                            value=5
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-Truncated-SVD-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.TruncatedSVD.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Truncated-SVD-settings",
    target="vis-Truncated-SVD-btn",
    trigger="click"
)

vis_tsvd_settings_keys = {
    'vis-Truncated-SVD-n-components': 'n_components',
    'vis-Truncated-SVD-solver': 'algorithm',
    'vis-Truncated-SVD-n-iter': 'n_iter',
    'vis-Truncated-SVD-random-state': 'random_state'
}




# Incremental PCA

vis_ipca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Incremental PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-IPCA-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Whiten"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=False,
                            id="vis-IPCA-whiten",
                            inline=True
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.decomposition.IncrementalPCA.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.IncrementalPCA.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-IPCA-settings",
    target="vis-IPCA-btn",
    trigger="click"
)

vis_ipca_settings_keys = {
    'vis-IPCA-n-components': 'n_components',
    'vis-IPCA-whiten': 'whiten'
}


# Kernel PCA

vis_kpca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Kernel PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-KPCA-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-KPCA-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html#sklearn.decomposition.KernelPCA",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html#sklearn.decomposition.KernelPCA",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-KPCA-settings",
    target="vis-KPCA-btn",
    trigger="click"
)

vis_kpca_settings_keys = {
    'vis-KPCA-n-components': 'n_components',
    'vis-KPCA-random-state': 'random_state'
}


# TSNE

vis_TSNE_settings = dbc.Popover(
    [
        dbc.PopoverHeader("TSNE Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-TSNE-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Learning Rate", html_for="slider"),
                        dcc.Slider(
                            id="vis-TSNE-learning_rate",
                            min=10, max=1000, step=1,
                            value=200,
                            marks={i: str(i) for i in range(10, 1001, 100)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-TSNE-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.TSNE.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-TSNE-settings",
    target="vis-TSNE-btn",
    trigger="click"
)

vis_tsne_settings_keys = {
    'vis-TSNE-n-components': 'n_components',
    'vis-TSNE-n-components': 'learning_rate',
    'vis-TSNE-random-state': 'random_state'
}


# MDS

vis_MDS_settings = dbc.Popover(
    [
        dbc.PopoverHeader("MDS Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-MDS-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),

                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-MDS-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.MDS.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-MDS-settings",
    target="vis-MDS-btn",
    trigger="click"
)

vis_mds_settings_keys = {
    'vis-MDS-n-components': 'n_components',
    'vis-MDS-random-state': 'random_state'
}


vis_isomap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Isomap Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-Isomap-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Eigen Solver"),
                        dbc.Select(
                            id="vis-Isomap-eigen-solver",
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "arpack", "value": "arpack"},
                                {"label": "dense", "value": "dense"}
                            ],
                            value="auto"
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.Isomap.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Isomap-settings",
    target="vis-Isomap-btn",
    trigger="click"
)

vis_isomap_settings_keys = {
    'vis-Isomap-n-components': 'n_components',
    'vis-Isomap-eigen-solver': 'eigen_solver'
}


vis_spectral_embedding_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Spectral Embedding Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-Spectral-Embedding-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Eigen Solver"),
                        dbc.Select(
                            id="vis-Spectral-Embedding-eigen-solver",
                            options=[
                                {"label": "arpack", "value": "arpack"},
                                {"label": "lobpcg", "value": "lobpcg"},
                                {"label": "amg", "value": "amg"},
                                
                            ],
                            value="arpack"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-Spectral-Embedding-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.manifold.SpectralEmbedding.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.manifold.SpectralEmbedding.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Spectral-Embedding-settings",
    target="vis-Spectral-Embedding-btn",
    trigger="click"
)

vis_spectral_embedding_settings_keys = {
    'vis-Spectral-Embedding-n-components': 'n_components',
    'vis-Spectral-Embedding-eigen-solver': 'eigen_solver',
    'vis-Spectral-Embedding-random-state': 'random_state'
}



vis_feature_agglomeration_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Feature Agglomeration Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-Feature-Agglomeration-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.cluster.FeatureAgglomeration.html#sklearn.cluster.FeatureAgglomeration.fit_transform",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.cluster.FeatureAgglomeration.html#sklearn.cluster.FeatureAgglomeration.fit_transform",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Feature-Agglomeration-settings",
    target="vis-Feature-Agglomeration-btn",
    trigger="click"
)

vis_feature_agglomeration_settings_keys = {
    'vis-Feature-Agglomeration-n-components': 'n_clusters'
}


vis_diffusion_map_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Diffusion Map Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-Diffusion-Map-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("pydiffmap.readthedocs.io/en/master/reference/diffusion_map.html",
                               href="https://pydiffmap.readthedocs.io/en/master/reference/diffusion_map.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Diffusion-Map-settings",
    target="vis-Diffusion-Map-btn",
    trigger="click"
)

vis_diffusion_map_settings_keys = {
    'vis-Diffusion-Map-n-components': 'n_evecs'
}


vis_umap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("UMAP Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="vis-UMAP-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.UMAP.fit_transform",
                               href="https://umap-learn.readthedocs.io/en/latest/api.html#umap.umap_.UMAP.fit_transform",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-UMAP-settings",
    target="vis-UMAP-btn",
    trigger="click"
)

vis_umap_settings_keys = {
    'vis-UMAP-n-components': 'n_components'
}


vis_settings = [vis_pca_settings, vis_tsvd_settings,vis_ipca_settings,vis_kpca_settings,vis_TSNE_settings,vis_MDS_settings, 
                vis_isomap_settings, vis_spectral_embedding_settings, vis_feature_agglomeration_settings, vis_diffusion_map_settings,
                vis_umap_settings]
vis_settings_keys = {
    'vis-PCA': vis_pca_settings_keys,
    'vis-Truncated-SVD': vis_tsvd_settings_keys,
    'vis-IPCA': vis_ipca_settings_keys,
    'vis-KPCA': vis_kpca_settings_keys,
    'vis-TSNE': vis_tsne_settings_keys,
    'vis-MDS': vis_mds_settings_keys,
    'vis-Isomap': vis_isomap_settings_keys,
    'vis-Spectral-Embedding': vis_spectral_embedding_settings_keys,
    'vis-Feature-Agglomeration': vis_feature_agglomeration_settings_keys,
    'vis-Diffusion-Map': vis_diffusion_map_settings_keys,
    'vis-UMAP': vis_umap_settings_keys
}







