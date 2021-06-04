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

dim_pca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-PCA-n-components",
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
                            id="dim-PCA-whiten",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("SVD Solver"),
                        dbc.Select(
                            id="dim-PCA-solver",
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
                            id="dim-PCA-random-state",
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
    id="dim-PCA-settings",
    target="dim-PCA-btn",
    trigger="click"
)

dim_pca_settings_keys = {
    'dim-PCA-n-components': 'n_components',
    'dim-PCA-whiten': 'whiten',
    'dim-PCA-solver': 'svd_solver',
    'dim-PCA-random-state': 'random_state'
}


dim_tsvd_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Truncated SVD Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-Truncated-SVD-n-components",
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
                            id="dim-Truncated-SVD-solver",
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
                            id="dim-Truncated-SVD-n-iter",
                            type="number",
                            value=5
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="dim-Truncated-SVD-random-state",
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
    id="dim-Truncated-SVD-settings",
    target="dim-Truncated-SVD-btn",
    trigger="click"
)

dim_tsvd_settings_keys = {
    'dim-Truncated-SVD-n-components': 'n_components',
    'dim-Truncated-SVD-solver': 'algorithm',
    'dim-Truncated-SVD-n-iter': 'n_iter',
    'dim-Truncated-SVD-random-state': 'random_state'
}




# Incremental PCA

dim_ipca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Incremental PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-IPCA-n-components",
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
                            id="dim-IPCA-whiten",
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
    id="dim-IPCA-settings",
    target="dim-IPCA-btn",
    trigger="click"
)

dim_ipca_settings_keys = {
    'dim-IPCA-n-components': 'n_components',
    'dim-IPCA-whiten': 'whiten'
}


# Kernel PCA

dim_kpca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Kernel PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-KPCA-n-components",
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
                            id="dim-KPCA-random-state",
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
    id="dim-KPCA-settings",
    target="dim-KPCA-btn",
    trigger="click"
)

dim_kpca_settings_keys = {
    'dim-KPCA-n-components': 'n_components',
    'dim-KPCA-random-state': 'random_state'
}


# TSNE

dim_TSNE_settings = dbc.Popover(
    [
        dbc.PopoverHeader("TSNE Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-TSNE-n-components",
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
                            id="dim-TSNE-learning_rate",
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
                            id="dim-TSNE-random-state",
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
    id="dim-TSNE-settings",
    target="dim-TSNE-btn",
    trigger="click"
)

dim_tsne_settings_keys = {
    'dim-TSNE-n-components': 'n_components',
    'dim-TSNE-n-components': 'learning_rate',
    'dim-TSNE-random-state': 'random_state'
}


# MDS

dim_MDS_settings = dbc.Popover(
    [
        dbc.PopoverHeader("MDS Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-MDS-n-components",
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
                            id="dim-MDS-random-state",
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
    id="dim-MDS-settings",
    target="dim-MDS-btn",
    trigger="click"
)

dim_mds_settings_keys = {
    'dim-MDS-n-components': 'n_components',
    'dim-MDS-random-state': 'random_state'
}


dim_isomap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Isomap Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-Isomap-n-components",
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
                            id="dim-Isomap-eigen-solver",
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
    id="dim-Isomap-settings",
    target="dim-Isomap-btn",
    trigger="click"
)

dim_isomap_settings_keys = {
    'dim-Isomap-n-components': 'n_components',
    'dim-Isomap-eigen-solver': 'eigen_solver'
}


dim_spectral_embedding_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Spectral Embedding Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-Spectral-Embedding-n-components",
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
                            id="dim-Spectral-Embedding-eigen-solver",
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
                            id="dim-Spectral-Embedding-random-state",
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
    id="dim-Spectral-Embedding-settings",
    target="dim-Spectral-Embedding-btn",
    trigger="click"
)

dim_spectral_embedding_settings_keys = {
    'dim-Spectral-Embedding-n-components': 'n_components',
    'dim-Spectral-Embedding-eigen-solver': 'eigen_solver',
    'dim-Spectral-Embedding-random-state': 'random_state'
}



dim_settings = [dim_pca_settings, dim_tsvd_settings,dim_ipca_settings,dim_kpca_settings,dim_TSNE_settings,dim_MDS_settings, 
                dim_isomap_settings, dim_spectral_embedding_settings ]
dim_settings_keys = {
    'dim-PCA': dim_pca_settings_keys,
    'dim-Truncated-SVD': dim_tsvd_settings_keys,
    'dim-IPCA': dim_ipca_settings_keys,
    'dim-KPCA': dim_kpca_settings_keys,
    'dim-TSNE': dim_tsne_settings_keys,
    'dim-MDS': dim_mds_settings_keys,
    'dim-Isomap': dim_isomap_settings_keys,
    'dim-Spectral-Embedding': dim_spectral_embedding_settings_keys
    
}







