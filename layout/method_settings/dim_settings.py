import dash_bootstrap_components as dbc
from dash import dcc, html


dim_pca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-PCA-n-components",
                            min=3, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                        html.A("scikit-learn.org/stable/modules/generated/"
                               "sklearn.decomposition.PCA.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.decomposition.PCA.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-PCA-settings",
    target="dim-PCA-btn",
    trigger="legacy"
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
                    [
                        dbc.Label("Number of Iterations"),
                        dbc.Input(
                            id="dim-Truncated-SVD-n-iter",
                            type="number",
                            value=5
                        )
                    ]
                ),
                dbc.Form(
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
                        html.A("scikit-learn.org/stable/modules/generated/"
                               "sklearn.decomposition.TruncatedSVD.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.decomposition."
                               "TruncatedSVD.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-Truncated-SVD-settings",
    target="dim-Truncated-SVD-btn",
    trigger="legacy"
)

dim_tsvd_settings_keys = {
    'dim-Truncated-SVD-n-components': 'n_components',
    'dim-Truncated-SVD-solver': 'algorithm',
    'dim-Truncated-SVD-n-iter': 'n_iter',
    'dim-Truncated-SVD-random-state': 'random_state'
}


dim_kpca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Kernel PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-kPCA-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Kernel"),
                        dbc.Select(
                            id="dim-kPCA-kernel",
                            options=[
                                {"label": "linear", "value": "linear"},
                                {"label": "poly", "value": "poly"},
                                {"label": "rbf", "value": "rbf"},
                                {"label": "sigmoid", "value": "sigmoid"},
                                {"label": "cosine", "value": "cosine"}
                            ],
                            value="linear"
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Gamma (only for rbf, poly, sigmoid)"),
                        dbc.Input(
                            id="dim-kPCA-gamma",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Degree (only for poly)",
                                  html_for="slider"),
                        dcc.Slider(
                            id="dim-kPCA-degree",
                            min=1, max=10, step=1,
                            value=3,
                            marks={i: str(i) for i in range(1, 11)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Coef0 (only for poly, sigmoid)"),
                        dbc.Input(
                            id="dim-kPCA-coef0",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Eigen Solver"),
                        dbc.Select(
                            id="dim-kPCA-solver",
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "dense", "value": "dense"},
                                {"label": "arpack", "value": "arpack"}
                            ],
                            value="auto"
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="dim-kPCA-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/"
                               "sklearn.decomposition.KernelPCA.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.decomposition."
                               "KernelPCA.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-Kernel-PCA-settings",
    target="dim-Kernel-PCA-btn",
    trigger="legacy"
)

dim_kpca_settings_keys = {
    'dim-kPCA-n-components': 'n_components',
    'dim-kPCA-kernel': 'kernel',
    'dim-kPCA-gamma': 'gamma',
    'dim-kPCA-degree': 'degree',
    'dim-kPCA-coef0': 'coef0',
    'dim-kPCA-solver': 'eigen_solver',
    'dim-kPCA-random-state': 'random_state'
}


dim_mds_settings = dbc.Popover(
    [
        dbc.PopoverHeader("MDS Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
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
                dbc.Form(
                    [
                        dbc.Label("Metric"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="dim-MDS-metric",
                            inline=True
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("N Init", html_for="slider"),
                        dcc.Slider(
                            id="dim-MDS-ninit",
                            min=1, max=30, step=1,
                            value=4,
                            marks={i: str(i) for i in range(1, 31, 3)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Max Iter", html_for="slider"),
                        dcc.Slider(
                            id="dim-MDS-max-iter",
                            min=10, max=1000, step=10,
                            value=300,
                            marks={i: str(i) for i in
                                   [10] + list(range(100, 1001, 100))},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
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
                        html.A("scikit-learn.org/stable/modules/generated/"
                               "sklearn.manifold.MDS.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.manifold.MDS.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-MDS-settings",
    target="dim-MDS-btn",
    trigger="legacy"
)

dim_mds_settings_keys = {
    'dim-MDS-n-components': 'n_components',
    'dim-MDS-metric': 'metric',
    'dim-MDS-ninit': 'n_init',
    'dim-MDS-max-iter': 'max_iter',
    'dim-MDS-random-state': 'random_state'
}


dim_diffmap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Diffusion Map Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-diffmap-n-evecs",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("alpha"),
                        dbc.Input(
                            id="dim-diffmap-alpha",
                            type="number",
                            value=0.5
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. Nearest Neighbors", html_for="slider"),
                        dcc.Slider(
                            id="dim-diffmap-n-neigh",
                            min=2, max=100, step=1,
                            value=64,
                            marks={2: '2', 64: '64', 100: '100'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("epsilon"),
                        dbc.Input(
                            id="dim-diffmap-epsilon",
                            type="text",
                            value='bgh'
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
                            id="dim-diffmap-metric",
                            options=[
                                {"label": "euclidean", "value": "euclidean"},
                                {"label": "manhattan", "value": "manhattan"},
                                {"label": "chebyshev", "value": "chebyshev"},
                                {"label": "minkowski", "value": "minkowski"},
                                {"label": "canberra", "value": "canberra"},
                                {"label": "braycurtis", "value": "braycurtis"},
                                {"label": "mahalanobis",
                                 "value": "mahalanobis"},
                                {"label": "wminkowski", "value": "wminkowski"},
                                {"label": "seuclidean", "value": "seuclidean"},
                                {"label": "cosine", "value": "cosine"},
                                {"label": "correlation",
                                 "value": "correlation"},
                                {"label": "haversine", "value": "haversine"},
                                {"label": "hamming", "value": "hamming"},
                                {"label": "jaccard", "value": "jaccard"},
                                {"label": "dice", "value": "dice"},
                                {"label": "russelrao", "value": "russelrao"},
                                {"label": "kulsinski", "value": "kulsinski"},
                                {"label": "ll_dirichlet", "value":
                                 "ll_dirichlet"},
                                {"label": "hellinger", "value": "hellinger"},
                                {"label": "rogerstanimoto",
                                    "value": "rogerstanimoto"},
                                {"label": "sokalmichener",
                                    "value": "sokalmichener"},
                                {"label": "sokalsneath",
                                 "value": "sokalsneath"},
                                {"label": "yule", "value": "yule"}
                            ],
                            value="euclidean"
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("pydiffmap.readthedocs.io/en/master"
                               "/reference/diffusion_map.html",
                               href="https://pydiffmap.readthedocs.io/en/"
                               "master/reference/diffusion_map.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-Diffmap-settings",
    target="dim-Diffmap-btn",
    trigger="legacy"
)

dim_diffmap_settings_keys = {
    'dim-diffmap-n-evecs': 'n_evecs',
    'dim-diffmap-alpha': 'alpha',
    'dim-diffmap-n-neigh': 'k',
    'dim-diffmap-epsilon': 'epsilon',
    'dim-diffmap-metric': 'metric'
}


dim_umap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("UMAP Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-UMAP-n-components",
                            min=2, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. of neighbors", html_for="slider"),
                        dcc.Slider(
                            id="dim-UMAP-n-neighbors",
                            min=2, max=100, step=1,
                            value=15,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
                            id="dim-UMAP-metric",
                            options=[
                                {"label": "euclidean", "value": "euclidean"},
                                {"label": "manhattan", "value": "manhattan"},
                                {"label": "chebyshev", "value": "chebyshev"},
                                {"label": "minkowski", "value": "minkowski"},
                                {"label": "canberra", "value": "canberra"},
                                {"label": "braycurtis", "value": "braycurtis"},
                                {"label": "mahalanobis",
                                 "value": "mahalanobis"},
                                {"label": "wminkowski", "value": "wminkowski"},
                                {"label": "seuclidean", "value": "seuclidean"},
                                {"label": "cosine", "value": "cosine"},
                                {"label": "correlation",
                                 "value": "correlation"},
                                {"label": "haversine", "value": "haversine"},
                                {"label": "hamming", "value": "hamming"},
                                {"label": "jaccard", "value": "jaccard"},
                                {"label": "dice", "value": "dice"},
                                {"label": "russelrao", "value": "russelrao"},
                                {"label": "kulsinski", "value": "kulsinski"},
                                {"label": "ll_dirichlet", "value":
                                 "ll_dirichlet"},
                                {"label": "hellinger", "value": "hellinger"},
                                {"label": "rogerstanimoto",
                                    "value": "rogerstanimoto"},
                                {"label": "sokalmichener",
                                    "value": "sokalmichener"},
                                {"label": "sokalsneath",
                                 "value": "sokalsneath"},
                                {"label": "yule", "value": "yule"}
                            ],
                            value="euclidean"
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. of epochs"),
                        dbc.Input(
                            id="dim-UMAP-n-epochs",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Learning Rate"),
                        dbc.Input(
                            id="dim-UMAP-lr",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Initialization"),
                        dbc.Select(
                            id="dim-UMAP-init",
                            options=[
                                {"label": "spectral", "value": "spectral"},
                                {"label": "random", "value": "random"}
                            ],
                            value="spectral"
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Minimum Distance"),
                        dbc.Input(
                            id="dim-UMAP-min-dist",
                            placeholder="Leave empty if None",
                            type="number",
                            value=0.1
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Spread"),
                        dbc.Input(
                            id="dim-UMAP-spread",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="dim-UMAP-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("umap-learn.readthedocs.io/en/latest/api.html",
                               href="https://umap-learn.readthedocs.io/en/"
                               "latest/api.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-UMAP-settings",
    target="dim-UMAP-btn",
    trigger="legacy"
)

dim_umap_settings_keys = {
    'dim-UMAP-n-neighbors': 'n_neighbors',
    'dim-UMAP-n-components': 'n_components',
    'dim-UMAP-metric': 'metric',
    'dim-UMAP-n-epochs': 'n_epochs',
    'dim-UMAP-lr': 'learning_rate',
    'dim-UMAP-init': 'init',
    'dim-UMAP-min-dist': 'min_dist',
    'dim-UMAP-spread': 'spread',
    'dim-UMAP-random-state': 'random_state'
}

dim_cistopic_settings = dbc.Popover(
    [
        dbc.PopoverHeader("cisTopic Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="dim-cisTopic-n-components",
                            min=10, max=100, step=1,
                            value=40,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Iterations", html_for="slider"),
                        dcc.Slider(
                            id="dim-cisTopic-iterations",
                            min=50, max=500, step=10,
                            value=150,
                            marks={i: str(i) for i in range(50, 501, 50)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("htmlpreview.github.io/?https://github.com"
                               "/aertslab/cisTopic/blob/master/"
                               "vignettes/WarpLDA_10X_workflow.html",
                               href="http://htmlpreview.github.io/?https"
                               "://github.com/aertslab/cisTopic"
                               "/blob/master/vignettes/"
                               "WarpLDA_10X_workflow.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="dim-cisTopic-settings",
    target="dim-cisTopic-btn",
    trigger="legacy"
)

dim_cistopic_settings_keys = {
    'dim-cisTopic-n-components': 'topics',
    'dim-cisTopic-iterations': 'iterations'
}


dim_settings = [
    dim_pca_settings,
    dim_tsvd_settings,
    dim_kpca_settings,
    dim_mds_settings,
    dim_umap_settings,
    dim_diffmap_settings,
    dim_cistopic_settings
]

dim_settings_keys = {
    'dim-PCA': dim_pca_settings_keys,
    'dim-Truncated-SVD': dim_tsvd_settings_keys,
    'dim-Kernel-PCA': dim_kpca_settings_keys,
    'dim-MDS': dim_mds_settings_keys,
    'dim-UMAP': dim_umap_settings_keys,
    'dim-Diffmap': dim_diffmap_settings_keys,
    'dim-cisTopic': dim_cistopic_settings_keys
}
