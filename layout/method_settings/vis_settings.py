import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


vis_pca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("PCA Settings"),
        dbc.PopoverBody(
            [
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
    'vis-Truncated-SVD-solver': 'algorithm',
    'vis-Truncated-SVD-n-iter': 'n_iter',
    'vis-Truncated-SVD-random-state': 'random_state'
}


vis_kpca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Kernel PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Kernel"),
                        dbc.Select(
                            id="vis-kPCA-kernel",
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
                dbc.FormGroup(
                    [
                        dbc.Label("Gamma (only for rbf, poly, sigmoid)"),
                        dbc.Input(
                            id="vis-kPCA-gamma",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Degree (only for poly)",
                                  html_for="slider"),
                        dcc.Slider(
                            id="vis-kPCA-degree",
                            min=1, max=10, step=1,
                            value=3,
                            marks={i: str(i) for i in range(1, 11)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Coef0 (only for poly, sigmoid)"),
                        dbc.Input(
                            id="vis-kPCA-coef0",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Eigen Solver"),
                        dbc.Select(
                            id="vis-kPCA-solver",
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "dense", "value": "dense"},
                                {"label": "arpack", "value": "arpack"}
                            ],
                            value="auto"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-kPCA-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.KernelPCA.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="vis-Kernel-PCA-settings",
    target="vis-Kernel-PCA-btn",
    trigger="click"
)

vis_kpca_settings_keys = {
    'vis-kPCA-kernel': 'kernel',
    'vis-kPCA-gamma': 'gamma',
    'vis-kPCA-degree': 'degree',
    'vis-kPCA-coef0': 'coef0',
    'vis-kPCA-solver': 'eigen_solver',
    'vis-kPCA-random-state': 'random_state'
}


vis_mds_settings = dbc.Popover(
    [
        dbc.PopoverHeader("MDS Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Metric"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="vis-MDS-metric",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("N Init", html_for="slider"),
                        dcc.Slider(
                            id="vis-MDS-ninit",
                            min=1, max=30, step=1,
                            value=4,
                            marks={i: str(i) for i in range(1, 31, 3)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Max Iter", html_for="slider"),
                        dcc.Slider(
                            id="vis-MDS-max-iter",
                            min=10, max=1000, step=10,
                            value=300,
                            marks={i: str(i)
                                   for i in [10] + list(range(100, 1001, 100))},
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
    'vis-MDS-metric': 'metric',
    'vis-MDS-ninit': 'n_init',
    'vis-MDS-max-iter': 'max_iter',
    'vis-MDS-random-state': 'random_state'
}


vis_umap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("UMAP Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("No. of neighbors", html_for="slider"),
                        dcc.Slider(
                            id="vis-UMAP-n-neighbors",
                            min=2, max=100, step=1,
                            value=15,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
                            id="vis-UMAP-metric",
                            options=[
                                {"label": "euclidean", "value": "euclidean"},
                                {"label": "manhattan", "value": "manhattan"},
                                {"label": "chebyshev", "value": "chebyshev"},
                                {"label": "minkowski", "value": "minkowski"},
                                {"label": "canberra", "value": "canberra"},
                                {"label": "braycurtis", "value": "braycurtis"},
                                {"label": "mahalanobis", "value": "mahalanobis"},
                                {"label": "wminkowski", "value": "wminkowski"},
                                {"label": "seuclidean", "value": "seuclidean"},
                                {"label": "cosine", "value": "cosine"},
                                {"label": "correlation", "value": "correlation"},
                                {"label": "haversine", "value": "haversine"},
                                {"label": "hamming", "value": "hamming"},
                                {"label": "jaccard", "value": "jaccard"},
                                {"label": "dice", "value": "dice"},
                                {"label": "russelrao", "value": "russelrao"},
                                {"label": "kulsinski", "value": "kulsinski"},
                                {"label": "ll_dirichlet", "value": "ll_dirichlet"},
                                {"label": "hellinger", "value": "hellinger"},
                                {"label": "rogerstanimoto",
                                    "value": "rogerstanimoto"},
                                {"label": "sokalmichener",
                                    "value": "sokalmichener"},
                                {"label": "sokalsneath", "value": "sokalsneath"},
                                {"label": "yule", "value": "yule"}
                            ],
                            value="euclidean"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("No. of epochs"),
                        dbc.Input(
                            id="vis-UMAP-n-epochs",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Learning Rate"),
                        dbc.Input(
                            id="vis-UMAP-lr",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Initialization"),
                        dbc.Select(
                            id="vis-UMAP-init",
                            options=[
                                {"label": "spectral", "value": "spectral"},
                                {"label": "random", "value": "random"}
                            ],
                            value="spectral"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Minimum Distance"),
                        dbc.Input(
                            id="vis-UMAP-min-dist",
                            placeholder="Leave empty if None",
                            type="number",
                            value=0.1
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Spread"),
                        dbc.Input(
                            id="vis-UMAP-spread",
                            placeholder="Leave empty if None",
                            type="number",
                            value=1.0
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="vis-UMAP-random-state",
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
                               href="https://umap-learn.readthedocs.io/en/latest/api.html",
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
    'vis-UMAP-n-neighbors': 'n_neighbors',
    'vis-UMAP-metric': 'metric',
    'vis-UMAP-n-epochs': 'n_epochs',
    'vis-UMAP-lr': 'learning_rate',
    'vis-UMAP-init': 'init',
    'vis-UMAP-min-dist': 'min_dist',
    'vis-UMAP-spread': 'spread',
    'vis-UMAP-random-state': 'random_state'
}


vis_settings = [
    vis_umap_settings,
    vis_pca_settings,
    # vis_tsvd_settings,
    vis_kpca_settings,
    vis_mds_settings
]

vis_settings_keys = {
    'vis-UMAP': vis_umap_settings_keys,
    'vis-PCA': vis_pca_settings_keys,
    # 'vis-Truncated-SVD': vis_tsvd_settings_keys,
    'vis-Kernel-PCA': vis_kpca_settings_keys,
    'vis-MDS': vis_mds_settings_keys
}
