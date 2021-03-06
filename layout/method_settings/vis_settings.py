import dash_bootstrap_components as dbc
from dash import dcc, html


vis_pca_settings = dbc.Popover(
    [
        dbc.PopoverHeader("PCA Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
    trigger="legacy"
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
                dbc.Form(
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
                dbc.Form(
                    [
                        dbc.Label("Number of Iterations"),
                        dbc.Input(
                            id="vis-Truncated-SVD-n-iter",
                            type="number",
                            value=5
                        )
                    ]
                ),
                dbc.Form(
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
    trigger="legacy"
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
    trigger="legacy"
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
    trigger="legacy"
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
                dbc.Form(
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
    trigger="legacy"
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


vis_tsne_settings = dbc.Popover(
    [
        dbc.PopoverHeader("t-SNE Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("Perplexity"),
                        dcc.Slider(
                            id="vis-TSNE-perplexity",
                            min=5, max=50, step=1,
                            value=30,
                            marks={5: '5', 30: '30', 50: '50'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Early Exaggeration"),
                        dcc.Slider(
                            id="vis-TSNE-early-exaggeration",
                            min=1, max=100, step=1,
                            value=12,
                            marks={1: '1', 12: '12', 100: '100'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Learning Rate"),
                        dcc.Slider(
                            id="vis-TSNE-learning-rate",
                            min=10, max=1000, step=10,
                            value=200,
                            marks={10: '10', 200: '200', 1000: '1000'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. Iterations"),
                        dcc.Slider(
                            id="vis-TSNE-n-iter",
                            min=250, max=5000, step=10,
                            value=1000,
                            marks={250: '250', 1000: '1000', 5000: '5000'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. Iterations without Progress"),
                        dcc.Slider(
                            id="vis-TSNE-n-iter-wo-progress",
                            min=50, max=1000, step=50,
                            value=300,
                            marks={50: '50', 300: '300', 1000: '1000'},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
                            id="vis-TSNE-metric",
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
                dbc.Form(
                    [
                        dbc.Label("Init"),
                        dbc.RadioItems(
                            id="vis-TSNE-init",
                            options=[
                                {"label": "random", "value": "random"},
                                {"label": "pca", "value": "pca"},
                            ],
                            value="random",
                            inline=True
                        )
                    ]
                ),
                dbc.Form(
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
    trigger="legacy"
)

vis_tsne_settings_keys = {
    'vis-TSNE-perplexity': 'perplexity',
    'vis-TSNE-early-exaggeration': 'early_exaggeration',
    'vis-TSNE-learning-rate': 'learning_rate',
    'vis-TSNE-n-iter': 'n_iter',
    'vis-TSNE-n-iter-wo-progress': 'n_iter_without_progress',
    'vis-TSNE-metric': 'metric',
    'vis-TSNE-init': 'init',
    'vis-TSNE-random-state': 'random_state'
}


vis_diffmap_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Diffusion Map Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("alpha"),
                        dbc.Input(
                            id="vis-diffmap-alpha",
                            type="number",
                            value=0.5
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("No. Nearest Neighbors", html_for="slider"),
                        dcc.Slider(
                            id="vis-diffmap-n-neigh",
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
                            id="vis-diffmap-epsilon",
                            type="text",
                            value='bgh'
                        )
                    ]
                ),
                dbc.Form(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
                            id="vis-diffmap-metric",
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
    id="vis-Diffmap-settings",
    target="vis-Diffmap-btn",
    trigger="legacy"
)

vis_diffmap_settings_keys = {
    'vis-diffmap-alpha': 'alpha',
    'vis-diffmap-n-neigh': 'k',
    'vis-diffmap-epsilon': 'epsilon',
    'vis-diffmap-metric': 'metric'
}


vis_settings = [
    vis_umap_settings,
    vis_tsne_settings,
    vis_pca_settings,
    vis_tsvd_settings,
    vis_kpca_settings,
    vis_mds_settings,
    vis_diffmap_settings
]

vis_settings_keys = {
    'vis-UMAP': vis_umap_settings_keys,
    'vis-TSNE': vis_tsne_settings_keys,
    'vis-PCA': vis_pca_settings_keys,
    'vis-Truncated-SVD': vis_tsvd_settings_keys,
    'vis-Kernel-PCA': vis_kpca_settings_keys,
    'vis-MDS': vis_mds_settings_keys,
    'vis-Diffmap': vis_diffmap_settings_keys,
}
