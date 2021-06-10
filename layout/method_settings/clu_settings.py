import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


clu_leiden_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Leiden Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Partition Type"),
                        dbc.Select(
                            id="clu-Leiden-partition-type",
                            options=[
                                {"label": "RBConfigurationVertexPartition",
                                    "value": "RBConfigurationVertexPartition"},
                                {"label": "ModularityVertexPartition",
                                    "value": "ModularityVertexPartition"},
                                {"label": "RBERVertexPartition",
                                    "value": "RBERVertexPartition"},
                                {"label": "CPMVertexPartition",
                                    "value": "CPMVertexPartition"},
                                {"label": "SurpriseVertexPartition",
                                    "value": "SurpriseVertexPartition"}
                            ],
                            value="RBConfigurationVertexPartition"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Iterations"),
                        dbc.Input(
                            id="clu-Leiden-n-iter",
                            placeholder="Set to -1 to run until convergence",
                            type="number",
                            value=-1
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label(
                            "Resolution (only for RBConfiguration, RBER, CPM)",
                            html_for="slider"),
                        dcc.Slider(
                            id="clu-Leiden-resolution",
                            min=0.01, max=10, step=0.01,
                            value=1,
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Max Community Size"),
                        dbc.Input(
                            id="clu-Leiden-max-comm",
                            placeholder="Set to 0 to allow any size",
                            type="number",
                            value=0
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Graph Construction Method"),
                        dbc.Select(
                            id="clu-Leiden-graph-method",
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "approximate neighbors",
                                 "value": "approximate"},
                                {"label": "full (slow for large datasets)",
                                 "value": "full"}
                            ],
                            value="auto"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label(
                            "Number of Neighbors for the Graph",
                            html_for="slider"),
                        dcc.Slider(
                            id="clu-Leiden-nneigh",
                            min=1, max=100, step=1,
                            value=15,
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Construct Weighted Graph"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="clu-Leiden-weights",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Construct Directed Graph"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="clu-Leiden-directed",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="clu-Leiden-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=""
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("leidenalg.readthedocs.io/en/stable/reference"
                               ".html#reference",
                               href="https://leidenalg.readthedocs.io/en/"
                               "stable/reference.html#reference",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="clu-Leiden-settings",
    target="clu-Leiden-btn",
    trigger="click"
)

clu_leiden_settings_keys = {
    'clu-Leiden-partition-type': 'partition_type',
    'clu-Leiden-n-iter': 'n_iterations',
    'clu-Leiden-resolution': 'resolution_parameter',
    'clu-Leiden-max-comm': 'max_comm_size',
    'clu-Leiden-graph-method': 'graph_method',
    'clu-Leiden-nneigh': 'n_neighbors',
    'clu-Leiden-weights': 'use_weights',
    'clu-Leiden-directed': 'directed',
    'clu-Leiden-random-state': 'seed'
}


clu_kmeans_settings = dbc.Popover(
    [
        dbc.PopoverHeader("KMeans Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Clusters"),
                        dbc.Input(
                            id="clu-KMeans-n-clusters",
                            type="text",
                            value='[4, 6, 8, 10, 12, 14, 16]'
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Initialization"),
                        dbc.Select(
                            options=[
                                {"label": "k-means++", "value": "k-means++"},
                                {"label": "random", "value": "random"}
                            ],
                            id="clu-KMeans-init",
                            value='k-means++'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("No. of Initializations", html_for="slider"),
                        dcc.Slider(
                            id="clu-KMeans-n-init",
                            min=1, max=100, step=1,
                            value=10,
                            marks={i: str(i)
                                   for i in list(range(10, 101, 10))},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Max Iterations", html_for="slider"),
                        dcc.Slider(
                            id="clu-KMeans-max-iter",
                            min=10, max=1000, step=10,
                            value=300,
                            marks={i: str(i) for i in
                                   [10] + list(range(100, 1001, 100))},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Algorithm"),
                        dbc.Select(
                            options=[
                                {"label": "auto", "value": "auto"},
                                {"label": "full", "value": "full"},
                                {"label": "elkan", "value": "elkan"},
                            ],
                            id="clu-KMeans-algorithm",
                            value='auto'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="clu-KMeans-random-state",
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
                               "sklearn.cluster.KMeans.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.cluster.KMeans.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="clu-KMeans-settings",
    target="clu-KMeans-btn",
    trigger="click"
)


clu_kmeans_settings_keys = {
    'clu-KMeans-n-clusters': 'n_clusters',
    'clu-KMeans-init': 'init',
    'clu-KMeans-n-init': 'n_init',
    'clu-KMeans-max-iter': 'max_iter',
    'clu-KMeans-algorithm': 'algorithm',
    'clu-KMeans-random-state': 'random_state'
}


clu_kmedoids_settings = dbc.Popover(
    [
        dbc.PopoverHeader("KMedoids Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Clusters"),
                        dbc.Input(
                            id="clu-KMedoids-n-clusters",
                            type="text",
                            value='[4, 6, 8, 10, 12, 14, 16]'
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Metric"),
                        dbc.Select(
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
                                {"label": "sqeuclidean",
                                 "value": "sqeuclidean"},
                                {"label": "nan_euclidean",
                                    "value": "nan_euclidean"},
                                {"label": "cosine", "value": "cosine"},
                                {"label": "correlation",
                                 "value": "correlation"},
                                {"label": "hamming", "value": "hamming"},
                                {"label": "haversine", "value": "haversine"},
                                {"label": "jaccard", "value": "jaccard"},
                                {"label": "dice", "value": "dice"},
                                {"label": "russellrao", "value": "russelrao"},
                                {"label": "kulsinski", "value": "kulsinski"},
                                {"label": "rogerstanimoto",
                                 "value": "rogerstanimoto"},
                                {"label": "sokalmichener",
                                 "value": "sokalmichener"},
                                {"label": "sokalsneath",
                                 "value": "sokalsneath"},
                                {"label": "yule", "value": "yule"}
                            ],
                            id="clu-KMedoids-metric",
                            value='euclidean'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Method"),
                        dbc.Select(
                            options=[
                                {"label": "alternate", "value": "alternate"},
                                {"label": "pam", "value": "pam"},
                            ],
                            id="clu-KMedoids-method",
                            value='alternate'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Initialization"),
                        dbc.Select(
                            options=[
                                {"label": "build", "value": "build"},
                                {"label": "random", "value": "random"},
                                {"label": "heuristic", "value": "heuristic"},
                                {"label": "k-medoids++",
                                 "value": "k-medoids++"},
                            ],
                            id="clu-KMedoids-init",
                            value='build'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Max Iterations", html_for="slider"),
                        dcc.Slider(
                            id="clu-KMedoids-max-iter",
                            min=10, max=1000, step=10,
                            value=300,
                            marks={i: str(i) for i in
                                   [10] + list(range(100, 1001, 100))},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="clu-KMedoids-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=None
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn-extra.readthedocs.io/en/stable/"
                               "generated/sklearn_extra.cluster.KMedoids.html",
                               href="https://scikit-learn-extra.readthedocs.io"
                               "/en/stable/generated/sklearn_extra.cluster"
                               ".KMedoids.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="clu-KMedoids-settings",
    target="clu-KMedoids-btn",
    trigger="click"
)


clu_kmedoids_settings_keys = {
    'clu-KMedoids-n-clusters': 'n_clusters',
    'clu-KMedoids-metric': 'metric',
    'clu-KMedoids-method': 'method',
    'clu-KMedoids-init': 'init',
    'clu-KMedoids-max-iter': 'max_iter',
    'clu-KMedoids-random-state': 'random_state',
}


clu_spectral_clustering_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Spectral Clustering Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Clusters"),
                        dbc.Input(
                            id="clu-SpectralClustering-n-clusters",
                            type="text",
                            value='[4, 6, 8, 10, 12, 14, 16]'
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Eigen Solver"),
                        dbc.Select(
                            options=[
                                {"label": "arpack", "value": "arpack"},
                                {"label": "lobpcg", "value": "lobpcg"},
                                {"label": "amg", "value": "amg"}
                            ],
                            id="clu-SpectralClustering-eigen-solver",
                            value="arpack"
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("No. of Components (0 if n_clusters)",
                                  html_for="slider"),
                        dcc.Slider(
                            id="clu-SpectralClustering-n-components",
                            min=0, max=100, step=1,
                            value=0,
                            marks={i: str(i) for i in range(0, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                # dbc.FormGroup(
                #     [
                #         dbc.Label("Gamma (only for rbf, poly, "
                #                   "sigmoid, laplacian, chi2)"),
                #         dbc.Input(
                #             id="clu-SpectralClustering-gamma",
                #             type="number",
                #             value=1.0
                #         )
                #     ]
                # ),
                # dbc.FormGroup(
                #     [
                #         dbc.Label("Affinity"),
                #         dbc.Select(
                #             id="clu-SpectralClustering-affinity",
                #             options=[
                #                 {"label": "rbf", "value": "rbf"},
                #                 {"label": "nearest_neighbors",
                #                     "value": "nearest_neighbors"},
                #                 {"label": "additive_chi2",
                #                  "value": "additive_chi2"},
                #                 {"label": "chi2", "value": "chi2"},
                #                 {"label": "linear", "value": "linear"},
                #                 {"label": "poly", "value": "poly"},
                #                 {"label": "polynomial",
                #                   "value": "polynomial"},
                #                 {"label": "laplacian", "value": "laplacian"},
                #                 {"label": "sigmoid", "value": "sigmoid"},
                #                 {"label": "cosine", "value": "cosine"},
                #             ],
                #             value="rbf"
                #         )
                #     ]
                # ),
                dbc.FormGroup(
                    [
                        dbc.Label("No. of Neighbors", html_for="slider"),
                        dcc.Slider(
                            id="clu-SpectralClustering-n-neighbors",
                            min=2, max=100, step=1,
                            value=10,
                            marks={i: str(i) for i in range(10, 101, 10)},
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Assign Labels"),
                        dbc.Select(
                            id="clu-SpectralClustering-assign-labels",
                            options=[
                                {"label": "kmeans", "value": "kmeans"},
                                {"label": "discretize", "value": "discretize"}
                            ],
                            value="kmeans"
                        )
                    ]
                ),
                # dbc.FormGroup(
                #     [
                #         dbc.Label("Degree (only for poly)",
                #                   html_for="slider"),
                #         dcc.Slider(
                #             id="clu-SpectralClustering-degree",
                #             min=1, max=10, step=1,
                #             value=3,
                #             marks={i: str(i) for i in range(1, 11)},
                #             tooltip={'always_visible': True}
                #         )
                #     ]
                # ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="clu-SpectralClustering-random-state",
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
                               "sklearn.cluster.SpectralClustering.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.cluster."
                               "SpectralClustering.html",
                               target="_blank")
                    ],
                    className="small"
                )

            ]
        )
    ],
    id="clu-Spectral-Clustering-settings",
    target="clu-Spectral-Clustering-btn",
    trigger="click"
)


clu_spectral_clustering_settings_keys = {
    'clu-SpectralClustering-n-clusters': 'n_clusters',
    'clu-SpectralClustering-eigen-solver': 'eigen_solver',
    'clu-SpectralClustering-n-components': 'n_components',
    # 'clu-SpectralClustering-gamma': 'gamma',
    # 'clu-SpectralClustering-affinity': 'affinity',
    'clu-SpectralClustering-n-neighbors': 'n_neighbors',
    'clu-SpectralClustering-assign-labels': 'assign_labels',
    # 'clu-SpectralClustering-degree': 'degree',
    'clu-SpectralClustering-random-state': 'random_state'
}


clu_agglomerative_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Agglomerative Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Clusters"),
                        dbc.Input(
                            id="clu-Agglomerative-n-clusters",
                            type="text",
                            value='[4, 6, 8, 10, 12, 14, 16]'
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Affinity"),
                        dbc.Select(
                            options=[
                                {"label": "euclidean", "value": "euclidean"},
                                {"label": "manhattan", "value": "manhattan"},
                                {"label": "cosine", "value": "cosine"},
                            ],
                            id="clu-Agglomerative-affinity",
                            value='euclidean'
                        ),
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Linkage"),
                        dbc.Select(
                            options=[
                                {"label": "ward", "value": "ward"},
                                {"label": "complete", "value": "complete"},
                                {"label": "average", "value": "average"},
                                {"label": "single", "value": "single"},
                            ],
                            id="clu-Agglomerative-linkage",
                            value='ward'
                        ),
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/"
                               "sklearn.cluster.AgglomerativeClustering.html",
                               href="https://scikit-learn.org/stable/modules/"
                               "generated/sklearn.cluster.Agglomerative"
                               "Clustering.html",
                               target="_blank")
                    ],
                    className="small"
                )
            ]
        )
    ],
    id="clu-Agglomerative-settings",
    target="clu-Agglomerative-btn",
    trigger="click"
)


clu_agglomerative_settings_keys = {
    'clu-Agglomerative-n-clusters': 'n_clusters',
    'clu-Agglomerative-affinity': 'affinity',
    'clu-Agglomerative-linkage': 'linkage'
}


clu_settings = [
    clu_leiden_settings,
    clu_kmeans_settings,
    clu_kmedoids_settings,
    clu_spectral_clustering_settings,
    clu_agglomerative_settings
]

clu_settings_keys = {
    'clu-Leiden': clu_leiden_settings_keys,
    'clu-KMeans': clu_kmeans_settings_keys,
    'clu-KMedoids': clu_kmedoids_settings_keys,
    'clu-Spectral-Clustering': clu_spectral_clustering_settings_keys,
    'clu-Agglomerative': clu_agglomerative_settings_keys
}
