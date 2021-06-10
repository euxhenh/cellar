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
                        dbc.Label("Resolution", html_for="slider"),
                        dcc.Slider(
                            id="clu-Leiden-resolution",
                            min=0.01, max=10, step=0.01,
                            value=1,
                            tooltip={'always_visible': True}
                        )
                    ]
                ),
                html.P(
                    [
                        "Details: ",
                        html.A("leidenalg.readthedocs.io/en/stable/reference.html#reference",
                               href="https://leidenalg.readthedocs.io/en/stable/reference.html#reference",
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
    'clu-Leiden-resolution': 'resolution_parameter'
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
                            value='(3,5,1)'
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Can be a single integer, or a list of comma separated values, "
                    "or a tuple specifying the range. e.g. (3,9,2) will start at 3"
                    "and end at 9 on increments of 2.",
                    placement='auto',
                    target="clu-KMeans-n-clusters",
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("Max Iter", html_for="slider"),
                        dcc.Slider(
                            id="clu-KMeans-max-iter",
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
                        dbc.Label("Tolerance"),
                        dbc.Input(
                            id="clu-KMeans-tol",
                            type="text",
                            value='1e-4'
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
                            value = 'auto'
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
                            value=None
                        )
                    ]
                ),
                
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.cluster.KMeans.html",
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
    'clu-KMeans-random-state': 'random_state',
    'clu-KMeans-max-iter': 'max_iter',
    'clu-KMeans-tol': 'tol',
    'clu-KMeans-algorithm': 'algorithm'
}





clu_KMedoids_settings = dbc.Popover(
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
                            value='(3,5,1)'
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Can be a single integer, or a list of comma separated values, "
                    "or a tuple specifying the range. e.g. (3,9,2) will start at 3"
                    "and end at 9 on increments of 2.",
                    placement='auto',
                    style={"fontFamily":"courier"},
                    target="clu-KMedoids-n-clusters",
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("Max Iter", html_for="slider"),
                        dcc.Slider(
                            id="clu-KMedoids-max-iter",
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
                        dbc.Label("Metric"),
                        dbc.Select(
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
                            id="clu-KMedoids-metric",
                            value = 'euclidean'
                        ),
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("Initialization"),
                        dbc.Select(
                            options=[
                                {"label": "random", "value": "random"},
                                {"label": "heuristic", "value": "heuristic"},
                                {"label": "k-medoids++", "value": "k-medoids++"},
                                ],
                            id="clu-KMedoids-init",
                            value = 'heuristic'
                        ),
                    ]
                ),
                dbc.Tooltip(
                    "Specify medoid initialization method. 'random' selects n_clusters"
                    "elements from the dataset. 'heuristic' picks the n_clusters points"
                    "with the smallest sum distance to every other point. 'k-medoids++'"
                    "follows an approach based on k-means++_, and in general, gives initial"
                    "medoids which are more separated than those generated by the other methods.",
                    placement='auto',
                    style={"fontFamily":"courier"},
                    target="clu-KMedoids-init",
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
                )
                
                #html.P(
                #    [
                #        "Details: ",
                #        html.A("scikit-learn-extra.readthedocs.io/en/stable/generated/sklearn_extra.cluster.KMedoids.html",
                #               href="https://scikit-learn-extra.readthedocs.io/en/stable/generated/sklearn_extra.cluster.KMedoids.html",
                #               target="_blank")
                #    ],
                #    className="small"
                #)
                
            ]
        )
    ],
    id="clu-KMedoids-settings",
    target="clu-KMedoids-btn",
    trigger="click"
)


clu_KMedoids_settings_keys = {
    'clu-KMedoids-n-clusters': 'n_clusters',
    'clu-KMedoids-random-state': 'random_state',
    'clu-KMedoids-max-iter': 'max_iter',
    'clu-KMedoids-init': 'init',
    'clu-KMedoids-metric': 'metric'
}






clu_SpectralClustering_settings = dbc.Popover(
    [
        dbc.PopoverHeader("SpectralClustering Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Number of Clusters"),
                        dbc.Input(
                            id="clu-SpectralClustering-n-clusters",
                            type="text",
                            value='(3,5,1)'
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Can be a single integer, or a list of comma separated values, "
                    "or a tuple specifying the range. e.g. (3,9,2) will start at 3"
                    "and end at 9 on increments of 2.",
                    placement='auto',
                    style={"fontFamily":"courier"},
                    target="clu-SpectralClustering-n-clusters",
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("No. of components", html_for="slider"),
                        dcc.Slider(
                            id="clu-SpectralClustering-n-components",
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
                            
                            options=[
                                {"label": "arpack", "value": "arpack"},
                                {"label": "logpcg", "value": "logpcg"},
                                {"label": "amg", "value": "amg"}
                            ],
                            id="clu-SpectralClustering-eigen-solver",
                            value="arpack"
                        )
                    ]
                ),
                
                
                                

                dbc.FormGroup(
                    [
                        dbc.Label("Gamma"),
                        dbc.Input(
                            id="clu-SpectralClustering-gamma",
                            type="text",
                            value='1.0'
                        )
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("Affinity"),
                        dbc.Select(
                            id="clu-SpectralClustering-affinity",
                            options=[
                                {"label": "nearest_neighbors", "value": "nearest_neighbors"},
                                {"label": "rbf", "value": "rbf"}
                            ],
                            value="rbf"
                        )
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("No. of neighbors"),
                        dbc.Input(
                            id="clu-SpectralClustering-n_neighbors",
                            type="number",
                            value=10
                        )
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("eigen_tol"),
                        dbc.Input(
                            id="clu-SpectralClustering-eigen_tol",
                            type="text",
                            value='0.0'
                        )
                    ]
                ),

                
                dbc.FormGroup(
                    [
                        dbc.Label("Assign Labels"),
                        dbc.Select(
                            id="clu-SpectralClustering-assign_labels",
                            options=[
                                {"label": "kmeans", "value": "kmeans"},
                                {"label": "discretize", "value": "discretize"}
                            ],
                            value="kmeans"
                        )
                    ]
                ),
                
                
                dbc.FormGroup(
                    [
                        dbc.Label("degree"),
                        dbc.Input(
                            id="clu-SpectralClustering-degree",
                            type="text",
                            value='3.0'
                        )
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("coef0"),
                        dbc.Input(
                            id="clu-SpectralClustering-coef0",
                            type="text",
                            value='1.0'
                        )
                    ]
                ),

                
                
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="clu-SpectralClustering-random-state",
                            placeholder="Leave empty if None",
                            type="number",
                            value=None
                        )
                    ]
                ),
                
                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.cluster.SpectralClustering.html",
                               target="_blank")
                    ],
                    className="small"
                )
                
            ]
        )
    ],
    id="clu-SpectralClustering-settings",
    target="clu-SpectralClustering-btn",
    trigger="click"
)


clu_SpectralClustering_settings_keys = {
    'clu-SpectralClustering-n-clusters': 'n_clusters',
    'clu-SpectralClustering-eigen-solver':'eigen_solver',
    'clu-SpectralClustering-random-state': 'random_state',
    'clu-SpectralClustering-coef0': 'coef0',
    'clu-SpectralClustering-degree': 'degree',
    'clu-SpectralClustering-assign_labels': 'assign_labels',
    'clu-SpectralClustering-eigen_tol': 'eigen_tol',
    'clu-SpectralClustering-n_neighbors': 'n_neighbors',
    'clu-SpectralClustering-affinity': 'affinity',
    'clu-SpectralClustering-gamma': 'gamma',
    'clu-SpectralClustering-n-components':'n_components'
    
    
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
                            value='(3,5,1)'
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Can be a single integer, or a list of comma separated values, "
                    "or a tuple specifying the range. e.g. (3,9,2) will start at 3"
                    "and end at 9 on increments of 2.",
                    placement='auto',
                    target="clu-Agglomerative-n-clusters",
                ),
                
                

                dbc.FormGroup(
                    [
                        dbc.Label("Affinity"),
                        dbc.Select(
                            options=[
                                {"label": "euclidean", "value": "euclidean"},
                                {"label": "l1", "value": "l1"},
                                {"label": "l2", "value": "l2"},
                                {"label": "manhattan", "value": "manhattan"},
                                {"label": "cosine", "value": "cosine"},
                                ],
                            id="clu-Agglomerative-affinity",
                            value = 'euclidean'
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
                            value = 'ward'
                        ),
                    ]
                ),
                
                dbc.FormGroup(
                    [
                        dbc.Label("Distance Threshold"),
                        dbc.Input(
                            id="clu-Agglomerative-distance_threshold",
                            placeholder="Leave empty if None",
                            type="text",
                            value=None
                        )
                    ]
                ),
                

                html.P(
                    [
                        "Details: ",
                        html.A("scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html",
                               href="https://scikit-learn.org/stable/modules/generated/sklearn.cluster.AgglomerativeClustering.html",
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
    'clu-Agglomerative-linkage': 'linkage',
    'clu-Agglomerative-distance_threshold': 'distance_threshold'
}



clu_settings = [clu_leiden_settings,clu_kmeans_settings,clu_KMedoids_settings,
                clu_SpectralClustering_settings,clu_agglomerative_settings]
clu_settings_keys = {
    'clu-Leiden': clu_leiden_settings_keys,
    'clu-KMeans': clu_kmeans_settings_keys,
    'clu-KMedoids': clu_KMedoids_settings_keys,
    'clu-SpectralClustering': clu_SpectralClustering_settings_keys,
    'clu-Agglomerative':clu_agglomerative_settings_keys
}
