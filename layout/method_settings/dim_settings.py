import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


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



dim_settings = [dim_pca_settings, dim_tsvd_settings]
dim_settings_keys = {
    'dim-PCA': dim_pca_settings_keys,
    'dim-Truncated-SVD': dim_tsvd_settings_keys
}