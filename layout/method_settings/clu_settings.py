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
                                {"label": "SignificanceVertexPartition",
                                    "value": "SignificanceVertexPartition"},
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


clu_settings = [clu_leiden_settings]
clu_settings_keys = {
    'clu-Leiden': clu_leiden_settings_keys,
}
