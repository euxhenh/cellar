import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html


ssclu_leiden_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Constrained Leiden Settings"),
        dbc.PopoverBody(
            [
                dbc.FormGroup(
                    [
                        dbc.Label("Merge Subsets"),
                        dbc.Checklist(
                            options=[],
                            id="ssclu-Leiden-main-checklist",
                            inline=True,
                            className="mb-2",
                            value=[]
                        ),
                        dbc.Checklist(
                            options=[],
                            id="ssclu-Leiden-side-checklist",
                            inline=True,
                            className="no-display mb-2",
                            value=[]
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Partition Type"),
                        dbc.Select(
                            id="ssclu-Leiden-partition-type",
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
                            id="ssclu-Leiden-n-iter",
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
                            id="ssclu-Leiden-resolution",
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
                            id="ssclu-Leiden-max-comm",
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
                            id="ssclu-Leiden-graph-method",
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
                            id="ssclu-Leiden-nneigh",
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
                            value=False,
                            id="ssclu-Leiden-weights",
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
                            value=False,
                            id="ssclu-Leiden-directed",
                            inline=True
                        )
                    ]
                ),
                dbc.FormGroup(
                    [
                        dbc.Label("Random State"),
                        dbc.Input(
                            id="ssclu-Leiden-random-state",
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
    id="ssclu-Constrained-Leiden-settings",
    target="ssclu-Constrained-Leiden-btn",
    trigger="click"
)

ssclu_leiden_settings_keys = {
    'ssclu-Leiden-main-checklist': 'main_constraints',
    'ssclu-Leiden-side-checklist': 'side_constraints',
    'ssclu-Leiden-partition-type': 'partition_type',
    'ssclu-Leiden-n-iter': 'n_iterations',
    'ssclu-Leiden-resolution': 'resolution_parameter',
    'ssclu-Leiden-max-comm': 'max_comm_size',
    'ssclu-Leiden-graph-method': 'graph_method',
    'ssclu-Leiden-nneigh': 'n_neighbors',
    'ssclu-Leiden-weights': 'use_weights',
    'ssclu-Leiden-directed': 'directed',
    'ssclu-Leiden-random-state': 'seed'
}


ssclu_settings = [
    ssclu_leiden_settings
]

ssclu_settings_keys = {
    'ssclu-Constrained-Leiden': ssclu_leiden_settings_keys
}
