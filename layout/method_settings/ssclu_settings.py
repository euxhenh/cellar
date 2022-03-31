import dash_bootstrap_components as dbc
from dash import dcc, html


ssclu_leiden_settings = dbc.Popover(
    [
        dbc.PopoverHeader("Constrained Leiden Settings"),
        dbc.PopoverBody(
            [
                dbc.Form(
                    [
                        dbc.Label("Constrained Subsets"),
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
                    ],
                    id="ssclu-Leiden-constrained-subsets-form"
                ),
                dbc.Tooltip(
                    "The membership for the points belonging to a "
                    "constrained subset will be set to the same value.",
                    target="ssclu-Leiden-constrained-subsets-form"
                ),
                dbc.Form(
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
                dbc.Tooltip(
                    "Check Leidenalg docs for a description of each.",
                    target="ssclu-Leiden-partition-type"
                ),
                dbc.Form(
                    [
                        dbc.Label("Number of Iterations"),
                        dbc.Input(
                            id="ssclu-Leiden-n-iter",
                            placeholder="Set to -1 to run until convergence",
                            type="number",
                            value=4
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Maximum number of iterations for the leiden algorithm. "
                    "If set to -1, will run until convergence.",
                    target="ssclu-Leiden-n-iter"
                ),
                dbc.Form(
                    [
                        dbc.Label(
                            "Resolution (only for RBConfiguration, RBER, CPM)",
                            html_for="slider",
                            id="ssclu-Leiden-tooltip-temp-resolution"),
                        dcc.Slider(
                            id="ssclu-Leiden-resolution",
                            min=0.01, max=10, step=0.01,
                            value=1,
                            tooltip={'always_visible': True},
                            marks=None,
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Higher resolution results in more clusters.",
                    target="ssclu-Leiden-tooltip-temp-resolution"
                ),
                # dbc.Form(
                #     [
                #         dbc.Label("Max Community Size"),
                #         dbc.Input(
                #             id="ssclu-Leiden-max-comm",
                #             placeholder="Set to 0 to allow any size",
                #             type="number",
                #             value=0
                #         )
                #     ]
                # ),
                # dbc.Tooltip(
                #     "Maximum number of nodes in a community. Set to 0 "
                #     "to allow any size.",
                #     target="ssclu-Leiden-max-comm"
                # ),
                dbc.Form(
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
                dbc.Tooltip(
                    "Method to use for finding KNN. 'auto' uses approximate "
                    "neighbors if n_samples > 5000, otherwise computes the "
                    "exact KNN.",
                    target="ssclu-Leiden-graph-method"
                ),
                dbc.Form(
                    [
                        dbc.Label(
                            "Number of Neighbors for the Graph",
                            html_for="slider",
                            id="ssclu-Leiden-tooltip-temp-nneigh"),
                        dcc.Slider(
                            id="ssclu-Leiden-nneigh",
                            min=1, max=100, step=1,
                            value=15,
                            tooltip={'always_visible': True},
                            marks=None,
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Number of k-Nearest Neighbors to compute.",
                    target="ssclu-Leiden-tooltip-temp-nneigh"
                ),
                dbc.Form(
                    [
                        dbc.Label("Fix Membership for Constraints"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="ssclu-Leiden-fix-membership",
                            inline=True
                        )
                    ]
                ),
                dbc.Tooltip(
                    "Fix labels for all the points in the "
                    "constrained subsets.",
                    target="ssclu-Leiden-fix-membership"
                ),
                dbc.Form(
                    [
                        dbc.Label("Use Cached Neighbors"),
                        dbc.RadioItems(
                            options=[
                                {"label": "True", "value": True},
                                {"label": "False", "value": False},
                            ],
                            value=True,
                            id="ssclu-Leiden-cached-neigh",
                            inline=True
                        )
                    ]
                ),
                dbc.Tooltip(
                    "If the KNN were already computed using the same settings, "
                    "will not recompute them if True. It only makes sense to "
                    "set this to False when using approximate neighbors and "
                    "a new random seed is desired.",
                    target="ssclu-Leiden-cached-neigh"
                ),
                # dbc.Form(
                #     [
                #         dbc.Label("Construct Weighted Graph"),
                #         dbc.RadioItems(
                #             options=[
                #                 {"label": "True", "value": True},
                #                 {"label": "False", "value": False},
                #             ],
                #             value=False,
                #             id="ssclu-Leiden-weights",
                #             inline=True
                #         )
                #     ]
                # ),
                dbc.Form(
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
                dbc.Tooltip(
                    "If set to True, then the neighbors graph will be directed",
                    target="ssclu-Leiden-directed"
                ),
                # dbc.Form(
                #     [
                #         dbc.Label("Random State"),
                #         dbc.Input(
                #             id="ssclu-Leiden-random-state",
                #             placeholder="Leave empty if None",
                #             type="number",
                #             value=""
                #         )
                #     ]
                # ),
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
    trigger="legacy"
)

ssclu_leiden_settings_keys = {
    'ssclu-Leiden-main-checklist': 'main_constraints',
    'ssclu-Leiden-side-checklist': 'side_constraints',
    'ssclu-Leiden-partition-type': 'partition_type',
    'ssclu-Leiden-n-iter': 'n_iterations',
    'ssclu-Leiden-resolution': 'resolution_parameter',
    # 'ssclu-Leiden-max-comm': 'max_comm_size',
    'ssclu-Leiden-graph-method': 'graph_method',
    'ssclu-Leiden-nneigh': 'n_neighbors',
    'ssclu-Leiden-fix-membership': 'fix_membership',
    'ssclu-Leiden-cached-neigh': 'use_cached_neigh',
    # 'ssclu-Leiden-weights': 'use_weights',
    'ssclu-Leiden-directed': 'directed'
    # 'ssclu-Leiden-random-state': 'seed'
}


ssclu_settings = [
    ssclu_leiden_settings
]

ssclu_settings_keys = {
    'ssclu-Constrained-Leiden': ssclu_leiden_settings_keys
}
