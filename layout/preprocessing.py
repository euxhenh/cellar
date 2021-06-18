import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .misc import empty_prep_figure

row_1 = dbc.Row([
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Filter Cells", className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Counts", html_for="slider"),
                                            dcc.RangeSlider(
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 2000],
                                                marks={i: str(i) for i in
                                                       range(0, 3001, 500)},
                                                # tooltip={
                                                #   'always_visible': True}
                                            )
                                        ]
                                    )
                                )
                            ],
                            no_gutters=True
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Genes", html_for="slider"),
                                            dcc.RangeSlider(
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 2000],
                                                marks={i: str(i) for i in
                                                       range(0, 3001, 500)},
                                                # tooltip={
                                                #   'always_visible': True}
                                            )
                                        ]
                                    )
                                )
                            ],
                            no_gutters=True
                        ),
                        dbc.Button(
                            "Run", block=True, color='secondary',
                            outline=True,
                            className="mb-3"
                        ),
                        html.P(
                            [
                                "Details: ",
                                html.A("scanpy.readthedocs.io/en/stable/api/"
                                       "scanpy.pp.filter_cells.html#scanpy.pp."
                                       "filter_cells",
                                       href="https://scanpy.readthedocs.io/en/"
                                       "stable/api/scanpy.pp.filter_cells."
                                       "html#scanpy.pp.filter_cells",
                                       target="_blank")
                            ],
                            className="small"
                        )
                    ]
                )
            ],
            className="prep-card-long"
        ),
        width=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Filter Genes", className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Counts", html_for="slider"),
                                            dcc.RangeSlider(
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 2000],
                                                marks={i: str(i) for i in
                                                       range(0, 3001, 500)},
                                                # tooltip={
                                                #     'always_visible': True}
                                            )
                                        ]
                                    )
                                )
                            ],
                            no_gutters=True
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label(
                                                "Cells", html_for="slider"),
                                            dcc.RangeSlider(
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 2000],
                                                marks={i: str(i) for i in
                                                       range(0, 3001, 500)},
                                                # tooltip={
                                                #     'always_visible': True}
                                            )
                                        ]
                                    )
                                )
                            ],
                            no_gutters=True
                        ),
                        dbc.Button(
                            "Run", block=True, color='secondary',
                            outline=True,
                            className="mb-3"
                        ),
                        html.P(
                            [
                                "Details: ",
                                html.A("scanpy.readthedocs.io/en/stable/api/"
                                       "scanpy.pp.filter_genes.html#scanpy.pp."
                                       "filter_genes",
                                       href="https://scanpy.readthedocs.io/en/"
                                       "stable/api/scanpy.pp.filter_genes."
                                       "html#scanpy.pp.filter_genes",
                                       target="_blank")
                            ],
                            className="small"
                        )
                    ]
                )
            ],
            className="prep-card-long"
        ),
        width=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Highly Variable Genes",
                                className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label("Flavor"),
                                            dbc.Select(
                                                options=[
                                                    {'label': 'seurat',
                                                        'value': 'seurat'},
                                                    {'label': 'cell_ranger',
                                                        'value': 'cell_ranger'},
                                                    {'label': 'seurat_v3',
                                                        'value': 'seurat_v3'},
                                                ],
                                                value="seurat"
                                            )
                                        ]
                                    )
                                ),
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Label("No. Top Genes"),
                                            dbc.Input(
                                                "No. top genes",
                                                type="number",
                                                placeholder="None"
                                            )
                                        ]
                                    )
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Label("Mean"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Min",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        type="number",
                                                        value=0.0125,
                                                        placeholder="None")
                                                ]
                                            )
                                        ),
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Max",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        type="number",
                                                        value=3,
                                                        placeholder="None")
                                                ]
                                            )
                                        )
                                    ],
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Label("Dispersion"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Min",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        type="number",
                                                        value=0.5,
                                                        placeholder="None")
                                                ]
                                            )
                                        ),
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Max",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        type="text",
                                                        value="inf",
                                                        placeholder="None")
                                                ]
                                            )
                                        )
                                    ],
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.Button(
                            "Run", block=True, color='secondary',
                            outline=True,
                            className="mb-3"
                        ),
                        html.P(
                            [
                                "Details: ",
                                html.A("scanpy.readthedocs.io/en/stable/api/"
                                       "scanpy.pp.highly_variable_genes.html#"
                                       "scanpy.pp.highly_variable_genes",
                                       href="https://scanpy.readthedocs.io/en/"
                                       "stable/api/scanpy.pp.highly_"
                                       "variable_genes.html#scanpy.pp."
                                       "highly_variable_genes",
                                       target="_blank")
                            ],
                            className="small"
                        )
                    ]
                )
            ],
            className="prep-card-long"
        ),
        width=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Misc", className="card-title"),
                        dbc.FormGroup(
                            [
                                dbc.Label("Normalize Total"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Target Sum", addon_type="prepend"),
                                                    dbc.Input(
                                                        type="number",
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        ),
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Max Fraction",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        type="number",
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        )
                                    ],
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Label("Log1p"),
                                dbc.RadioItems(
                                    options=[
                                        {"label": "True", "value": True},
                                        {"label": "False", "value": False},
                                    ],
                                    value=True,
                                    inline=True
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Label("Scale"),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Zero Center",
                                                        addon_type="prepend",
                                                        className="mr-2"),
                                                    dbc.RadioItems(
                                                        options=[
                                                            {"label": "True",
                                                                "value": True},
                                                            {"label": "False",
                                                                "value": False},
                                                        ],
                                                        value=True,
                                                        # inline=True
                                                    )
                                                ]
                                            )
                                        ),
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Max Val.",
                                                        addon_type="prepend",
                                                        className="mr-2"),
                                                    dbc.Input(
                                                        type="number",
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        )
                                    ],
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.Button(
                            "Run", block=True, color='secondary',
                            outline=True,
                            className="mb-3"
                        ),
                        html.P(
                            [
                                "Details: ",
                                html.A("https://scanpy.readthedocs.io/en/"
                                       "stable/api/index.html",
                                       href="https://scanpy.readthedocs.io/"
                                       "en/stable/api/index.html",
                                       target="_blank")
                            ],
                            className="small"
                        )
                    ]
                )
            ],
            className="prep-card-long"
        ),
        width=3
    )
], no_gutters=True)


atac = dbc.Row([
    dbc.Col([
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Misc", className="card-title"),
                    ]
                )
            ]
        )
    ], width=3)
], no_gutters=True)


prep = html.Div([
    row_1
], className="mb-3", id="prep-row")


prep_tabs = dbc.Tabs(
    [
        dbc.Tab(prep, label="Preprocessing", tab_id="prep"),
        dbc.Tab([], label="scATAC-seq", tab_id="atacseq"),
    ],
    id="prep-tabs"
)
