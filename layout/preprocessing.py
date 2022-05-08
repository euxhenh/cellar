import dash_bootstrap_components as dbc
from dash import html, dcc

row_1 = dbc.Row([
    dbc.Col(
        dbc.Card([
            dbc.CardBody([
                dcc.Loading(
                    html.H5(
                        "Filter Cells",
                        className="card-title",
                        id="temp-prep-h5"
                    ),
                    fullscreen=True
                ),
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            dbc.Col(dbc.Checkbox(
                                "prep-filter-cells-counts-checkbox",
                                className="mr-3 prep-ckb",
                            ), width='auto'),
                            dbc.Col(dbc.Label(
                                "Counts", html_for="slider"
                            ), width='auto'),
                        ], align="baseline", className="g-0"),
                        dcc.RangeSlider(
                            id="prep-filter-cells-counts-slider",
                            min=0,
                            max=10000,
                            step=50,
                            value=[100, 3000],
                            marks={i: str(i) for i in range(0, 10001, 2000)},
                            tooltip={
                                'always_visible': True,
                                'placement': 'bottom'
                            }
                        )
                   ])
                ], className="g-0 mb-2"),
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            dbc.Col(dbc.Checkbox(
                                "prep-filter-cells-genes-checkbox",
                                className="mr-3 prep-ckb",
                                value=True
                            ), width='auto'),
                            dbc.Col(dbc.Label("Genes", html_for="slider"), width='auto'),
                        ], align="baseline", className="g-0"),
                        dcc.RangeSlider(
                            id="prep-filter-cells-genes-slider",
                            min=0,
                            max=10000,
                            step=50,
                            value=[100, 3000],
                            marks={i: str(i) for i in range(0, 10001, 2000)},
                            tooltip={
                                'always_visible': True,
                                'placement': 'bottom'
                            }
                        )
                    ])
                ], className="g-0 mb-2"),
            ])
        ], className="prep-card-long"),
        sm=12, md=6, lg=3
    ),
    dbc.Col(
        dbc.Card([
            dbc.CardBody([
                html.H5("Filter Genes", className="card-title"),
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            dbc.Col(dbc.Checkbox(
                                "prep-filter-genes-counts-checkbox",
                                className="mr-3 prep-ckb"
                            ), width='auto'),
                            dbc.Col(dbc.Label("Counts", html_for="slider"), width='auto'),
                        ], align="baseline", className="g-0"),
                        dcc.RangeSlider(
                            id="prep-filter-genes-counts-slider",
                            min=0,
                            max=10000,
                            step=50,
                            value=[100, 3000],
                            marks={i: str(i) for i in range(0, 10001, 2000)},
                            tooltip={
                                'always_visible': True,
                                'placement': 'bottom'
                            }
                        )
                    ])
                ], className="g-0 mb-2"),
                dbc.Row([
                    dbc.Col([
                        dbc.Row([
                            dbc.Col(dbc.Checkbox(
                                "prep-filter-genes-cells-checkbox",
                                className="mr-3 prep-ckb",
                                value=True
                            ), width='auto'),
                            dbc.Col(dbc.Label("Cells", html_for="slider"), width='auto'),
                        ], align="baseline", className="g-0",),
                        dcc.RangeSlider(
                            id="prep-filter-genes-cells-slider",
                            min=0,
                            max=10000,
                            step=50,
                            value=[100, 3000],
                            marks={i: str(i) for i in range(0, 10001, 2000)},
                            tooltip={
                                'always_visible': True,
                                'placement': 'bottom'
                            }
                        )
                    ]
                    )
                ], className="g-0"),
            ])
        ], className="prep-card-long"
        ),
        sm=12, md=6, lg=3
    ),
    dbc.Col(
        dbc.Card([
            dbc.CardBody([
                dbc.Row([
                    dbc.Col(dbc.Checkbox(
                        "prep-high-var-checkbox",
                        className="mr-3 prep-ckb",
                        value=False
                    ), width='auto'),
                    dbc.Col(html.H5("Highly Variable Genes", className="card-title"), width='auto')
                ], align="baseline", className="g-0"),
                dbc.Row([
                    dbc.Col([
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
                            id="prep-high-var-flavor",
                            value="seurat"
                        )
                    ]),
                    dbc.Col(
                        dbc.Form([
                            dbc.Label("No. Top Genes"),
                            dbc.Input(
                                id="prep-high-var-top-genes",
                                type="number",
                                placeholder="None"
                            )
                        ])
                    )
                ], className="g-0 mb-2"),
                dbc.Form(
                    [
                        dbc.Label("Mean"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.InputGroup(
                                        [
                                            dbc.InputGroupText("Min"),
                                            dbc.Input(
                                                type="number",
                                                value=0.0125,
                                                placeholder="None",
                                                id="prep-high-var-mean-min"
                                            )
                                        ]
                                    )
                                ),
                                dbc.Col(
                                    dbc.InputGroup(
                                        [
                                            dbc.InputGroupText("Max"),
                                            dbc.Input(
                                                type="number",
                                                value=3,
                                                placeholder="None",
                                                id="prep-high-var-mean-max"
                                            )
                                        ]
                                    )
                                )
                            ],
                            className="g-0",
                        )
                    ], className="mb-2"
                ),
                dbc.Form(
                    [
                        dbc.Label("Dispersion"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.InputGroup(
                                        [
                                            dbc.InputGroupText("Min"),
                                            dbc.Input(
                                                id="prep-high-var-disp-min",
                                                type="number",
                                                value=0.5,
                                                placeholder="None")
                                        ]
                                    )
                                ),
                                dbc.Col(
                                    dbc.InputGroup(
                                        [
                                            dbc.InputGroupText("Max"),
                                            dbc.Input(
                                                id="prep-high-var-disp-max",
                                                type="text",
                                                value="inf",
                                                placeholder="None")
                                        ]
                                    )
                                )
                            ],
                            className="g-0",
                        )
                    ]
                ),
            ]
            )
        ],
            className="prep-card-long"
        ),
        sm=12, md=6, lg=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Misc", className="card-title"),
                        dbc.Form(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(dbc.Checkbox(
                                            "prep-normalize-total-checkbox",
                                            className="mr-3 prep-ckb",
                                            value=True
                                        ), width='auto'),
                                        dbc.Col(dbc.Label("Normalize Total"), width='auto'),
                                    ],
                                    align="baseline",
                                    className="g-0",
                                ),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText(
                                                        "Target Sum"),
                                                    dbc.Input(
                                                        id="prep-norm-target",
                                                        type="number",
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        ),
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText(
                                                        "Max Fraction"),
                                                    dbc.Input(
                                                        id="prep-norm-maxf",
                                                        type="number",
                                                        value=0.05,
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        )
                                    ],
                                    className="g-0",
                                )
                            ], className="mb-2"
                        ),
                        dbc.Form(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(dbc.Checkbox(
                                            "prep-log1p-checkbox",
                                            className="mr-3 prep-ckb",
                                            value=True
                                        ), width='auto'),
                                        dbc.Col(dbc.Label("Log1p"), width='auto'),
                                    ],
                                    # justify="between",
                                    align="baseline",
                                    className="g-0",
                                )
                            ], className="mb-2"
                        ),
                        dbc.Form(
                            [
                                dbc.Row(
                                    [
                                        dbc.Col(dbc.Checkbox(
                                            "prep-scale-checkbox",
                                            className="mr-3 prep-ckb",
                                            value=True
                                        ), width='auto'),
                                        dbc.Col(dbc.Label("Scale"), width='auto'),
                                    ],
                                    # justify="between",
                                    align="baseline",
                                    className="g-0",
                                ),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupText(
                                                        "Zero Center",
                                                        className="me-2"),
                                                    dbc.RadioItems(
                                                        id="prep-scale-zero",
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
                                                    dbc.InputGroupText(
                                                        "Max Val.",
                                                        className="mr-2"),
                                                    dbc.Input(
                                                        id="prep-scale-max",
                                                        type="number",
                                                        placeholder="None"
                                                    )
                                                ]
                                            )
                                        )
                                    ],
                                    className="g-0",
                                )
                            ], className="mb-2"
                        ),
                        html.P(
                            [
                                "Details: ",
                                html.A("https://scanpy.readthedocs.io/en/"
                                       "stable/api.html#module-scanpy.pp",
                                       href="https://scanpy.readthedocs.io/"
                                       "en/stable/api.html#module-scanpy.pp",
                                       target="_blank")
                            ],
                            className="small"
                        )
                    ]
                )
            ],
            className="prep-card-long"
        ),
        sm=12, md=6, lg=3
    )
], className="g-0")


atac_row = dbc.Row([
    dbc.Col([
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        dcc.Loading(
                            html.H5(
                                "Bin Operation",
                                className="card-title",
                                id="temp-prep-h5-atac"
                            ),
                            fullscreen=True
                        ),
                        dbc.RadioItems(
                            options=[
                                {"label": "Sum", "value": "sum"},
                                {"label": "Mean", "value": "mean"},
                            ],
                            value="sum",
                            inline=True,
                            className="mb-3",
                            id="atac-operation"
                        ),
                        html.H5("Extend Gene Length", className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Extend"),
                                        dbc.Input(
                                            "atac-extend-input",
                                            type="text", value="5x"
                                        )
                                    ]
                                ), width=6),
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Max Extend"),
                                        dbc.Input(
                                            "atac-max-extend-input",
                                            type="text", value="5000"
                                        )
                                    ]
                                ), width=6, className="mb-3")
                            ]
                        ),
                        html.H5("Orientation", className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.Label("Consider Strand Orientation"),
                                    width=6
                                ),
                                dbc.Col(
                                    dbc.RadioItems(
                                        options=[
                                            {"label": "True", "value": True},
                                            {"label": "False", "value": False},
                                        ],
                                        value=False,
                                        inline=True,
                                        id="atac-strand-bool"
                                    ),
                                    width=6
                                )
                            ]
                        ),
                        dbc.Row(
                            [
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Op Extend"),
                                        dbc.Input(
                                            "atac-op-extend-input",
                                            type="text", value="1x"
                                        )
                                    ]
                                ), width=6),
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Op Max Extend"),
                                        dbc.Input(
                                            "atac-max-op-extend-input",
                                            type="text", value="1000"
                                        )
                                    ]
                                ), width=6)
                            ], className="mb-2"
                        ),
                        html.P("Note: This process may take a long time (10-30 mins).")
                    ]
                )
            ],
            className="prep-card-long"
        )
    ], width=3, xs=12, sm=12, md=3, lg=3)
], className="g-0")


cite_seq_row = dbc.Row([
    dbc.Col(
        dbc.Card([
            dbc.CardBody([
                dcc.Loading(
                    html.H5(
                        "CITE-seq Protein expression",
                        className="card-title",
                        id="temp-prep-h5-cite")
                ),
                html.H6("Must have a 'protein.X' key under adata.obsm!"),
                dbc.Form([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-normalize-total-checkbox",
                            className="mr-3 prep-ckb",
                            value=True
                        ),
                        dbc.Label("Normalize Total"),
                    ], align="baseline", className="g-0"),
                    dbc.Row([
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupText("Target Sum"),
                                dbc.Input(
                                    id="cite-prep-norm-target",
                                    type="number",
                                    placeholder="None"
                                )
                            ])
                        ),
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupText("Max Fraction"),
                                dbc.Input(
                                    id="cite-prep-norm-maxf",
                                    type="number",
                                    value=0.05,
                                    placeholder="None"
                                )
                            ])
                        )
                    ], className="g-0")
                ]),
                dbc.Form([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-log1p-checkbox",
                            className="mr-3 prep-ckb",
                            value=True
                        ),
                        dbc.Label("Log1p"),
                    ], align="baseline", className="g-0")
                ]),
                dbc.Form([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-scale-checkbox",
                            className="mr-3 prep-ckb",
                            value=True
                        ),
                        dbc.Label("Scale"),
                    ], align="baseline", className="g-0"),
                    dbc.Row([
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupText(
                                    "Zero Center",
                                    className="mr-2"),
                                dbc.RadioItems(
                                    id="cite-prep-scale-zero",
                                    options=[
                                        {"label": "True",
                                            "value": True},
                                        {"label": "False",
                                            "value": False},
                                    ],
                                    value=True
                                )
                            ])
                        ),
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupText(
                                    "Max Val.",
                                    className="mr-2"),
                                dbc.Input(
                                    id="cite-prep-scale-max",
                                    type="number",
                                    placeholder="None"
                                )
                            ])
                        )
                    ], className="g-0",)
                ]),
                html.P([
                    "Details: ",
                    html.A("https://scanpy.readthedocs.io/en/"
                           "stable/api/index.html",
                           href="https://scanpy.readthedocs.io/"
                           "en/stable/api/index.html",
                           target="_blank")
                ], className="small")
            ])
        ], className="prep-card-long"),
        width=3, xs=12, sm=12, md=3, lg=3
    )
])


prep = html.Div([
    row_1,
    dbc.Row(
        dbc.Col(
            dbc.Button(
                "Run",
                outline=False,
                color='primary',
                id="prep-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2 d-grid gap-2"
    )
], className="mb-3", id="prep-row")


atac = html.Div([
    atac_row,
    dbc.Row(
        dbc.Col(
            dbc.Button(
                "Run",
                outline=False,
                color='primary',
                id="prep-atac-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2 d-grid gap-2"
    )
], className="mb-3", id="prep-atac-row")


cite = html.Div([
    cite_seq_row,
    dbc.Row(
        dbc.Col(
            dbc.Button(
                "Run",
                outline=False,
                color='primary',
                id="prep-cite-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2 d-grid gap-2"
    )
], className="mb-3", id="prep-cite-row")


prep_tabs = dbc.Tabs(
    [
        dbc.Tab(prep, label="Preprocessing", tab_id="prep"),
        dbc.Tab(atac, label="scATAC-seq", tab_id="atacseq")
        # dbc.Tab(cite, label="CITE-seq", tab_id="citeseq"),
    ],
    id="prep-tabs"
)
