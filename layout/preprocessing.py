import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

row_1 = dbc.Row([
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        dcc.Loading(
                            html.H5(
                                "Filter Cells",
                                className="card-title",
                                id="temp-prep-h5"
                            ),
                            fullscreen=True
                        ),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Checkbox(
                                                        "prep-filter-cells-counts-checkbox",
                                                        className="mr-3 prep-ckb"
                                                    ),
                                                    dbc.Label(
                                                        "Counts", html_for="slider"),
                                                ],
                                                # justify="between",
                                                align="baseline",
                                                no_gutters=True
                                            ),
                                            # dbc.Label(
                                            # "Counts", html_for="slider"),
                                            dcc.RangeSlider(
                                                id="prep-filter-cells-counts-slider",
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 3000],
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
                                            dbc.Row(
                                                [
                                                    dbc.Checkbox(
                                                        "prep-filter-cells-genes-checkbox",
                                                        className="mr-3 prep-ckb",
                                                        checked=True
                                                    ),
                                                    dbc.Label(
                                                        "Genes", html_for="slider"),
                                                ],
                                                # justify="between",
                                                align="baseline",
                                                no_gutters=True
                                            ),
                                            dcc.RangeSlider(
                                                id="prep-filter-cells-genes-slider",
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 3000],
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
        width=3, xs=12, sm=12, md=3, lg=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Filter Genes",
                                className="card-title"),
                        dbc.Row(
                            [
                                dbc.Col(
                                    dbc.FormGroup(
                                        [
                                            dbc.Row(
                                                [
                                                    dbc.Checkbox(
                                                        "prep-filter-genes-counts-checkbox",
                                                        className="mr-3 prep-ckb"
                                                    ),
                                                    dbc.Label(
                                                        "Counts", html_for="slider"),
                                                ],
                                                # justify="between",
                                                align="baseline",
                                                no_gutters=True
                                            ),
                                            dcc.RangeSlider(
                                                id="prep-filter-genes-counts-slider",
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 3000],
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
                                            dbc.Row(
                                                [
                                                    dbc.Checkbox(
                                                        "prep-filter-genes-cells-checkbox",
                                                        className="mr-3 prep-ckb",
                                                        checked=True
                                                    ),
                                                    dbc.Label(
                                                        "Cells", html_for="slider"),
                                                ],
                                                # justify="between",
                                                align="baseline",
                                                no_gutters=True
                                            ),
                                            dcc.RangeSlider(
                                                id="prep-filter-genes-cells-slider",
                                                min=0,
                                                max=3000,
                                                step=50,
                                                value=[100, 3000],
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
        width=3, xs=12, sm=12, md=3, lg=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        dbc.Row(
                            [
                                dbc.Checkbox(
                                    "prep-high-var-checkbox",
                                    className="mr-3 prep-ckb",
                                    checked=False
                                ),
                                html.H5("Highly Variable Genes",
                                        className="card-title")
                            ],
                            # justify="between",
                            align="baseline",
                            no_gutters=True
                        ),
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
                                                id="prep-high-var-flavor",
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
                                                id="prep-high-var-top-genes",
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
                                                        placeholder="None",
                                                        id="prep-high-var-mean-min"
                                                    )
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
                                                        placeholder="None",
                                                        id="prep-high-var-mean-max"
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
                                                    dbc.InputGroupAddon(
                                                        "Max",
                                                        addon_type="prepend"),
                                                    dbc.Input(
                                                        id="prep-high-var-disp-max",
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
        width=3, xs=12, sm=12, md=3, lg=3
    ),
    dbc.Col(
        dbc.Card(
            [
                dbc.CardBody(
                    [
                        html.H5("Misc", className="card-title"),
                        dbc.FormGroup(
                            [
                                dbc.Row(
                                    [
                                        dbc.Checkbox(
                                            "prep-normalize-total-checkbox",
                                            className="mr-3 prep-ckb",
                                            checked=True
                                        ),
                                        dbc.Label(
                                            "Normalize Total"),
                                    ],
                                    # justify="between",
                                    align="baseline",
                                    no_gutters=True
                                ),
                                dbc.Row(
                                    [
                                        dbc.Col(
                                            dbc.InputGroup(
                                                [
                                                    dbc.InputGroupAddon(
                                                        "Target Sum", addon_type="prepend"),
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
                                                    dbc.InputGroupAddon(
                                                        "Max Fraction",
                                                        addon_type="prepend"),
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
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Row(
                                    [
                                        dbc.Checkbox(
                                            "prep-log1p-checkbox",
                                            className="mr-3 prep-ckb",
                                            checked=True
                                        ),
                                        dbc.Label("Log1p"),
                                    ],
                                    # justify="between",
                                    align="baseline",
                                    no_gutters=True
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Row(
                                    [
                                        dbc.Checkbox(
                                            "prep-scale-checkbox",
                                            className="mr-3 prep-ckb",
                                            checked=True
                                        ),
                                        dbc.Label("Scale"),
                                    ],
                                    # justify="between",
                                    align="baseline",
                                    no_gutters=True
                                ),
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
                                                    dbc.InputGroupAddon(
                                                        "Max Val.",
                                                        addon_type="prepend",
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
                                    no_gutters=True
                                )
                            ]
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
        width=3, xs=12, sm=12, md=3, lg=3
    )
], no_gutters=True)


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
                                        dbc.InputGroupAddon(
                                            "Extend", addon_type="append"),
                                        dbc.Input(
                                            "atac-extend-input",
                                            type="text", value="5x"
                                        )
                                    ]
                                ), width=6),
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupAddon(
                                            "Max Extend", addon_type="append"),
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
                                        dbc.InputGroupAddon(
                                            "Op Extend", addon_type="append"),
                                        dbc.Input(
                                            "atac-op-extend-input",
                                            type="text", value="1x"
                                        )
                                    ]
                                ), width=6),
                                dbc.Col(dbc.InputGroup(
                                    [
                                        dbc.InputGroupAddon(
                                            "Op Max Extend",
                                            addon_type="append"),
                                        dbc.Input(
                                            "atac-max-op-extend-input",
                                            type="text", value="1000"
                                        )
                                    ]
                                ), width=6)
                            ]
                        )
                    ]
                )
            ],
            className="prep-card-long"
        )
    ], width=3, xs=12, sm=12, md=3, lg=3)
], no_gutters=True)


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
                dbc.FormGroup([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-normalize-total-checkbox",
                            className="mr-3 prep-ckb",
                            checked=True
                        ),
                        dbc.Label("Normalize Total"),
                    ], align="baseline", no_gutters=True),
                    dbc.Row([
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupAddon(
                                    "Target Sum", addon_type="prepend"),
                                dbc.Input(
                                    id="cite-prep-norm-target",
                                    type="number",
                                    placeholder="None"
                                )
                            ])
                        ),
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupAddon(
                                    "Max Fraction",
                                    addon_type="prepend"),
                                dbc.Input(
                                    id="cite-prep-norm-maxf",
                                    type="number",
                                    value=0.05,
                                    placeholder="None"
                                )
                            ])
                        )
                    ], no_gutters=True)
                ]),
                dbc.FormGroup([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-log1p-checkbox",
                            className="mr-3 prep-ckb",
                            checked=True
                        ),
                        dbc.Label("Log1p"),
                    ], align="baseline", no_gutters=True)
                ]),
                dbc.FormGroup([
                    dbc.Row([
                        dbc.Checkbox(
                            "cite-prep-scale-checkbox",
                            className="mr-3 prep-ckb",
                            checked=True
                        ),
                        dbc.Label("Scale"),
                    ], align="baseline", no_gutters=True),
                    dbc.Row([
                        dbc.Col(
                            dbc.InputGroup([
                                dbc.InputGroupAddon(
                                    "Zero Center",
                                    addon_type="prepend",
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
                                dbc.InputGroupAddon(
                                    "Max Val.",
                                    addon_type="prepend",
                                    className="mr-2"),
                                dbc.Input(
                                    id="cite-prep-scale-max",
                                    type="number",
                                    placeholder="None"
                                )
                            ])
                        )
                    ], no_gutters=True)
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
                block=True,
                id="prep-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2"
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
                block=True,
                id="prep-atac-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2"
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
                block=True,
                id="prep-cite-run-btn"
            ),
            width=3
        ),
        justify='center',
        className="mt-2"
    )
], className="mb-3", id="prep-cite-row")


prep_tabs = dbc.Tabs(
    [
        dbc.Tab(prep, label="Preprocessing", tab_id="prep"),
        dbc.Tab(atac, label="scATAC-seq", tab_id="atacseq"),
        dbc.Tab(cite, label="CITE-seq", tab_id="citeseq"),
    ],
    id="prep-tabs"
)
