import dash_bootstrap_components as dbc
import dash_html_components as html
import dash_core_components as dcc
import dash_table
import gseapy

from .misc import empty_analysis_figure, empty_spatial_figure


def get_de_card(prefix):
    de_card = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.H5("DE Analysis", className="card-title"),
                    dbc.Col(
                        [
                            dbc.Row(
                                dbc.InputGroup(
                                    [
                                        dbc.Select(
                                            options=[],
                                            id=prefix + "-de-cluster-select"
                                        ),
                                        dbc.InputGroupAddon(
                                            dbc.Button(
                                                "Find DE Genes",
                                                id=prefix + "-find-de-genes-btn",
                                            ),
                                            addon_type="append"
                                        )
                                    ]
                                ),
                                className="mb-2",
                                no_gutters=True
                            ),
                            dbc.Row(
                                dbc.Col(
                                    [
                                        dcc.Loading(
                                            [
                                                html.H6(
                                                    "",
                                                    id=prefix + "-de-analysis-title",
                                                    className="card-subtitle m-2"
                                                ),
                                                dash_table.DataTable(
                                                    id=prefix + "-de-table",
                                                    page_size=10,
                                                    export_format='none',
                                                    style_table={
                                                        'overflowX': 'auto'}
                                                ),
                                            ],
                                            type="circle"
                                        )
                                    ],
                                    width=12
                                ),
                                no_gutters=True
                            )
                        ]
                    )
                ]
            )
        ]
    )
    return de_card


def get_feature_card(prefix):
    feature_card = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.H5("Feature Visualization",
                            className="card-title"),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Loading(
                                    dcc.Dropdown(
                                        options=[],
                                        multi=True,
                                        id=prefix + '-feature-list',
                                        placeholder="Search gene by symbol",
                                    ), type="circle"),
                                width=8
                            ),
                            dbc.Col(
                                dbc.DropdownMenu(
                                    [
                                        dbc.DropdownMenuItem(
                                            "Plot Expression",
                                            id=prefix + "-expression"
                                        ),
                                        dbc.DropdownMenuItem(
                                            "Heatmap",
                                            id=prefix + "-heatmap"
                                        ),
                                        dbc.DropdownMenuItem(
                                            "Violin Plot",
                                            id=prefix + "-violin-plot"
                                        )
                                    ],
                                    label="Plotting"
                                ),
                                width=2
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Clear PLOT",
                                    id=prefix + "-clear-expression-btn"
                                ),
                                width=2
                            )
                        ],
                        justify='between',
                        no_gutters=True
                    ),
                    dcc.Loading(
                        dcc.Graph(
                            id=prefix + '-analysis-plot',
                            figure=empty_analysis_figure,
                            config={'autosizable': True,
                                    'displaylogo': False,
                                    'displayModeBar': False}
                        ),
                        type="circle"
                    )
                ],
                style={'overflow': 'auto'}
            )
        ]
    )

    return feature_card


def get_features_carddeck(prefix):
    features = dbc.CardDeck(
        [
            get_de_card(prefix),
            get_feature_card(prefix)
        ],
        className="mb-2 analysis-container"
    )

    return features


to_show = [
    "Cell Type",
    "GO_Biological_Process_2018",
    "GO_Cellular_Component_2018",
    "GO_Molecular_Function_2018",
    "KEGG_2021_Human",
    "KEGG_2019_Mouse",
    "MSigDB_Computational",
    "MSigDB_Hallmark_2020",
    "WikiPathway_2021_Human",
    "WikiPathways_2019_Mouse",
    "Elsevier_Pathway_Collection",
    "OMIM_Disease",
    "OMIM_Expanded",
    "TRANSFAC_and_JASPAR_PWMs",
    "Reactome_2016",
    "BioPlex_2017",
    "CORUM",
    "huMAP",
    "Panther_2016"
]


def get_enrich_card(prefix):
    enrich = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.H5("Enrichment Analysis", className="card-title"),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Dropdown(
                                    options=[
                                        {'label': i, 'value': 'gene-set-' + i}
                                        for i in to_show
                                    ],
                                    value='gene-set-' + to_show[0],
                                    id=prefix + '-gene-set-dropdown',
                                    className="mr-2",
                                ),
                                width=6
                            ),
                            dbc.Col(
                                dbc.InputGroup(
                                    [
                                        dbc.InputGroupText(
                                            "Top DE genes to use"
                                        ),
                                        dbc.Input(
                                            type="number",
                                            min=1, step=1,
                                            value=50,
                                            id=prefix + "-de-genes-enrich-no"
                                        )
                                    ],
                                    className="mr-2",
                                ),
                                width=4
                            ),
                            dbc.Col(
                                dbc.Button(
                                    "Run",
                                    id=prefix + "-run-enrich-btn",
                                    block=True
                                ),
                                width=2
                            )
                        ],
                        no_gutters=True
                    ),
                    dbc.Col(
                        dcc.Loading(
                            [
                                dbc.Row(
                                    dbc.Col(
                                        html.H6(
                                            "",
                                            id=prefix + "-enrich-title",
                                            className="card-subtitle m-2"
                                        ),
                                        width=12
                                    ),
                                    no_gutters=True
                                ),
                                dbc.Row(
                                    dbc.Col(
                                        dash_table.DataTable(
                                            id=prefix + "-enrich-table",
                                            page_size=10,
                                            export_format='none',
                                            style_table={'overflowX': 'auto'},
                                            style_cell_conditional=[
                                                {
                                                    'if': {'column_id': 'Genes'},
                                                    'textAlign': 'left'
                                                }
                                            ]
                                        ),
                                        width=12
                                    ),
                                    no_gutters=True
                                )
                            ],
                            type="circle"
                        )
                    )
                ],
                className="overflow-auto"
            )
        ],
        className="mb-2 analysis-container"
    )

    return enrich


def get_spatial_card(prefix):
    spatial = dbc.Card(
        [
            dbc.CardBody(
                [
                    html.H5("Spatial Data", className="card-title"),
                    dbc.Col(
                        [
                            dbc.Row(
                                [
                                    dcc.Upload(
                                        id=prefix + "-upload-spatial",
                                        children=[
                                            dbc.Button("Upload")
                                        ],
                                        className="mr-2"
                                    ),
                                    dbc.Button(
                                        "Generate Tile",
                                        id=prefix + "-generate-tile-btn"
                                    )
                                ],
                                no_gutters=True
                            ),
                            dbc.Row(
                                [
                                    dcc.Loading(
                                        children=[
                                            html.Div(id=prefix + '-tile')
                                        ],
                                        type="circle"
                                    )
                                ]
                            )
                        ],
                        width=12
                    )
                ]
            )
        ]
    )
    return spatial


def get_analysis_tabs(prefix):
    analysis_tabs = dbc.Tabs(
        [
            dbc.Tab(
                dbc.Col(
                    [
                        dbc.Row(dbc.Col(
                            get_features_carddeck(prefix), width=12),
                            no_gutters=True),
                        dbc.Row(dbc.Col(
                            get_enrich_card(prefix), width=12),
                            no_gutters=True)
                    ],
                    className="p-0",
                    width=12
                ),
                label="Analysis"
            ),
            dbc.Tab(get_spatial_card(prefix),
                    label="Spatial Data", style={'height': '400px'})
        ]
    )

    return analysis_tabs


main_analysis_tabs = get_analysis_tabs("main")
side_analysis_tabs = get_analysis_tabs("side")
