import os

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import dash_table

from .misc import (empty_analysis_figure, empty_colocalization_figure,
                   empty_spatial_figure)


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
                                        dbc.Tooltip(
                                            "Subset to find DE genes for",
                                            target=prefix + "-de-cluster-select",
                                            delay=500,
                                            placement='top'
                                        ),
                                        dbc.Select(
                                            options=[],
                                            id=prefix + "-de-cluster-select2",
                                            className="mr-1"
                                        ),
                                        dbc.Tooltip(
                                            "Subset to compare against",
                                            target=prefix + "-de-cluster-select2",
                                            delay=500,
                                            placement='top'
                                        ),
                                        dbc.Button(
                                            html.I(className="fas fa-cog"),
                                            id=prefix + "-de-analysis-btn",
                                            color='primary',
                                            outline=True
                                        ),
                                        dbc.Popover(
                                            [
                                                dbc.PopoverHeader(
                                                    "DE Settings"),
                                                dbc.PopoverBody([
                                                    dbc.FormGroup(
                                                        [
                                                            dbc.Label(
                                                                "Alpha"),
                                                            dbc.Input(
                                                                id=prefix + "-de-analysis-alpha",
                                                                type="number",
                                                                value=0.05
                                                            ),
                                                            dbc.Tooltip(
                                                                "q-Value Threshold",
                                                                target=prefix + "-de-analysis-alpha",
                                                                delay=500
                                                            )
                                                        ]
                                                    ),
                                                    dbc.FormGroup(
                                                        [
                                                            dbc.Label(
                                                                "FoldChange Threshold"),
                                                            dbc.Input(
                                                                id=prefix + "-de-analysis-fc",
                                                                type="number",
                                                                value=1
                                                            ),
                                                            dbc.Tooltip(
                                                                "Genes with FC less than this value will be filtered. " +
                                                                "If set to 0, no filtering will be applied.",
                                                                target=prefix + "-de-analysis-fc",
                                                                delay=500
                                                            )
                                                        ]
                                                    ),
                                                    dbc.FormGroup(
                                                        [
                                                            dbc.Label(
                                                                "Is Data Logged?"),
                                                            dbc.RadioItems(
                                                                options=[
                                                                    {"label": "True",
                                                                        "value": True},
                                                                    {"label": "False",
                                                                        "value": False},
                                                                ],
                                                                value=True,
                                                                id=prefix + "-de-analysis-logged",
                                                                inline=True
                                                            ),
                                                            dbc.Tooltip(
                                                                "If your data has been log-transformed " +
                                                                "during the preprocessing step, set this to True.",
                                                                target=prefix + "-de-analysis-logged",
                                                                delay=500
                                                            )
                                                        ]
                                                    ),
                                                ])
                                            ],
                                            target=prefix + "-de-analysis-btn",
                                            trigger="legacy"
                                        ),
                                        dbc.InputGroupAddon(
                                            dbc.Button(
                                                "Find DE Genes",
                                                id=prefix+"-find-de-genes-btn",
                                                color='primary'
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
                            className="card-title", id=prefix + "-feat-title"),
                    dbc.Row(
                        [
                            dbc.Col([
                                dbc.Button(
                                    dbc.Row(
                                        [
                                            html.Span("DE")
                                        ],
                                        align='baseline',
                                        justify='center',
                                        no_gutters=True
                                    ),
                                    id=prefix + "-paste-de-genes",
                                    color='secondary'
                                ),
                                dbc.Tooltip(
                                    "Paste DE genes from the current page",
                                    target=prefix + "-paste-de-genes",
                                    delay=500
                                )
                            ], width=1),
                            dbc.Col([
                                dcc.Loading(
                                    dcc.Dropdown(
                                        options=[],
                                        multi=True,
                                        id=prefix + '-feature-list',
                                        placeholder="Search gene by symbol",
                                    ), type="circle")
                            ], width=7, id=prefix + "-feature-list-col"),
                            dbc.Tooltip(
                                "Select one or multiple features " +
                                "for visualization purposes. This input " +
                                "will automatically parse any selected " +
                                "gene from the DE table. To understand " +
                                "what is shown when multiple genes are " +
                                "selected, please check the documentation.",
                                target=prefix + "-feat-title",
                                delay=1000
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
                                    label="Plotting",
                                    id=prefix + "-plotting-menu"
                                ),
                                width=2
                            ),
                            dbc.Col([
                                dbc.Button(
                                    "Clear",
                                    id=prefix + "-clear-expression-btn",
                                    color='primary'
                                ),
                                dbc.Tooltip(
                                    "Clears the expression view from the plot",
                                    target=prefix + "-clear-expression-btn",
                                    delay=500
                                )
                            ], width=2)
                        ],
                        justify='between',
                        no_gutters=True,
                        className="mb-2"
                    ),
                    dbc.Row(
                        dbc.Col([
                            dcc.RangeSlider(
                                min=0,
                                max=0,
                                step=None,
                                id=prefix + "-feature-rangeslider",
                                marks={}
                            ),

                        ], id=prefix + "-feature-rangeslider-col"),
                        justify='between',
                        no_gutters=True
                    ),
                    dbc.Tooltip(
                        "Clip expression values to this range. " +
                        "Can be useful to filter outliers or modify " +
                        "the colorbar scale. Only works for a " +
                        "single gene.",
                        target=prefix + "-feature-rangeslider-col",
                        delay=500
                    ),
                    dcc.Loading(
                        dcc.Graph(
                            id=prefix + '-analysis-plot',
                            figure=empty_analysis_figure,
                            config={'autosizable': True,
                                    'displaylogo': False,
                                    'displayModeBar': True}
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
                                        ),
                                        dbc.Tooltip(
                                            "Will only use this many DE " +
                                            "genes from the DE table. DE " +
                                            "genes are sorted by Fold-Change.",
                                            target=prefix + "-de-genes-enrich-no",
                                            delay=500
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
                                    block=True,
                                    color='primary'
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
                                                    'if': {
                                                        'column_id': 'Genes'},
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
                                    dcc.Loading(
                                        id=prefix + "-buf-load",
                                        fullscreen=True
                                    ),
                                    dcc.Dropdown(
                                        options=[
                                            {'label': 'CODEX',
                                             'value': 'spatial-codex'},
                                            {'label': '10X Genomics',
                                             'value': 'spatial-10x'}
                                        ],
                                        className="mw-300",
                                        placeholder="Select spatial data type",
                                        id=prefix + '-spatial-type-dropdown'
                                    ),
                                    # dbc.Tooltip(
                                    #     "Please check the documentation " +
                                    #     "for file format details.",
                                    #     target=prefix + '-spatial-type-dropdown'
                                    # ),
                                    dcc.Upload(
                                        id=prefix + "-upload-spatial",
                                        children=[
                                            dbc.Button("Upload")
                                        ],
                                        className="mr-2",
                                        max_size=150*1024*1024  # 150 MB
                                    ),
                                    dbc.Col(
                                        dcc.Loading(
                                            dcc.Dropdown(
                                                options=[],
                                                multi=True,
                                                id=prefix + '-feature-list-spatial',
                                                placeholder="Leave empty for cluster view",
                                            ),
                                            type="circle"
                                        ),
                                        width=4,
                                        id=prefix + "-feature-list-spatial-col"
                                    ),
                                    dbc.Button(
                                        "Generate Tile",
                                        color='primary',
                                        className="mr-2",
                                        id=prefix + "-generate-tile-btn"
                                    ),
                                    dbc.Tooltip(
                                        "Generate a spatial tile using " +
                                        "the uploaded file. Note: there is " +
                                        "no need to upload any supplementary " +
                                        "files for server datasets and " +
                                        "hitting this button will generate " +
                                        "the tile automatically.",
                                        target=prefix + "-generate-tile-btn",
                                        delay=500
                                    ),
                                    dcc.Loading(
                                        dcc.Download(
                                            id=prefix + "-download-tile-buf")
                                    ),
                                    dbc.Button(
                                        "Download Tile",
                                        color='primary',
                                        id=prefix + "-download-tile-btn"
                                    ),
                                    dbc.Tooltip(
                                        "This will download the tile as a " +
                                        "high resolution .png image file.",
                                        target=prefix + "-download-tile-btn",
                                        delay=500
                                    )
                                ],
                                no_gutters=True
                            ),
                            dbc.Row([
                                dbc.Col(
                                    dcc.Loading(
                                        children=[
                                            # html.Div(id=prefix + '-tile',)
                                            dcc.Graph(
                                                id=prefix + '-tile',
                                                figure=empty_spatial_figure,
                                                config={
                                                    'autosizable': True,
                                                    'displaylogo': False,
                                                    # 'staticPlot': True,
                                                    # 'toImageButtonOptions': {
                                                    #     'format': 'png'
                                                    # },
                                                    'modeBarButtonsToRemove': [
                                                        'toImage',
                                                    ]
                                                }
                                            )
                                        ],
                                        type="circle"
                                    ),
                                    align='center'
                                )
                            ], align='center', no_gutters=True)
                        ],
                        width=12
                    )
                ]
            )
        ],
        className="mb-2 analysis-container"
    )
    return spatial


def get_spatial_cluster_score_card(prefix):
    spatial_cluster_score_card = dbc.Card([
        dbc.CardBody([
            html.H5([
                "Cluster Co-Localization Score",
                dbc.Button(
                    "Compute Scores",
                    color='light',
                    className="mr-2 ml-4",
                    size='sm',
                    outline=False,
                    id=prefix + "-generate-cluster-scores-btn"
                ),
                dbc.Button(
                    html.I(className="fas fa-cog"),
                    id=prefix + "-spatial-cluster-scores-cog",
                    color='light',
                    # outline=True
                ),
                dbc.Popover(
                    [
                        dbc.PopoverHeader("Settings"),
                        dbc.PopoverBody([
                            dbc.FormGroup([
                                dbc.Label("# Neighbors"),
                                dcc.Slider(
                                    id=prefix + "-spatial-cluster-scores-nneigh",
                                    min=3, max=15, step=1,
                                    value=3,
                                    marks={i: str(i) for i in range(3, 16, 1)},
                                    tooltip={'always_visible': True}
                                )
                            ])
                        ])
                    ],
                    target=prefix + "-spatial-cluster-scores-cog",
                    trigger="legacy"
                ),
            ], className="card-title"),
            dcc.Loading(
                children=[dcc.Graph(
                    id=prefix + '-cluster-scores',
                    figure=empty_colocalization_figure,
                    config={
                        'autosizable': True,
                        'displaylogo': False,
                        'displayModeBar': True,
                        'toImageButtonOptions': {
                            'format': 'png'
                        }
                    }
                )],
                type="circle"
            )
        ], className="overflow-auto")
    ])
    return spatial_cluster_score_card


def get_spatial_protein_score_card(prefix):
    spatial_protein_score_card = dbc.Card([
        dbc.CardBody([
            html.H5([
                "Protein Co-Localization Score",
                dbc.Button(
                    "Compute Scores",
                    color='light',
                    className="mr-2 ml-4",
                    size='sm',
                    outline=False,
                    id=prefix + "-generate-protein-scores-btn"
                ),
                dbc.Button(
                    html.I(className="fas fa-cog"),
                    id=prefix + "-spatial-protein-scores-cog",
                    color='light'
                ),
                dbc.Popover(
                    [
                        dbc.PopoverHeader("Settings"),
                        dbc.PopoverBody([
                            dbc.FormGroup([
                                dbc.Label("# Neighbors"),
                                dcc.Slider(
                                    id=prefix + "-spatial-protein-scores-nneigh",
                                    min=3, max=15, step=1,
                                    value=5,
                                    marks={i: str(i) for i in range(3, 16, 1)},
                                    tooltip={'always_visible': True}
                                )
                            ])
                        ])
                    ],
                    target=prefix + "-spatial-protein-scores-cog",
                    trigger="legacy"
                ),
            ], className="card-title"),
            dcc.Loading(
                children=[dcc.Graph(
                    id=prefix + '-protein-scores',
                    figure=empty_colocalization_figure,
                    config={
                        'autosizable': True,
                        'displaylogo': False,
                        'displayModeBar': True,
                        'toImageButtonOptions': {
                            'format': 'png'
                        }
                    }
                )],
                type="circle"
            )
        ], className="overflow-auto")
    ])
    return spatial_protein_score_card


def get_spatial_scores_carddeck(prefix):
    scores = dbc.CardDeck(
        [
            get_spatial_cluster_score_card(prefix),
            get_spatial_protein_score_card(prefix)
        ],
        className="mb-2 analysis-container"
    )

    return scores


def get_asct_tables(prefix):
    tables = {
        "Brain": "ASCT-B_Allen_Brain.csv",
        "Lymph Node": "ASCT-B_NIH_Lymph_Node.csv",
        "Bone Marrow / Blood / Pelvis": "ASCT-B_VH_BM_Blood_Pelvis.csv",
        "Heart": "ASCT-B_VH_Heart.csv",
        "Large Intestine": "ASCT-B_VH_Intestine_Large.csv",
        "Kidney": "ASCT-B_VH_Kidney.csv",
        "Lung": "ASCT-B_VH_Lung.csv",
        "Skin": "ASCT-B_VH_Skin.csv",
        "Spleen": "ASCT-B_VH_Spleen.csv",
        "Thymus": "ASCT-B_VH_Thymus.csv",
        "Vasculature": "ASCT-B_VH_Vasculature.csv"
    }

    asct_table = dbc.Card([
        dbc.CardBody([
            html.H5("Anatomical Structures, Cell Types, and Biomarkers Tables",
                    className="card-title"),
            dbc.Select(
                id=prefix + "-asct-select",
                options=[{
                    "label": key, "value": value
                } for key, value in tables.items()]
            ),
            dash_table.DataTable(
                id=prefix + "-asct-table",
                page_size=10,
                export_format='csv',
                style_table={
                    'overflowX': 'auto'
                },
                style_cell={'padding': '5px'}
            ),
        ])
    ], className="mb-2 analysis-container")

    return asct_table


def get_analysis_tabs(prefix):
    analysis_tabs = dbc.Tabs(
        [
            dbc.Tab(
                dbc.Col([
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
            dbc.Tab(
                dbc.Col([
                    dbc.Row(dbc.Col(
                        get_spatial_card(prefix), width=12),
                        no_gutters=True),
                    dbc.Row(dbc.Col(
                        get_spatial_scores_carddeck(prefix), width=12),
                        no_gutters=True)
                ],
                    className="p-0",
                    width=12
                ),
                label="Spatial Data"
            ),
            dbc.Tab(get_asct_tables(prefix), label="ASCT-B")
        ]
    )

    return analysis_tabs


main_analysis_tabs = get_analysis_tabs("main")
side_analysis_tabs = get_analysis_tabs("side")
