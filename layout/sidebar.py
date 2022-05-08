import dash_bootstrap_components as dbc
from dash import html, dash_table, dcc
from controller.methods import (clu_list, dim_list, lbt_list, ssclu_list,
                                vis_list, intg_list)

from .method_settings.dim_settings import dim_settings
from .method_settings.clu_settings import clu_settings
from .method_settings.vis_settings import vis_settings
from .method_settings.ssclu_settings import ssclu_settings
from .method_settings.lbt_settings import lbt_settings
from .method_settings.intg_settings import intg_settings


def get_cog_btn(id, display='none'):
    return dbc.Button(
        html.I(className="fas fa-cog"),
        id=id,
        color='primary',
        outline=True,
        style={'display': display}
    )


dataset_shape = dbc.Row([
    dbc.Col([
        dbc.InputGroup([
            dbc.InputGroupText("No. Samples"),
            dbc.InputGroupText("N/A", className="shape-text", id="no-samples"),
        ], size='sm'),
    ], width=5,),
    dbc.Col([
        dbc.InputGroup([
            dbc.InputGroupText("No. Features"),
            dbc.InputGroupText("N/A", className="shape-text", id="no-features"),
        ], size='sm'),
    ], width=5,),
    dbc.Col([
        dbc.Button(
            html.I(className="far fa-lightbulb"),
            id="run-time-btn",
            color="info",
            outline=True,
            size='sm',
        ),
        dbc.Popover(
            [
                dbc.PopoverHeader("Run-Time Estimates"),
                dbc.PopoverBody([
                    "No dataset loaded."
                ], id="run-time-popover-body")
            ],
            id="run-time-popover",
            target="run-time-btn",
            trigger="hover"
        )
    ], width=1,)
], justify='between', className="mb-2 g-0")

dim_block = dbc.Card(
    [
        dbc.Button([
            html.Span("Dimensionality Reduction"),
            html.I(className="fas fa-chevron-down ms-2"),
        ],
            id="dim-reduce-collapse-btn",
            color='light',
            outline=False,
            className="text-nowrap"
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    [
                        dbc.Form(
                            [
                                dbc.Label("Obtain Embeddings via"),
                                dbc.InputGroup(
                                    [
                                        dbc.Select(
                                            options=[
                                                {
                                                    'label': m['label'],
                                                    'value': m['value']
                                                } for m in dim_list],
                                            required=True,
                                            id="dim-methods-select",
                                            value=dim_list[0]['value']
                                        ),
                                        *([
                                            get_cog_btn(
                                                i['value'] + '-btn')
                                            for i in dim_list
                                        ]),
                                        *dim_settings
                                    ]
                                )
                            ], className="mb-3"
                        ),
                        dbc.Form(
                            [
                                dbc.Label("Obtain 2D Embeddings via"),
                                dbc.InputGroup(
                                    [
                                        dbc.Select(
                                            options=[
                                                {
                                                    'label': m['label'],
                                                    'value': m['value']
                                                } for m in vis_list],
                                            required=True,
                                            id="vis-methods-select",
                                            value=vis_list[0]['value']
                                        ),
                                        *([
                                            get_cog_btn(i['value'] + '-btn')
                                            for i in vis_list
                                        ]),
                                        *vis_settings
                                    ]
                                )
                            ], className="mb-3"
                        ),
                        html.Div(
                            dbc.Button(
                                "Run",
                                id="dim-run-btn",
                                color='primary',
                                outline=False
                            ),
                            className="d-grid gap-2",
                        ),
                    ]
                )
            ],
            is_open=True,
            id="dim-reduce-collapse"
        )
    ],
    className="shadow-sm panel mb-2"
)

clu_block = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Form(
                    [
                        dbc.Label("Method to use"),
                        dbc.InputGroup(
                            [
                                dbc.Select(
                                    options=[
                                        {
                                            'label': m['label'],
                                            'value': m['value']
                                        } for m in clu_list],
                                    required=True,
                                    id="clu-methods-select",
                                    value=clu_list[0]['value']
                                ),
                                *([
                                    get_cog_btn(
                                        i['value'] + '-btn')
                                    for i in clu_list
                                ]),
                                *clu_settings
                            ]
                        )
                    ]
                )
            ]
        )
    ]
)

ssclu_block = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Form(
                    [
                        dbc.Label("Method to use"),
                        dbc.InputGroup(
                            [
                                dbc.Select(
                                    options=[
                                        {
                                            'label': m['label'],
                                            'value': m['value']
                                        } for m in ssclu_list],
                                    required=True,
                                    id="ssclu-methods-select",
                                    value=ssclu_list[0]['value']
                                ),
                                *([
                                    get_cog_btn(
                                        i['value'] + '-btn')
                                    for i in ssclu_list
                                ]),
                                *ssclu_settings
                            ]
                        )
                    ]
                )
            ]
        )
    ]
)

lbt_block = dbc.Card(
    [
        dbc.CardBody(
            [
                dbc.Form(
                    [
                        dbc.Label("Method to use"),
                        dbc.InputGroup(
                            [
                                dbc.Select(
                                    options=[
                                        {
                                            'label': m['label'],
                                            'value': m['value']
                                        } for m in lbt_list],
                                    required=True,
                                    id="lbt-methods-select",
                                    value=lbt_list[0]['value']
                                ),
                                *([
                                    get_cog_btn(
                                        i['value'] + '-btn')
                                    for i in lbt_list
                                ]),
                                *lbt_settings
                            ]
                        )
                    ]
                )
            ]
        )
    ]
)

label_tabs = dbc.Tabs(
    [
        dbc.Tab(clu_block, label="Unsupervised", tab_id="clu"),
        dbc.Tab(ssclu_block, label="Semi-Supervised", tab_id="ssclu"),
        dbc.Tab(lbt_block, label="Transfer", tab_id="lbt"),
    ],
    id="label-tabs"
)

label_block = dbc.Card(
    [
        dbc.Button([
            html.Span("Clustering"),
            html.I(className="fas fa-chevron-down ms-2"),
        ],
            id="clustering-collapse-btn",
            color='light',
            outline=False,
            className="text-nowrap"
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    [
                        label_tabs,
                        dbc.Button(
                            "Run",
                            id="label-run-btn",
                            color='primary',
                            outline=False
                        )
                    ],
                    className="mt-2 d-grid gap-2"
                )
            ],
            is_open=True,
            id="clustering-collapse"
        )
    ],
    className="shadow-sm panel mb-2",
)

annotations_block = dbc.Card(
    [
        dbc.Button([
            html.Span("Annotations"),
            html.I(className="fas fa-chevron-down ms-2"),
        ],
            id="annotation-collapse-btn",
            color='light',
            outline=False,
            className="text-nowrap"
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    dbc.Col(
                        [
                            dbc.Row(
                                dbc.InputGroup(
                                    [
                                        dbc.InputGroupText(
                                            dbc.Select(
                                                options=[],
                                                id="main-annotation-select"
                                            ),
                                            style={'padding': '0'},
                                            id="main-annotation-addon"
                                        ),
                                        dbc.InputGroupText(
                                            dbc.Select(
                                                options=[],
                                                id="side-annotation-select"
                                            ),
                                            style={'padding': '0'},
                                            id="side-annotation-addon",
                                            className="no-display",
                                        ),
                                        dbc.Input(
                                            value="",
                                            placeholder="Annotation",
                                            type="text",
                                            id="annotation-input"
                                        ),
                                        dbc.Tooltip(
                                            "Annotate the selected cluster",
                                            target="annotation-input"
                                        ),
                                        dbc.Button(
                                            "Store",
                                            id="annotation-store-btn",
                                            color='primary'
                                        )
                                    ]
                                ),
                                className="mb-2 g-0",
                            ),
                            dbc.Row(
                                dbc.Col(
                                    dash_table.DataTable(
                                        id="main-annotation-table",
                                        page_size=10,
                                        export_format='none',
                                        style_table={
                                            'overflowY': 'auto'
                                        },
                                        data=[
                                            {
                                                "cluster_id": "N/A",
                                                "annotation": "N/A"
                                            }
                                        ],
                                        columns=[
                                            {"name": "ID", "id": "cluster_id"},
                                            {"name": "Annotation",
                                                "id": "annotation"}
                                        ]
                                    )
                                ),
                                justify='center',
                                id="main-annotation-table-row",
                                className="mb-2 g-0"
                            ),
                            dbc.Row(
                                dbc.Col(
                                    dash_table.DataTable(
                                        id="side-annotation-table",
                                        page_size=10,
                                        export_format='none',
                                        style_table={
                                            'overflowY': 'auto'
                                        },
                                        data=[
                                            {
                                                "cluster_id": "N/A",
                                                "annotation": "N/A"
                                            }
                                        ],
                                        columns=[
                                            {"name": "ID", "id": "cluster_id"},
                                            {"name": "Annotation",
                                                "id": "annotation"}
                                        ]
                                    )
                                ),
                                justify='center',
                                id="side-annotation-table-row",
                                className="no-display mb-2 g-0"
                            )
                        ]
                    )
                )
            ],
            is_open=True,
            id="annotation-collapse"
        )
    ],
    className="shadow-sm panel mb-2",
)

tools_block = dbc.Card(
    [
        dbc.Button([
            html.Span("Tools"),
            html.I(className="fas fa-chevron-down ms-2"),
        ],
            id="tools-collapse-btn",
            color='light',
            outline=False,
            className="text-nowrap"
        ),
        dbc.Collapse(
            dcc.Loading(dbc.CardBody(
                [
                    dbc.Row(
                        [
                            dbc.InputGroup(
                                [
                                    dbc.Input(
                                        value="",
                                        placeholder="Subset Name",
                                        type="text",
                                        id="subset-name-input"
                                    ),
                                    dbc.Tooltip(
                                        "Assign the selected group " +
                                        "of cells to a subset. To select " +
                                        "a group of cells click and drag " +
                                        "on the plot using the lasso tool.",
                                        target="subset-name-input"
                                    ),
                                    dbc.Button(
                                        "Store",
                                        id="subset-name-store-btn",
                                        color='primary'
                                    ),
                                ]
                            )
                        ],
                        className="mb-2 g-0",
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dcc.Loading(
                                    dbc.InputGroup([
                                        dbc.Select(
                                            options=[],
                                            placeholder="Merge subset",
                                            id="main-subset-select"
                                        ),
                                        dbc.Select(
                                            options=[],
                                            placeholder="Merge subset",
                                            className="no-display",
                                            id="side-subset-select"
                                        ),
                                        dbc.Tooltip(
                                            "Cells belonging to this subset " +
                                            "will merge into a single " +
                                            "cluster. This is useful only " +
                                            "for user defined subsets.",
                                            target="main-subset-select"
                                        ),
                                        dbc.Tooltip(
                                            "Cells belonging to this subset " +
                                            "will merge into a single " +
                                            "cluster. This is useful only " +
                                            "for user defined subsets.",
                                            target="side-subset-select"
                                        ),
                                        dbc.Button(
                                            "Merge",
                                            id="merge-subset-btn",
                                            color='primary'
                                        )
                                    ])
                                )
                            )
                        ],
                        className="g-0",
                    ),
                    dbc.Form([
                        dbc.InputGroup([
                            dbc.InputGroupText([
                                "Integrate"
                            ]),
                            dbc.Select(
                                options=[
                                    {
                                        'label': m['label'],
                                        'value': m['value']
                                    } for m in intg_list],
                                required=True,
                                id="intg-methods-select",
                                value=intg_list[0]['value']
                            ),
                            *([
                                get_cog_btn(
                                    i['value'] + '-btn')
                                for i in intg_list
                            ]),
                            dbc.Button(
                                "Run", id="run-integrate-btn",
                                color="primary"),
                            html.Span(id="main-buf-integrate"),
                            html.Span(id="side-buf-integrate"),
                            *intg_settings
                        ])
                    ], className="mt-2")
                ]
            )),
            is_open=True,
            id="tools-collapse"
        )
    ],
    className="shadow-sm panel mb-2",
)

session_block = dbc.Card(
    [
        dbc.Button([
            html.Span("Session"),
            html.I(className="fas fa-chevron-down ms-2"),
        ],
            id="session-collapse-btn",
            color='light',
            outline=False,
            className="text-nowrap"
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    [
                        html.Div(dbc.Button(
                            "Export Session",
                            id='export-session-btn',
                            color='primary',
                            outline=False
                        ), className="d-grid gap-2 mb-2"),
                        dbc.Tooltip(
                            "Export the dataset and the analysis as an " +
                            ".h5ad file. This make take a few seconds.",
                            target="export-session-btn"
                        ),
                        html.Div(dbc.Button(
                            "Export Annotations",
                            id='export-annotations-btn',
                            color='primary',
                            outline=False
                        ), className="d-grid gap-2"),
                        dbc.Tooltip(
                            "Export the cell IDs and their annotations as " +
                            "a .csv file.",
                            target="export-annotations-btn"
                        ),
                        dcc.Loading(
                            dcc.Download("export-session-d"),
                            fullscreen=True
                        ),
                        dcc.Loading(
                            dcc.Download("export-annotations-d"),
                            fullscreen=True
                        )
                    ]
                )
            ],
            is_open=True,
            id="session-collapse"
        )
    ],
    className="shadow-sm panel mb-2",
)


sidebar = [
    dataset_shape,
    dim_block,
    label_block,
    annotations_block,
    tools_block,
    session_block
]
