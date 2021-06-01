import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from controller.methods import (clu_list, dim_list, lbt_list, ssclu_list,
                                vis_list)

from .method_settings.dim_settings import dim_settings
from .method_settings.clu_settings import clu_settings


def get_cog_btn(id, display='none'):
    return dbc.Button(
        html.I(className="fas fa-cog"),
        id=id,
        color='primary',
        outline=True,
        style={'display': display}
    )


dataset_shape = dbc.Row(
    [
        dbc.Col(
            [
                dbc.InputGroup(
                    [
                        dbc.InputGroupAddon(
                            "No. Samples", addon_type="prepend"),
                        dbc.InputGroupText(
                            "N/A", className="shape-text", id="no-samples"),
                    ]
                )
            ]
        ),
        dbc.Col(
            [
                dbc.InputGroup(
                    [
                        dbc.InputGroupAddon(
                            "No. Features", addon_type="prepend"),
                        dbc.InputGroupText(
                            "N/A", className="shape-text", id="no-features"),
                    ]
                )
            ]
        )
    ],
    no_gutters=True,
    justify='between',
    className="panel mb-2"
)

dim_block = dbc.Card(
    [
        dbc.Button(
            "Dimensionality Reduction",
            id="dim-reduce-collapse-btn",
            color='primary'
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    [
                        dbc.FormGroup(
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
                                        dbc.InputGroupAddon(
                                            [
                                                get_cog_btn(
                                                    i['value'] + '-btn')
                                                for i in dim_list
                                            ],
                                            addon_type="append"
                                        ),
                                        *dim_settings
                                    ]
                                )
                            ]
                        ),
                        dbc.FormGroup(
                            [
                                dbc.Label("Obtain 2D Embeddings via"),
                                dbc.InputGroup(
                                    [
                                        dbc.Select(
                                            options=vis_list,
                                            required=True,
                                            id="vis-methods-select",
                                            value=vis_list[0]['value']
                                        ),
                                        dbc.InputGroupAddon(
                                            get_cog_btn('vis-settings-btn'),
                                            addon_type="append"
                                        )
                                    ]
                                )
                            ]
                        )
                    ]
                ),
                dbc.CardFooter(
                    dbc.Button(
                        "Run",
                        block=True,
                        id="dim-run-btn",
                        color='secondary',
                        outline=True
                    ),
                    className="p-0"
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
                dbc.FormGroup(
                    [
                        dbc.Label("Obtain Labels via"),
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
                                dbc.InputGroupAddon(
                                    [
                                        get_cog_btn(i['value'] + '-btn')
                                        for i in clu_list
                                    ],
                                    addon_type="append"
                                ),
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
                dbc.FormGroup(
                    [
                        dbc.Label("Obtain Labels via"),
                        dbc.InputGroup(
                            [
                                dbc.Select(
                                    options=ssclu_list,
                                    required=True,
                                    id="ssclu-methods-select",
                                    value=ssclu_list[0]['value']
                                ),
                                dbc.InputGroupAddon(
                                    get_cog_btn('ssclu-settings-btn'),
                                    addon_type="append"
                                ),
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
                dbc.FormGroup(
                    [
                        dbc.Label("Obtain Labels via"),
                        dbc.InputGroup(
                            [
                                dbc.Select(
                                    options=lbt_list,
                                    required=True,
                                    id="lbt-methods-select",
                                    value=lbt_list[0]['value']
                                ),
                                dbc.InputGroupAddon(
                                    get_cog_btn('lbt-settings-btn'),
                                    addon_type="append"
                                )
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
        dbc.Tab(clu_block, label="Unsupervised"),
        dbc.Tab(ssclu_block, label="Semi-Supervised"),
        dbc.Tab(lbt_block, label="Label Transfer"),
    ]
)

label_block = dbc.Card(
    [
        dbc.Button(
            "Clustering",
            id="clustering-collapse-btn",
            color='primary'
        ),
        dbc.Collapse(
            [
                dbc.CardBody(label_tabs, className="p-0 mt-2"),
                dbc.CardFooter(
                    dbc.Button(
                        "Run",
                        block=True,
                        id="label-run-btn",
                        color='secondary',
                        outline=True
                    ),
                    className="p-0"
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
        dbc.Button(
            "Annotations",
            id="annotation-collapse-btn",
            color='primary'
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    dbc.Col(
                        [
                            dbc.InputGroup(
                                [
                                    dbc.InputGroupAddon(
                                        dbc.Select(
                                            options=[],
                                            id="main-annotation-select"
                                        ),
                                        id="main-annotation-addon",
                                        addon_type='prepend'
                                    ),
                                    dbc.InputGroupAddon(
                                        dbc.Select(
                                            options=[],
                                            id="side-annotation-select"
                                        ),
                                        id="side-annotation-addon",
                                        className="no-display",
                                        addon_type='prepend'
                                    ),
                                    dbc.Input(
                                        value="",
                                        placeholder="Annotation",
                                        type="text",
                                        id="annotation-input"
                                    ),
                                    dbc.InputGroupAddon(
                                        [
                                            dbc.Button(
                                                "Store",
                                                id="annotation-store-btn"
                                            )
                                        ],
                                        addon_type='append'
                                    )
                                ]
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

session_block = dbc.Card(
    [
        dbc.Button(
            "Session",
            id="session-collapse-btn",
            color='primary'
        ),
        dbc.Collapse(
            [
                dbc.CardBody(
                    [
                        dbc.Button(
                            "Import Session",
                            id='import-session-btn',
                            block=True,
                            color='secondary',
                            outline=True
                        ),
                        dbc.Button(
                            "Export Session",
                            id='export-session-btn',
                            block=True,
                            color='secondary',
                            outline=True
                        ),
                        dbc.Button(
                            "Export Annotations Only",
                            id='export-annotations-btn',
                            block=True,
                            color='secondary',
                            outline=True
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
    session_block
]
