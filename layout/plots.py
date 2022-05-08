import dash_bootstrap_components as dbc
from dash import html, dcc

from .misc import empty_figure
from app import dbroot


def get_plot_palette_popover(prefix):
    palette = dbroot.palettes[prefix]
    return dbc.Popover([
        dbc.PopoverHeader("Color Palette"),
        dbc.PopoverBody(
            [
                dbc.Row([
                    dbc.Col(
                        dbc.InputGroup(
                            [
                                dbc.InputGroupText(
                                    ("0" if i * 5 + j < 10 else "") +
                                    str(i * 5 + j),
                                    className="input-group-text",
                                    style={
                                        'background-color': palette[
                                            i * 5 + j]
                                    },
                                    id=prefix + f"-color-bg-{i* 5 + j}"
                                ),
                                dbc.Input(
                                    id=prefix + f"-color-input-{i* 5 + j}",
                                    debounce=True,
                                    maxLength=6,
                                    minLength=6,
                                    pattern='^(?:[0-9a-fA-F]{6})$',
                                    placeholder=str(palette[i * 5 + j])[1:]
                                )
                            ]
                        )
                    ) for j in range(5)
                ], className="") for i in range(10)
            ] + [
                dbc.Row(
                    [
                        dbc.Col([
                            dbc.Button(
                                "Clear All",
                                id=prefix + "-clear-palette-btn",
                                color='secondary',
                                className="me-3"
                            ),
                            dbc.Button(
                                "Apply Palette",
                                id=prefix + "-apply-palette-btn",
                                color='primary'
                            )
                        ], width=4)
                    ],
                    justify='center',
                    className="mt-3 g-0"
                )
            ]
        )
    ],
        className="color-popover",
        id=prefix + "-color-palette-popover",
        is_open=False,
        placement="bottom-start",
        target=prefix + "-change-color-palette-btn",
        trigger="legacy"
    )


def get_plot_download_popover(prefix):
    return dbc.Popover(
        [
            dbc.PopoverHeader("Download Plot"),
            dbc.PopoverBody(
                [
                    dbc.InputGroup(
                        [
                            dbc.InputGroupText("Format"),
                            dbc.Select(
                                options=[
                                    {"label": "svg", "value": "svg"},
                                    {"label": "png", "value": "png"},
                                    {"label": "jpg", "value": "jpg"},
                                    {"label": "webp", "value": "webp"},
                                    {"label": "pdf", "value": "pdf"}
                                ],
                                value="svg",
                                id=prefix + "-plot-download-format"
                            )
                        ],
                        className="mb-2"
                    ),
                    dbc.Row(
                        [
                            dbc.Col(
                                dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Width"),
                                        dbc.Input(
                                            value=1000,
                                            type="number",
                                            id=prefix + "-plot-download-width"
                                        )
                                    ]
                                ),
                                className="mr-2"
                            ),
                            dbc.Col(
                                dbc.InputGroup(
                                    [
                                        dbc.InputGroupText("Height"),
                                        dbc.Input(
                                            value=1000,
                                            type="number",
                                            id=prefix + "-plot-download-height"
                                        )
                                    ]
                                )
                            )
                        ],
                        className="mb-2 g-0"
                    ),
                    dbc.InputGroup(
                        [
                            dbc.InputGroupText("Scale"),
                            dbc.Input(
                                value=1,
                                type="number",
                                id=prefix + "-plot-download-scale"
                            )
                        ],
                        className="mb-2"
                    ),
                    dcc.Loading(
                        dcc.Download(id=prefix + "-download-plot-buf")
                    ),
                    html.Div(
                        dbc.Button(
                            "Download Plot",
                            color='primary',
                            id=prefix + "-plot-download-btn"
                        ),
                        className="d-grid gap-2"
                    )
                ]
            )
        ],
        id=prefix + "-plot-download-popover",
        is_open=False,
        placement="bottom-end",
        target=prefix + "-download-plot-popover-btn",
        trigger="legacy"
    )


def get_plot(prefix):
    title = "PLOT" if prefix == 'main' else "PLOT 2"

    plot = dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row(
                    [
                        dbc.Col(
                            title, className="ml-3", align='center',
                            width='auto'
                        ),
                        dbc.Col(
                            [
                                dbc.Button(
                                    html.I(
                                        className="fas fa-palette",
                                        style={'font-size': '16px'}
                                    ),
                                    id=prefix + "-change-color-palette-btn",
                                    className="mr-3",
                                    color="light"
                                ),
                                dbc.Tooltip(
                                    "Change Color Palette",
                                    target=prefix + "-change-color-palette-btn"
                                ),
                                get_plot_palette_popover(prefix),
                                dbc.Button(
                                    html.I(
                                        className="fas fa-save",
                                        style={'font-size': '16px'}
                                    ),
                                    id=prefix + "-download-plot-popover-btn",
                                    className="mr-3",
                                    color='light'
                                ),
                                dbc.Tooltip(
                                    "Save Plot to a File",
                                    target=prefix + "-download-plot-popover-btn"
                                ),
                                get_plot_download_popover(prefix),
                                dbc.Button(
                                    "Activate",
                                    id=prefix + "-activate-btn",
                                    disabled=True,
                                    color='light'
                                ),
                                dbc.Tooltip(
                                    "Make this the Active Plot",
                                    target=prefix + "-activate-btn"
                                ),
                            ],
                            width='auto'
                        )
                    ],
                    align='center',
                    justify='between',
                    className="g-0",
                ),
                id=prefix + '-cardheader',
                className="active-plot-bg" if prefix == 'main'
                else "non-active-plot-bg",
                style={'padding-top': '5px', 'padding-bottom': '5px'}
            ),
            dbc.CardBody(
                [
                    dcc.Loading(
                        children=[dcc.Graph(
                            id=prefix + '-plot',
                            figure=empty_figure,
                            loading_state={'component_name': prefix + '-plot'},
                            config={
                                'autosizable': True,
                                'displaylogo': False,
                                'edits': {
                                    'titleText': True,
                                    # 'legendText': True,
                                    'colorbarTitleText': True
                                },
                                'toImageButtonOptions': {
                                    'format': 'svg'
                                }
                            }
                        )],
                        type="circle"
                    )
                ]
            )
        ],
        className="shadow-sm"
    )

    return plot


plots = dbc.Row(
    [
        dbc.Col(
            get_plot("main"),
            md=12, lg=12,
            id="main-plot-col"
        ),
        dbc.Col(
            get_plot("side"),
            md=0, lg=0,
            id="side-plot-col",
            className="no-display"
        )
    ],
    className="mb-5 g-0",
)
