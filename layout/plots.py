import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .misc import empty_figure


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
                        no_gutters=True,
                        className="mb-2"
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
                    dbc.Button(
                        "Download Plot",
                        block=True,
                        color='primary',
                        id=prefix + "-plot-download-btn"
                    )
                ]
            )
        ],
        id=prefix + "-plot-download-popover",
        is_open=False,
        placement="bottom-end",
        target=prefix + "-download-plot-popover-btn",
        trigger="click"
    )


def get_plot(prefix):
    title = "PLOT" if prefix == 'main' else "PLOT 2"

    plot = dbc.Card(
        [
            dbc.CardHeader(
                dbc.Row(
                    [
                        dbc.Row(
                            title, className="ml-3", align='center'
                        ),
                        dbc.Row(
                            [
                                dbc.Button(
                                    html.I(
                                        className="fas fa-save",
                                        style={'font-size': '16px'}
                                    ),
                                    id=prefix + "-download-plot-popover-btn",
                                    className="mr-3",
                                    color='light'
                                ),
                                get_plot_download_popover(prefix),
                                dbc.Button(
                                    "Activate",
                                    id=prefix + "-activate-btn",
                                    disabled=True,
                                    color='light'
                                )
                            ]
                        )
                    ],
                    align='center',
                    justify='between',
                    no_gutters=True
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
            width=12,
            id="main-plot-col"
        ),
        dbc.Col(
            get_plot("side"),
            width=0,
            id="side-plot-col",
            className="no-display"
        )
    ],
    className="mb-5",
    no_gutters=True
)
