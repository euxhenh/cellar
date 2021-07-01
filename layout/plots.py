import dash_bootstrap_components as dbc
import dash_core_components as dcc
# import dash_html_components as html

from .misc import empty_figure


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
                            dbc.Button(
                                "Activate",
                                id=prefix + "-activate-btn",
                                disabled=True,
                                color='light'
                            )
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
                            config={'autosizable': True, 'displaylogo': False}
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
