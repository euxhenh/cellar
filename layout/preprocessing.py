import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

from .misc import empty_prep_figure

row_1 = dbc.Row(
    [
        dbc.Col(
            dbc.Card(
                [
                    dbc.CardBody(
                        [
                            html.H5("Filter Cells", className="card-title"),
                            dbc.Row(
                                [
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Min Counts",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Max Counts",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                ],
                                no_gutters=True
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Min Genes",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                value=200,
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Max Genes",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                value=3000,
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1")
                                ],
                                no_gutters=True
                            )
                        ]
                    )
                ]
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
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Min Counts",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Max Counts",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                ],
                                no_gutters=True
                            ),
                            dbc.Row(
                                [
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Min Cells",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                value=3,
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1"),
                                    dbc.Col(dbc.InputGroup(
                                        [
                                            dbc.InputGroupAddon(
                                                "Max Cells",
                                                addon_type="prepend",
                                                className="prep-addon"),
                                            dbc.Input(
                                                type="number",
                                                placeholder="None"),
                                        ]
                                    ), width=6, className="p-1")
                                ],
                                no_gutters=True
                            )
                        ]
                    )
                ]
            ),
            width=3
        ),
        dbc.Col(
            dbc.Card(
                [
                    dbc.CardBody(
                        [
                            html.H5("Stats", className="card-title"),
                            dcc.Graph(
                                figure=empty_prep_figure,
                                id="cell-gene-heatmap"
                            )
                        ]
                    )
                ]
            ),
            width=6
        ),
    ],
    no_gutters=True
)


prep = html.Div([
    row_1
], className="mb-3")
