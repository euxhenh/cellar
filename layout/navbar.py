import dash_bootstrap_components as dbc
from dash import html, dcc
from controller.cellar.utils.misc import get_server_dataset_dict
from controller.multiplexer import create_multiplexer
from gvars import (CELLAR_LOGO, DATA_PATH, DEMOS_URL, DOCS_URL, GITHUB_URL,
                   PAPER_URL, EMAIL, BUGS)

from .preprocessing import prep_tabs

dataset_dict = get_server_dataset_dict(DATA_PATH)

cellar_bar = dbc.Nav([
    dbc.NavItem(
        dbc.Button(
            [
                html.I(className="fas fa-filter me-2"),
                "Preprocessing"
            ],
            id="preprocessing-mode-btn",
            color='secondary',
            outline=False,
            className="d-flex align-items-center text-nowrap",
        ),
        className="me-2"
    ),
    dbc.NavItem(
        dbc.Button(
            [
                html.I(className="fas fa-clone me-2"),
                "Dual Mode"
            ],
            id="dual-mode-btn",
            color='secondary',
            outline=False,
            className="d-flex align-items-center text-nowrap",
        ),
        className="me-2"
    ),
    dbc.NavItem(
        dbc.Button(
            [
                html.I(className="fas fa-plus-square me-2"),
                "Load Data"
            ],
            id="change-dataset-btn",
            color='primary',
            outline=False,
            className="d-flex align-items-center text-nowrap",
            style={'font-size': '13px'},
        ),
        className="me-4"
    ),
    dbc.NavItem([
        dbc.Button(
            html.I(className="fas fa-bell"),
            outline=True,
            color="secondary",
            id="toggle-notifications-btn"
        ),
        dbc.Tooltip(
            "Notifications",
            target="toggle-notifications-btn",
        ),
        dbc.Toast(
            id="notifications-toast",
            header="Notifications",
            is_open=False,
            dismissable=True,
            style={
                "position": "fixed",
                "top": 70, "right": 20
            }
        ),
        *create_multiplexer('push-notification', 'data', 120),
        *create_multiplexer('notifications-toast', 'children', 2),
        *create_multiplexer('notifications-toast', 'is_open', 2),
        *create_multiplexer('notifications-toast', 'duration', 2),
        *create_multiplexer('main-other-features-change', 'data', 20),
        *create_multiplexer('side-other-features-change', 'data', 20)
    ]),
    # dbc.NavItem(
    #     [
    #         dbc.Button(
    #             html.I(className="fas fa-power-off"),
    #             outline=True,
    #             color='danger',
    #             id="shutdown-btn",
    #             style={
    #                 'font-size': '20px',
    #                 'position': 'absolute',
    #                 'right': 20
    #             }
    #         ),
    #         dbc.Tooltip(
    #             "End Session",
    #             target="shutdown-btn",
    #         )
    #     ]
    # ),
],
    # className="ml-auto flex-nowrap mt-3 mt-md-0 m-auto",
    navbar=True,
    # fill=True,
    # justified=True
)


navbar = dbc.Navbar(
    dbc.Container([
        html.A(
            dbc.Row(
                [
                    dbc.Col(html.Img(src=CELLAR_LOGO, height="30px")),
                    dbc.Col(dbc.NavbarBrand("Cellar", className="ms-2")),
                    dbc.Col(html.Div(className="vline")),
                    dbc.Col(html.Span(
                        "Cell Type Annotation Tool",
                        className="light-text",
                        style={'font-size': '12px'}))
                ],
                align="center",
                className="g-0 text-nowrap",
            ),
            href=PAPER_URL,
            target="_blank",
            style={"textDecoration": "none"},
        ),
        cellar_bar
    ]),
    color='primary',
    dark=True,
    className='mb-2 shadow-sm'
)


doc_bar = dbc.Nav(
    [
        dbc.NavItem(
            dbc.NavLink([
                html.I(className="fas fa-book me-2"),
                html.Span("Documentation"),
            ],
                external_link=True,
                href=DOCS_URL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink([
                html.I(className="fab fa-github me-2"),
                html.Span("GitHub"),
            ],
                external_link=True,
                href=GITHUB_URL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink([
                html.I(className="fab fa-youtube me-2"),
                html.Span(style={'width': '5px'}),
                html.Span("Demos"),
            ],
                external_link=True,
                href=DEMOS_URL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink([
                html.I(className="fas fa-envelope me-2"),
                html.Span(style={'width': '5px'}),
                html.Span("Email Us"),
            ],
                className="ml-5",
                external_link=True,
                href=EMAIL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink([
                html.I(className="fab fa-github me-2"),
                html.Span(style={'width': '5px'}),
                html.Span("Issues"),
            ],
                external_link=True,
                href=BUGS,
                target="_blank")
        ),
    ],
    className="ml-auto flex-nowrap mt-3 mt-md-0 m-auto",
    navbar=True
)


footer = dbc.Navbar(
    doc_bar,
    color='light',
    dark=False,
    className='m-auto shadow-sm'
)


dataset_bar = dbc.Row(
    [
        dbc.Col(
            [
                dcc.Upload(
                    id="upload-dataset",
                    children=[
                        dbc.Button("Upload", color="secondary")
                    ],
                    className="mr-2",
                    max_size=1040*1024*1024  # 1 GB
                ),
            ],
            width='auto'
        ),
        dbc.Col(
            dcc.Loading(
                dbc.Select(
                    options=[
                        {'label': dataset_dict[d], 'value': d}
                        for d in dataset_dict
                    ],
                    className='mw-800',
                    placeholder="Select dataset",
                    id='dataset-dropdown',
                    style={
                        'height': '32.5px',
                        'font-size': '0.8125rem'
                    }
                ),
                fullscreen=True
            ),
            width='auto'
        ),
        dbc.Col(
            dbc.Row([
                dbc.Col(html.Div(
                    dbc.Button(
                        "Load",
                        color="primary",
                        outline=False,
                        id='load-dataset-btn'
                    )
                ), className="me-2"),
                dbc.Col(dcc.Loading(id="dataset-load"))
            ], justify='center', align='center'),
            width='auto'
        )
    ],
    justify='center',
    align='center'
)


modes_bar = html.Div(
    [
        dbc.Collapse(
            dbc.Card(dbc.CardBody(dataset_bar), className="mb-2"),
            id="collapsible-dataset-bar",
        ),
        dbc.Collapse(
            prep_tabs,
            id="collapsible-prep"
        )
    ]
)
