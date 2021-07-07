import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from controller.cellar.utils.misc import get_server_dataset_dict
from gvars import (CELLAR_LOGO, DATA_PATH, DEMOS_URL, DOCS_URL, GITHUB_URL,
                   PAPER_URL)

from .preprocessing import prep_tabs

dataset_dict = get_server_dataset_dict(DATA_PATH)

cellar_bar = dbc.Nav(
    [
        dbc.NavItem(
            dbc.Button(
                dbc.Row(
                    [
                        html.Span(style={'width': '15px'}),
                        html.I(className="fas fa-filter"),
                        html.Span("Preprocessing", className="pl-3 pr-3"),
                    ],
                    align='baseline',
                    justify='center',
                    no_gutters=True
                ),
                id="preprocessing-mode-btn",
                block=True,
                color='secondary',
                outline=False,
                style={'font-size': '11px'}
            ),
            className="mr-2"
        ),
        dbc.NavItem(
            dbc.Button(
                dbc.Row(
                    [
                        html.Span(style={'width': '15px'}),
                        html.I(className="fas fa-clone"),
                        html.Span("Dual Mode", className="pl-3 pr-3"),
                    ],
                    align='baseline',
                    justify='center',
                    no_gutters=True
                ),
                id="dual-mode-btn",
                block=True,
                color='secondary',
                outline=False,
                style={'font-size': '11px'}
            ),
            className="mr-2"
        ),
        dbc.NavItem(
            dbc.Button(
                dbc.Row(
                    [
                        html.Span(style={'width': '15px'}),
                        html.I(className="fas fa-plus-square"),
                        html.Span("Load Data", className="pl-3 pr-3"),
                    ],
                    align='baseline',
                    justify='center',
                    no_gutters=True
                ),
                id="change-dataset-btn",
                block=True,
                color='primary',
                outline=False,
                style={'font-size': '13px'}
            ),
            className="mr-2"
        )
    ],
    className="ml-auto flex-nowrap mt-3 mt-md-0 m-auto",
    navbar=True
)


doc_bar = dbc.Nav(
    [
        dbc.NavItem(
            dbc.NavLink(
                dbc.Row(
                    [
                        html.I(className="fas fa-book"),
                        html.Span(style={'width': '5px'}),
                        html.Span("Documentation"),
                    ],
                    align='baseline',
                    no_gutters=True
                ),
                external_link=True,
                href=DOCS_URL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink(
                dbc.Row(
                    [
                        html.I(className="fab fa-github"),
                        html.Span(style={'width': '5px'}),
                        html.Span("GitHub"),
                    ],
                    align='baseline',
                    no_gutters=True
                ),
                external_link=True,
                href=GITHUB_URL,
                target="_blank")
        ),
        dbc.NavItem(
            dbc.NavLink(
                dbc.Row(
                    [
                        html.I(className="fab fa-youtube"),
                        html.Span(style={'width': '5px'}),
                        html.Span("Demo"),
                    ],
                    align='baseline',
                    no_gutters=True
                ),
                external_link=True,
                href=DEMOS_URL,
                target="_blank")
        )
    ],
    className="ml-auto flex-nowrap mt-3 mt-md-0 m-auto",
    navbar=True
)


navbar = dbc.Navbar(
    [
        dbc.Row(
            [
                html.A(
                    dbc.Row(
                        [
                            dbc.Col(html.Img(src=CELLAR_LOGO, height='30px')),
                            dbc.Col(dbc.NavbarBrand(
                                "Cellar", className="ml-2"))
                        ],
                        align='center',
                        no_gutters=True
                    ),
                    href=PAPER_URL,
                    target="_blank"
                ),
                html.Div(className="vline"),
                html.Span("Cell Type Annotation Tool", className="light-text")
            ],
            align='center',
            no_gutters=True,
            className="m-auto"
        ),
        cellar_bar
    ],
    color='primary',
    dark=True,
    className='mb-2 shadow-sm'
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
                        dbc.Button("Upload")
                    ],
                    className="mr-2",
                    max_size=1020*1024*1024  # 1 GB
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
            dbc.InputGroupAddon(
                [
                    dbc.Spinner(
                        children=[
                            html.Div(id="dataset-load"),
                        ],
                        type='grow',
                        fullscreen=True,
                        size='lg',
                        debounce=100
                    ),
                    dbc.Button(
                        "Load",
                        color="primary",
                        outline=False,
                        block=True,
                        id='load-dataset-btn'
                    )
                ],
                className='mw-150',
                addon_type="append"
            ),
            width='auto',
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
