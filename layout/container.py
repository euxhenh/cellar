import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
import plotly.graph_objects as go
from controller.methods import (clu_list, dim_list, lbt_list, ssclu_list,
                                vis_list)

# from .greeting import greeting_page
from .plots import plots
from .sidebar import sidebar
from .analysis import main_analysis_tabs, side_analysis_tabs


main_body = dbc.Row(
    [
        dbc.Col(
            sidebar,
            className="p-0 sidebar pr-2",
            width=3, xs=12, sm=12, md=3, lg=3
        ),
        dbc.Col(
            [
                plots,
                dbc.Row(
                    [
                        dbc.Col(
                            main_analysis_tabs,
                            # width=6,
                            className="block-display",
                            id="main-analysis-col"
                        ),
                        dbc.Col(
                            side_analysis_tabs,
                            # width=6,
                            className="no-display",
                            id="side-analysis-col"
                        )
                    ]
                )
            ],
            className="main-body",
            width=9, xs=12, sm=12, md=9, lg=9
        ),
    ],
    no_gutters=True
)

main_page = dbc.Col(
    [
        dcc.Store(id="active-plot", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="dual-mode", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="shape-signal-upload", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="shape-signal-upload-prep", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="shape-signal-upload-prep-atac", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="shape-signal-atoggle", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="main-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="main-cluster-merge-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-cluster-merge-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="main-subset-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-subset-list-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="merge-plot-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="feature-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="feature-list-signal-prep", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="feature-list-signal-prep-atac", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="annotation-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="main-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="data-loaded-plot-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="data-loaded-plot-signal-prep", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="data-loaded-plot-signal-prep-atac", storage_type="memory",
                  data=1, clear_data=True),
        main_body
    ],
    id='main-body',
    # className="no-display"
)


# container = dbc.Row([greeting_page, main_page], id='container')
container = dbc.Row([main_page], id='container')
