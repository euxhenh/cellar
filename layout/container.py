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
            className="p-0 sidebar"
        ),
        dbc.Col(
            [
                plots,
                dbc.Row(
                    [
                        dbc.Col(
                            main_analysis_tabs,
                            width=12,
                            className="block-display",
                            id="main-analysis-col"
                        ),
                        dbc.Col(
                            side_analysis_tabs,
                            width=12,
                            className="no-display",
                            id="side-analysis-col"
                        )
                    ]
                )
            ],
            width=9
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
        dcc.Store(id="shape-signal-atoggle", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="main-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-cluster-list-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="feature-list-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="annotation-signal", storage_type="memory",
                  data=1, clear_data=True),

        dcc.Store(id="main-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),
        dcc.Store(id="side-plot-signal-code", storage_type="memory",
                  data=1, clear_data=True),

        # *[dcc.Store(id='signal-' + m['value'], storage_type='memory',
        #             clear_data=True) for m in dim_list],
        # *[dcc.Store(id='signal-' + m['value'], storage_type='memory',
        #             clear_data=True) for m in clu_list],
        # nav_bar,
        main_body
    ],
    id='main-body',
    # className="no-display"
)


# container = dbc.Row([greeting_page, main_page], id='container')
container = dbc.Row([main_page], id='container')
