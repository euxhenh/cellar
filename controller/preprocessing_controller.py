import dash
import numpy as np
import plotly.express as px
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


@app.callback(
    Output("cell-gene-heatmap", "figure"),
    Input("preprocessing-mode-btn", "n_clicks")
)
def show_cell_gene_heatmap(n1):
    a = np.random.random((400, 500))

    fig = px.imshow(a)
    fig.update_layout(coloraxis_showscale=False)
    fig.update_xaxes(showticklabels=False)
    fig.update_yaxes(showticklabels=False)

    return fig
