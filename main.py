import dash_bootstrap_components as dbc

from app import app
from layout.container import container
from layout.navbar import navbar, modes_bar
from controller.data_loader import *
from controller.ui_controller import *
from controller.operations import *
from controller.update_plot import *
from controller.cluster_controller import *
from controller.annotations import *

app.layout = dbc.Container(
    [
        navbar,
        modes_bar,
        container
    ],
    fluid=True
)


if __name__ == "__main__":
    app.run_server(debug=True)
