import dash_bootstrap_components as dbc
import os

from app import app
from layout.container import container
from layout.navbar import navbar, modes_bar
from controller.data_loader import *
from controller.ui_controller import *
from controller.operations import *
from controller.update_plot import *
from controller.cluster_controller import *
from controller.annotations import *
from controller.analysis_controller import *
from controller.session import *
from controller.preprocessing_controller import *
from controller.data_uploader import *
from controller.subsets import *

app.layout = dbc.Container(
    [
        navbar,
        modes_bar,
        container
    ],
    fluid=True
)

if 'SHINYPROXY_PUBLIC_PATH' in os.environ:
    app.config.update({
        'routes_pathname_prefix': os.environ['SHINYPROXY_PUBLIC_PATH'],
        'requests_pathname_prefix': os.environ['SHINYPROXY_PUBLIC_PATH']
    })

if __name__ == "__main__":
    dev = False
    app.run_server(debug=dev, use_reloader=dev, host='0.0.0.0')
