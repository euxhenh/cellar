import dash
import numpy as np
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate


def get_de_cluster_update_func(prefix, an):
    def _func(s1):
        ctx = dash.callback_context
        if not ctx.triggered:
            raise PreventUpdate

        if an in dbroot.adatas:
            if 'adata' in dbroot.adatas[an]:
                if 'labels' not in dbroot.adatas[an]['adata'].obs:
                    return [], None, [], None

                unq_labels = np.unique(
                    dbroot.adatas[an]['adata'].obs['labels'].to_numpy())
                value = prefix + "-cluster" + str(unq_labels[0])

                return [{
                    "label": "Cluster " + str(i),
                    "value": prefix + "-cluster" + str(i)
                } for i in unq_labels], value, [{
                    "label": "Cluster " + str(i),
                    "value": prefix + "-cluster" + str(i)
                } for i in unq_labels], value

        return [], None, [], None
    return _func


# Loop over both plots
for prefix, an in zip(['main', 'side'], ['a1', 'a2']):
    app.callback(
        Output(prefix + '-de-cluster-select', "options"),
        Output(prefix + '-de-cluster-select', "value"),
        Output(prefix + '-annotation-select', "options"),
        Output(prefix + '-annotation-select', "value"),
        Input(prefix + '-cluster-list-signal', "data"),
        prevent_initial_call=True
    )(get_de_cluster_update_func(prefix, an))
