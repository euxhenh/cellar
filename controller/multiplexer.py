from typing import List, Optional

from dash import callback_context
from dash.dependencies import ALL, Input, Output
from dash.exceptions import PreventUpdate
from dash import dcc

from app import app, dbroot


def get_triggered() -> Optional[str]:
    ctx = callback_context
    if not ctx.triggered:
        return
    return ctx.triggered[0]['prop_id'].split('.')[0]


def create_multiplexer(
        output_component: str, output_field: str, count: int) -> List:
    store_component = f'{output_component}-{output_field}-store'

    @app.callback(
        # Output is an arbitrary dash component
        Output(output_component, output_field),
        # Input is always a dcc.Store()
        Input({'id': store_component, 'idx': ALL}, 'data'),
        prevent_initial_call=True
    )
    def multiplexer(_):
        triggered = get_triggered()
        if triggered is None:
            raise PreventUpdate

        inputs = callback_context.inputs
        for k, v in inputs.items():
            id_ = k.split('.')[0]
            if id_ == triggered:
                return v

        raise PreventUpdate

    return [
        dcc.Store({'id': store_component, 'idx': idx}, data=None)
        for idx in range(count)]


def MultiplexerOutput(output_component: str, output_field: str):
    store_component = f'{output_component}-{output_field}-store'
    dbroot.MULTIPLEXER_OUTPUTS[
        store_component] = idx = dbroot.MULTIPLEXER_OUTPUTS.setdefault(
        store_component, -1) + 1
    return Output({'id': store_component, 'idx': idx}, 'data')
