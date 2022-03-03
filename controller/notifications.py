import dash
from dash import html
import dash_bootstrap_components as dbc
import time
from app import app, dbroot
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from .multiplexer import MultiplexerOutput


def _prep_notification(text, icon="info"):
    # icon can be in ["primary", "secondary", "success",
    #          "warning", "danger", "info", "light", "dark"]
    notification = {
        'text': text,
        'icon': icon
    }
    return notification


def _prep_notification_card(key, value, add_mb=True):
    return dbc.Toast(
        html.P(value['text']),
        header=time.strftime("%H:%M:%S", time.localtime(key)),
        icon=value["icon"],
        className="shadow-none"
    )


@app.callback(
    MultiplexerOutput("notifications-toast", "children"),
    MultiplexerOutput("notifications-toast", "is_open"),
    MultiplexerOutput("notifications-toast", "duration"),
    Input("toggle-notifications-btn", "n_clicks"),
    State("notifications-toast", "is_open"),
    prevent_initial_call=True
)
def open_toast(n1, is_open):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    nots = []

    for key in list(reversed(sorted(list(dbroot.notifications.keys())))):
        nots.append(_prep_notification_card(key, dbroot.notifications[key]))

    nots = nots[:5]

    return nots, not is_open, None


@app.callback(
    MultiplexerOutput("notifications-toast", "children"),
    MultiplexerOutput("notifications-toast", "is_open"),
    MultiplexerOutput("notifications-toast", "duration"),

    Input("push-notification", "data"),
    prevent_initial_call=True
)
def push_notification(pn):
    ctx = dash.callback_context
    if not ctx.triggered or pn is None:
        raise PreventUpdate

    tt = time.time()
    text = _prep_notification_card(tt, pn, False)
    dbroot.notifications[tt] = pn

    return text, True, 20000
