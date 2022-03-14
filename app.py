import os

import BTrees.OOBTree
import dash
import dash_bootstrap_components as dbc
import ZODB
import ZODB.FileStorage

from log import setup_logger
from controller.cellar.utils.colors import PALETTE


# In-memory database. Used for storing AnnData objects
db = ZODB.DB(None)
connection = db.open()
dbroot = connection.root
dbroot.adatas = BTrees.OOBTree.BTree()
dbroot.MULTIPLEXER_OUTPUTS = BTrees.OOBTree.BTree()
dbroot.notifications = BTrees.OOBTree.BTree()
dbroot.palettes = BTrees.OOBTree.BTree()
dbroot.palettes['main'] = PALETTE.copy()
dbroot.palettes['side'] = PALETTE.copy()

# Bootstrap theme
external_stylesheets = [
    dbc.themes.SANDSTONE,
    {
        'href': 'https://use.fontawesome.com/releases/v5.8.1/css/all.css',
        'rel': 'stylesheet',
        'integrity': 'sha384-50oBUHEmvpQ+1lW4y57PTFmhCaXp0ML' +
        '5d60M1M7uH2+nqUivzIebhndOJK28anvf',
        'crossorigin': 'anonymous'
    }
]

"""
Dash requires the knowledge of the path used to access the app.
ShinyProxy makes this path available as an environment variable,
which we expose to dash below.
"""
if 'SHINYPROXY_PUBLIC_PATH' in os.environ:
    routes_pathname_prefix = os.environ['SHINYPROXY_PUBLIC_PATH']
    requests_pathname_prefix = os.environ['SHINYPROXY_PUBLIC_PATH']
else:
    routes_pathname_prefix = None
    requests_pathname_prefix = None

app = dash.Dash(
    external_stylesheets=external_stylesheets,
    title="Cellar",
    suppress_callback_exceptions=True,
    routes_pathname_prefix=routes_pathname_prefix,
    requests_pathname_prefix=requests_pathname_prefix,
    )
logger = setup_logger('Cellar')
