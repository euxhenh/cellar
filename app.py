import BTrees.OOBTree
import dash
import dash_bootstrap_components as dbc
import ZODB
import ZODB.FileStorage

from log import setup_logger

db = ZODB.DB(None)
connection = db.open()
dbroot = connection.root
dbroot.adatas = BTrees.OOBTree.BTree()


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


app = dash.Dash(external_stylesheets=external_stylesheets, title="Cellar")
logger = setup_logger('Cellar')
