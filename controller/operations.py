import dash
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from layout.method_settings.clu_settings import clu_settings_keys
from layout.method_settings.dim_settings import dim_settings_keys
from layout.method_settings.vis_settings import vis_settings_keys

from .cellar.utils.exceptions import InternalError
from .methods import clu_list, dim_list, vis_list, find_method


@app.callback(
    [Output(m['value'] + '-btn', 'style') for m in dim_list],
    [Output(m['value'] + '-settings', 'is_open') for m in dim_list],
    Input("dim-methods-select", "value")
)
def switch_dim_settings_button(val):
    is_opens = [False] * len(dim_list)
    styles = [{'display': 'none'} for _ in range(len(dim_list))]

    if val is None:
        styles[0]['display'] = 'block'
        return *styles, *is_opens

    for i, method in enumerate(dim_list):
        if method['value'] == val:
            break
    else:
        raise PreventUpdate

    styles[i]['display'] = 'block'
    return *styles, *is_opens


@app.callback(
    [Output(m['value'] + '-btn', 'style') for m in vis_list],
    [Output(m['value'] + '-settings', 'is_open') for m in vis_list],
    Input("vis-methods-select", "value")
)
def switch_vis_settings_button(val):
    is_opens = [False] * len(vis_list)
    styles = [{'display': 'none'} for _ in range(len(vis_list))]

    if val is None:
        styles[0]['display'] = 'block'
        return *styles, *is_opens

    for i, method in enumerate(vis_list):
        if method['value'] == val:
            break
    else:
        raise PreventUpdate

    styles[i]['display'] = 'block'
    return *styles, *is_opens




@app.callback(
    [Output(m['value'] + '-btn', 'style') for m in clu_list],
    [Output(m['value'] + '-settings', 'is_open') for m in clu_list],
    Input("clu-methods-select", "value")
)
def switch_clu_settings_button(val):
    is_opens = [False] * len(clu_list)
    styles = [{'display': 'none'} for _ in range(len(clu_list))]

    if val is None:
        styles[0]['display'] = 'block'
        return *styles, *is_opens

    for i, method in enumerate(clu_list):
        if method['value'] == val:
            break
    else:
        raise PreventUpdate

    styles[i]['display'] = 'block'
    return *styles, *is_opens


def _recur_search(var, key):
    """
    Given a mixture of lists, tuples, dicts as var, try to find
    a dictionary which has a key 'id' = key in it and return it.
    """
    if isinstance(var, tuple) or isinstance(var, list):
        for v in var:
            for res in _recur_search(v, key):
                yield res
    elif isinstance(var, dict):
        if 'id' in var:
            if var['id'] == key:
                yield var
        if 'children' in var:
            for res in _recur_search(var['children'], key):
                yield res
        if 'props' in var:
            for res in _recur_search(var['props'], key):
                yield res


def _search_settings(method_settings_keys, settings):
    """
    We maintain a dictionary of dimensionality reduction methods
    in dim_settings_keys where each key (method) stores another
    dictionary (md) holding that method's settings (parameters).
    The keys of md are component ids and the values are parameter
    names that will be passed to dim_reduce.

    For example, dim_settings_keys['dim-PCA'] holds a dictionary
        dim_pca_settings_keys = {
            'dim-PCA-n-components': 'n_components',
            'dim-PCA-whiten': 'whiten',
            'dim-PCA-solver': 'svd_solver',
            'dim-PCA-random-state': 'random_state'
        }
    where the keys (dim-PCA-key) is the widget id and the value is the
    parameter name to pass to sklearn's PCA.

    Parameters
    __________
    method_settings_keys: dict
        Dicionary holding setting id's and parameter names.

    settings: tuple of list of dicts of ...
        Holds all children in the method-settings elements. This
        is a mixture of lists, tuples, and dicts. We recursively search
        this element to find the children with id's as determined by
        dim_pca_settings_keys[dim_method]. By doing this we avoid having
        to write new Input elements into our callbacks every time we add
        a new setting. All that needs to be done is add the setting's
        id into the settings_keys dict and it will be parsed automatically.

    """
    kwargs = {}


    for key in method_settings_keys:
        child = next(_recur_search(settings, key))

        # if there exists a component with 'key' as its 'id'
        # then child should never be None. 'value' may be missing
        # if not manually specified when constructing the widget.
        if child is None or 'value' not in child:
            raise InternalError(f"'value' key not found in child.")

        kwargs[method_settings_keys[key]] = child['value']

    return kwargs


def get_filter(settings_keys, m_list, key):
    """
    Returns a function that searches for 'method' and its parameters
    in the popover components. That function itself returns the function
    that is used to run the chosen method with the chosen parameters.
    """
    def _func(adata, method, settings):
        if method not in settings_keys:
            raise InternalError(f"{method} not found in settings keys.")

        kwargs = _search_settings(settings_keys[method], settings)
        kwargs['key'] = key
        if (key=='x_emb_2d'):
            if 'n_components' in kwargs.keys():
                kwargs['n_components']=2
            if 'n_evecs' in kwargs.keys():
                kwargs['n_evecs']=2
            if 'n_clusters' in kwargs.keys():
                kwargs['n_clusters']=2
        return find_method(m_list, method)['func'](adata, **kwargs)
    return _func


# functions
dim_reduce_filter = get_filter(dim_settings_keys, dim_list, 'x_emb')
clu_filter = get_filter(clu_settings_keys, clu_list, 'labels')

vis_reduce_filter = get_filter(vis_settings_keys, vis_list, 'x_emb_2d')

