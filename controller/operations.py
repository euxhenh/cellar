import dash
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from layout.method_settings.clu_settings import clu_settings_keys
from layout.method_settings.dim_settings import dim_settings_keys
from layout.method_settings.ssclu_settings import ssclu_settings_keys
from layout.method_settings.vis_settings import vis_settings_keys
from layout.method_settings.lbt_settings import lbt_settings_keys
from layout.method_settings.intg_settings import intg_settings_keys

from .cellar.utils.exceptions import InternalError
from .methods import (clu_list, dim_list, find_method, lbt_list, ssclu_list,
                      vis_list, intg_list)


def get_button_switch_func(m_list):
    def _func(val):
        is_opens = [False] * len(m_list)
        styles = [{'display': 'none'} for _ in range(len(m_list))]

        if val is None:
            styles[0]['display'] = 'block'
            return *styles, *is_opens

        for i, method in enumerate(m_list):
            if method['value'] == val:
                break
        else:
            raise PreventUpdate

        styles[i]['display'] = 'block'
        return *styles, *is_opens
    return _func


for m_list, m_name in zip(
    [dim_list, clu_list, vis_list, ssclu_list, lbt_list, intg_list],
    ['dim', 'clu', 'vis', 'ssclu', 'lbt', 'intg']
):
    app.callback(
        [Output(m['value'] + '-btn', 'style') for m in m_list],
        [Output(m['value'] + '-settings', 'is_open') for m in m_list],
        Input(m_name + "-methods-select", "value")
    )(get_button_switch_func(m_list))


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
            raise InternalError("'value' key not found in child.")

        kwargs[method_settings_keys[key]] = child['value']

    return kwargs


def get_filter(settings_keys, m_list, key, x_to_use):
    """
    Returns a function that searches for 'method' and its parameters
    in the popover components. That function itself returns the function
    that is used to run the chosen method with the chosen parameters.
    """
    def _func(adata, method, settings, extras=None):
        if method not in settings_keys:
            raise InternalError(f"{method} not found in settings keys.")

        kwargs = _search_settings(settings_keys[method], settings)
        kwargs['key'] = key
        kwargs['x_to_use'] = x_to_use
        if extras is not None:
            kwargs['extras'] = extras
        return find_method(m_list, method)['func'](adata, **kwargs)
    return _func


# functions
dim_reduce_filter = get_filter(
    dim_settings_keys, dim_list, key='x_emb', x_to_use='x')
vis_filter = get_filter(
    vis_settings_keys, vis_list, key='x_emb_2d', x_to_use='x_emb')
clu_filter = get_filter(
    clu_settings_keys, clu_list, key='labels', x_to_use='x_emb')
ssclu_filter = get_filter(
    ssclu_settings_keys, ssclu_list, key='labels', x_to_use='x_emb')
lbt_filter = get_filter(
    lbt_settings_keys, lbt_list, key='labels', x_to_use='x_emb')
intg_filter = get_filter(
    intg_settings_keys, intg_list, key=None, x_to_use=None)
