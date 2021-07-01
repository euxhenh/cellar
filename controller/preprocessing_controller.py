import gc
import os
from ast import literal_eval

import dash
import numpy as np
import scanpy as sc
from BinToGene import BinToGene
from app import app, dbroot, logger
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate
from anndata import AnnData

from controller.cellar.core import read_adata
from controller.cellar.utils.exceptions import InvalidArgument


def _infer_bin_prefix(bin_name):
    prefix = ''
    for i in range(len(bin_name)):
        if bin_name[i].isalpha():
            prefix += bin_name[i]
        else:
            break
    return prefix


def correct_bin_names(bin_names):
    for i in range(len(bin_names)):
        bin_names[i] = bin_names[i].replace(
            ':', '_', bin_names[i].count(':')-1)
    return bin_names


def _validate_interval_extension(extend):
    if extend is None:
        return None

    extend = str(extend)

    if extend in ["0", "null", "None", "none", "-1", "0x"]:
        return None

    if extend[-1] == 'x':
        multiplier = extend[:-1]
        try:
            multiplier = literal_eval(multiplier)
            multiplier = float(multiplier)
        except Exception as e:
            logger.warn(str(e))
            raise InvalidArgument("Invalid extension parameter.")
        return str(multiplier) + 'x'
    else:
        try:
            multiplier = literal_eval(extend)
            multiplier = float(multiplier)
        except Exception as e:
            logger.warn(str(e))
            raise InvalidArgument("Invalid extension parameter.")
        return multiplier


@app.callback(
    Output("shape-signal-upload-prep", "data"),
    Output("feature-list-signal-prep", "data"),
    Output("temp-prep-h5", "className"),
    Output("data-loaded-plot-signal-prep", "data"),

    Input("prep-run-btn", "n_clicks"),

    # Checkboxes
    State("prep-filter-cells-counts-checkbox", "checked"),
    State("prep-filter-cells-genes-checkbox", "checked"),
    State("prep-filter-genes-counts-checkbox", "checked"),
    State("prep-filter-genes-cells-checkbox", "checked"),
    State("prep-high-var-checkbox", "checked"),
    State("prep-normalize-total-checkbox", "checked"),
    State("prep-log1p-checkbox", "checked"),
    State("prep-scale-checkbox", "checked"),

    # Parameters
    State("prep-filter-cells-counts-slider", "value"),
    State("prep-filter-cells-genes-slider", "value"),
    State("prep-filter-genes-counts-slider", "value"),
    State("prep-filter-genes-cells-slider", "value"),

    State("prep-high-var-flavor", "value"),
    State("prep-high-var-top-genes", "value"),
    State("prep-high-var-mean-min", "value"),
    State("prep-high-var-mean-max", "value"),
    State("prep-high-var-disp-min", "value"),
    State("prep-high-var-disp-max", "value"),

    State("prep-norm-target", "value"),
    State("prep-norm-maxf", "value"),
    State("prep-scale-zero", "value"),
    State("prep-scale-max", "value"),

    State("active-plot", "data"),
    prevent_initial_call=True
)
def run_prep(
        s1, cfcc, cfcg, cfgc, cfgc2, chvar, cnt, clog, cscale,
        sfcc, sfcg, sfgc, sfgc2,
        hvf, hvtg, hvmmin, hvmmax, hvdmin, hvdmax,
        nt, nmax, sz, smax,
        actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = "a1" if actp == 1 else "a2"

    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    try:
        if dbroot.adatas[an]['adata'].raw is not None:
            adata = dbroot.adatas[an]['adata'].to_memory(
            ).raw.to_adata().copy()
            logger.info("Found raw key for preprocessing.")
        else:
            adata = dbroot.adatas[an]['adata'].to_memory().copy()
            logger.info("No raw key found for preprocessing.")

        # copy old adata into raw
        adata.raw = adata

        if cfcc is not None and cfcc:
            sc.pp.filter_cells(adata, min_counts=sfcc[0])
            sc.pp.filter_cells(adata, max_counts=sfcc[1])
        if cfcg is not None and cfcg:
            sc.pp.filter_cells(adata, min_genes=sfcg[0])
            sc.pp.filter_cells(adata, max_genes=sfcg[1])
        if cfgc is not None and cfgc:
            sc.pp.filter_genes(adata, min_counts=sfgc[0])
            sc.pp.filter_genes(adata, max_counts=sfgc[1])
        if cfgc2 is not None and cfgc2:
            sc.pp.filter_genes(adata, min_cells=sfgc2[0])
            sc.pp.filter_genes(adata, max_cells=sfgc2[1])

        if chvar is not None and chvar:
            if hvdmax == 'inf':
                hvdmax = np.Inf
            if hvtg == "":
                hvtg = None
            sc.pp.highly_variable_genes(
                adata, flavor=hvf, n_top_genes=hvtg,
                min_mean=hvmmin, max_mean=hvmmax,
                min_disp=hvdmin, max_disp=hvdmax
            )

            adata = adata[:, adata.var.highly_variable]

        if cnt is not None and cnt:
            if nt == "":
                nt = None
            sc.pp.normalize_total(adata, target_sum=nt, max_fraction=nmax)

        if cscale is not None and cscale:
            if smax == "":
                smax = None
            sc.pp.scale(adata, zero_center=sz, max_value=smax)
    except Exception as e:
        logger.warn("An error occurred in preprocessing.")
        logger.warn(str(e))
        raise PreventUpdate

    filename = f'tmp/{an}/adata.h5ad'
    if os.path.exists(filename):
        os.remove(filename)
    adata.write(filename)

    del adata
    del dbroot.adatas[an]['adata']
    gc.collect()

    dbroot.adatas[an]['adata'] = read_adata(filename)

    return 1, 1, "", 1


@app.callback(
    Output("shape-signal-upload-prep-atac", "data"),
    Output("feature-list-signal-prep-atac", "data"),
    Output("temp-prep-h5-atac", "className"),
    Output("data-loaded-plot-signal-prep-atac", "data"),

    Input("prep-atac-run-btn", "n_clicks"),

    State("atac-operation", "value"),
    State("atac-extend-input", "value"),
    State("atac-max-extend-input", "value"),
    State("atac-strand-bool", "value"),
    State("atac-op-extend-input", "value"),
    State("atac-max-op-extend-input", "value"),

    State("active-plot", "data"),
    prevent_initial_call=True
)
def run_atac_prep(n1, aop, aei, amei, asb, aoei, amoei, actp):
    ctx = dash.callback_context
    if not ctx.triggered:
        raise PreventUpdate

    an = "a1" if actp == 1 else "a2"

    if an not in dbroot.adatas:
        raise PreventUpdate
    if 'adata' not in dbroot.adatas[an]:
        raise PreventUpdate

    try:
        if dbroot.adatas[an]['adata'].raw is not None:
            adata = dbroot.adatas[an]['adata'].to_memory(
            ).raw.to_adata().copy()
            logger.info("Found raw key for preprocessing.")
        else:
            adata = dbroot.adatas[an]['adata'].to_memory().copy()
            logger.info("No raw key found for preprocessing.")

        # copy old adata into raw
        adata_raw = adata

        btg = BinToGene(
            operation=aop,
            extend=_validate_interval_extension(aei),
            max_extend=_validate_interval_extension(amei),
            stream_direction=asb,
            op_extend=_validate_interval_extension(aoei),
            max_op_extend=_validate_interval_extension(amoei),
            n_jobs=4)

        prefix = _infer_bin_prefix(adata_raw.var_names.to_numpy()[0])

        counts, ids = btg.convert(
            adata_raw.X, adata_raw.var_names, prefix=prefix)

        adata = AnnData(counts)
        adata.var_names = ids
        adata.obs_names = adata_raw.obs_names.copy()
        adata.raw = adata_raw

        logger.info(f"Bin to gene matrix has shape {adata.shape}.")
    except Exception as e:
        logger.warn("An error occurred in preprocessing.")
        logger.warn(str(e))
        raise PreventUpdate

    filename = f'tmp/{an}/adata.h5ad'
    if os.path.exists(filename):
        os.remove(filename)
    adata.write(filename)

    del adata
    del dbroot.adatas[an]['adata']
    gc.collect()

    dbroot.adatas[an]['adata'] = read_adata(filename)

    return 1, 1, "", 1
