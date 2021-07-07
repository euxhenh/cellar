import os
import io
import tarfile
import shutil
import gc

import anndata
import scanpy as sc
from base64 import b64decode
from dash.dependencies import Input, Output, State
from dash.exceptions import PreventUpdate

from app import app, logger
from gvars import DATA_PATH
from controller.cellar.utils.misc import get_server_dataset_dict


@app.callback(
    Output("dataset-dropdown", "options"),
    Input('upload-dataset', 'contents'),
    State('upload-dataset', 'filename'),
    State('active-plot', "data"),
    prevent_initial_call=True
)
def upload_data(contents, filename, actp):
    # Don't allow / in filename
    filename = filename.split('/')[-1]

    an = "a1" if actp == 1 else 'a2'

    if filename in os.listdir(os.path.join(DATA_PATH, 'uploaded')):
        logger.warn(f"Dataset {filename} found in path. Skipping...")
        raise PreventUpdate

    content_type, content_string = contents.split(',')
    decoded = b64decode(content_string)

    if filename.endswith('.h5ad'):
        with open(os.path.join(DATA_PATH, 'uploaded', filename), "wb") as f:
            logger.info(f"Writing {filename} into uploaded directory.")
            f.write(io.BytesIO(decoded).read())
    elif filename.endswith('.tar.gz'):
        filename = filename[:-len('.tar.gz')] + ".h5ad"
        extract_path = f'tmp/{an}/d10x'

        logger.info("Extracting tar.gz file.")
        tar = tarfile.open(fileobj=io.BytesIO(decoded))
        if os.path.isdir(extract_path):
            shutil.rmtree(extract_path)
        tar.extractall(extract_path)
        tar.close()

        adata = sc.read_10x_mtx(extract_path, var_names='gene_symbols')
        adata.write(os.path.join(DATA_PATH, 'uploaded', filename))

        shutil.rmtree(extract_path)

        del adata
        gc.collect()
    elif filename.endswith('.csv'):
        fpath = os.path.join(DATA_PATH, 'uploaded', filename)
        with open(fpath, "wb") as f:
            logger.info(f"Writing {filename} into uploaded directory.")
            f.write(io.BytesIO(decoded).read())

        adata = anndata.read_csv(fpath)
        adata.write_h5ad(
            os.path.join(DATA_PATH, 'uploaded', filename[:-4] + '.h5ad'))
        del adata
        gc.collect()
    else:
        logger.warn(f"File format {filename} not implemented.")
        raise PreventUpdate

    dataset_dict = get_server_dataset_dict(DATA_PATH)

    return [
        {'label': dataset_dict[d], 'value': d}
        for d in dataset_dict
    ]
