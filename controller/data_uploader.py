import os
import io
import tarfile
import shutil
import gc
import glob

import anndata
import dash
import scanpy as sc
from base64 import b64decode
from dash.dependencies import Input, Output, State

from app import app, logger
from gvars import DATA_PATH
from .cellar.utils.misc import get_server_dataset_dict
from .cellar.utils.exceptions import UserError
from .multiplexer import MultiplexerOutput
from .notifications import _prep_notification


@app.callback(
    Output("dataset-dropdown", "options"),
    MultiplexerOutput("push-notification", "data"),

    Input('upload-dataset', 'contents'),
    State('upload-dataset', 'filename'),
    State('active-plot', "data"),
    prevent_initial_call=True
)
def upload_data(contents, filename, actp):
    # Don't allow / in filename
    filename = filename.split('/')[-1]

    an = "a1" if actp == 1 else 'a2'

    # if filename in os.listdir(os.path.join(DATA_PATH, 'uploaded')):
    #     error_msg = f"Dataset {filename} already exists. Skipping..."
    #     logger.warn(error_msg)
    #     return dash.no_update, _prep_notification(error_msg, "warning")

    content_type, content_string = contents.split(',')
    decoded = b64decode(content_string)

    if filename.endswith('.h5ad'):
        with open(os.path.join(DATA_PATH, 'uploaded', filename), "wb") as f:
            logger.info(f"Writing {filename} into uploaded directory.")
            try:
                f.write(io.BytesIO(decoded).read())
            except Exception as e:
                logger.error(str(e))
                error_msg = "Couldn't write h5ad file."
                return dash.no_update, _prep_notification(error_msg, "danger")
    elif filename.endswith('.tar.gz'):
        filename = filename[:-len('.tar.gz')] + ".h5ad"
        extract_path = f'tmp/{an}/d10x'

        logger.info("Extracting tar.gz file.")

        try:
            tar = tarfile.open(fileobj=io.BytesIO(decoded))
            if os.path.isdir(extract_path):
                shutil.rmtree(extract_path)
            tar.extractall(extract_path)
            tar.close()
        except Exception as e:
            logger.error(str(e))
            error_msg = "Couldn't extract tar.gz file."
            return dash.no_update, _prep_notification(error_msg, "danger")

        try:
            barcodes = glob.glob(f'tmp/{an}/d10x/**/barcodes.tsv', recursive=True)
            barcodes_gz = glob.glob(f'tmp/{an}/d10x/**/barcodes.tsv.gz', recursive=True)

            if len(barcodes) > 1 or len(barcodes_gz) > 1:
                raise UserError("More than one barcode file found.")
            if len(barcodes) == 0 and len(barcodes_gz) == 0:
                raise UserError("No barcodes.tsv file found.")
            barcodes = barcodes_gz[0] if len(barcodes) == 0 else barcodes[0]
            folder_name = os.path.dirname(barcodes)
            adata = sc.read_10x_mtx(folder_name, var_names='gene_symbols')
            adata.write(os.path.join(DATA_PATH, 'uploaded', filename))
        except UserError as e:
            logger.error(str(e))
            error_msg = str(e)
            return dash.no_update, _prep_notification(error_msg, "danger")
        except Exception as e:
            logger.error(str(e))
            error_msg = "Couldn't read Feature-Barcode-Matrix."
            return dash.no_update, _prep_notification(error_msg, "danger")

        shutil.rmtree(extract_path)

        del adata
        gc.collect()
    elif filename.endswith('.csv') or filename.endswith('.h5'):
        fpath = os.path.join(DATA_PATH, 'uploaded', filename)
        with open(fpath, "wb") as f:
            logger.info(f"Writing {filename} into uploaded directory.")
            try:
                f.write(io.BytesIO(decoded).read())
            except Exception as e:
                logger.error(str(e))
                error_msg = "Couldn't write file."
                return dash.no_update, _prep_notification(error_msg, "danger")

        try:
            if filename.endswith('.csv'):
                adata = sc.read_csv(fpath)
            else:
                adata = sc.read_10x_h5(fpath)
            adata.write_h5ad(
                os.path.join(DATA_PATH, 'uploaded', filename[:-4] + '.h5ad'))
        except Exception as e:
            logger.error(str(e))
            error_msg = "Couldn't read file."
            return dash.no_update, _prep_notification(error_msg, "danger")

        del adata
        gc.collect()
    else:
        error_msg = f"File format .{filename.split('.')[-1]} not implemented."
        logger.warn(error_msg)
        return dash.no_update, _prep_notification(error_msg, "warning")

    dataset_dict = get_server_dataset_dict(DATA_PATH)

    return [{'label': dataset_dict[d], 'value': d} for d in dataset_dict],\
        _prep_notification("Dataset uploaded.", "info")
