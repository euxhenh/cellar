import os
import tifffile
import pandas as pd
import numpy as np
import matplotlib

from .exceptions import InvalidArgument
from .colors import palette_to_rgb

from app import logger


def get_name_index(x, y, im_names):
    """
    Given integer coordinates x and y, and a list of image filenames,
    return the index in the list that corresponds to tile (x, y).
    """
    name = f'R001_X00{x+1}_Y00{y+1}.tif'
    if name not in im_names:
        raise InvalidArgument("Invalid image filenames encountered.")
    return im_names.index(name)


def generate_tile(
        path_to_tiff, path_to_df, adata=None, palette=None, savepath=None):
    """
    Given a path to image files, construct a codex tile and return figure.
    Parameters
    __________
    path_to_tiff: string, path to a folder with the Tiff Image files
    path_to_df: string, path to data.csv
    adata: AnnData object containing labels
    palette: list of strings containing colors
    savepath: path to save the generated tile
    Returns
    _______
    The generated tile as a numpy array (n, m, channels)
    """
    im_names = os.listdir(path_to_tiff)
    im_names = [i for i in im_names if i[-3:] == 'tif']
    ims = [tifffile.imread(os.path.join(path_to_tiff, i)) for i in im_names]
    data = pd.read_csv(path_to_df)

    # Whether to use colors or not
    has_labels = False
    if adata is not None:
        if 'labels' in adata.obs:
            has_labels = True
            if np.min(adata.obs['labels']) < 0:
                raise InvalidArgument("Negative labels found.")

    # Begin Grid Construction
    grid = []
    grid_x_len = np.max(data['tile_x'])  # assuming min = 0
    grid_y_len = np.max(data['tile_y'])  # assuming min = 0
    palr, palb, palg = palette_to_rgb(palette)

    for y in range(grid_y_len+1):  # columns first
        grid_row = []
        for x in range(grid_x_len+1):
            i = get_name_index(x, y, im_names)

            my_cells = data[(data['tile_x'] == x) & (data['tile_y'] == y)]
            # same z for entire tile
            top_z = my_cells.loc[my_cells.index[0]]['z']
            tile = (ims[i][top_z, 0]).astype(int)  # filled cells

            if has_labels:
                cell_ids = my_cells.loc[my_cells.index]['id']
                cell_rids = my_cells.loc[my_cells.index]['rid']
                labels = adata.obs['labels'].to_numpy()[cell_rids] + 1
                labels = np.hstack((0, np.array(labels)))
                tile = labels[tile]
            else:
                tile[tile != 0] = -1  # if no labels, fill cells with white

            tile[ims[i][top_z, 2] != 0] = 0  # blackout the cell boundaries

            # Color tile and append
            grid_row.append(np.array([palr[tile], palb[tile], palg[tile]]))

        grid.append(np.concatenate(grid_row, axis=-1))
        logger.info(f'Finished tile row {y}')

    grid = np.moveaxis(np.hstack(grid), 0, -1)

    if savepath is not None:
        matplotlib.image.imsave(savepath, grid.astype(np.uint8))

    return grid.astype(np.uint8)
