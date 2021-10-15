import os
import tifffile
import pandas as pd
import numpy as np
import matplotlib
import json
import matplotlib.pyplot as plt

from .exceptions import InvalidArgument
from .colors import palette_to_rgb
from skimage import draw
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

    # In case no z-planes are provided, we add an extra dimension
    for i in range(len(ims)):
        if len(ims[i].shape) == 3:
            ims[i] = np.expand_dims(ims[i], 0)
    # Tiff files are now assumed to have shape (z-planes, 4, h, w)
    # The 4 channels stand for
    # 0) cell seg.filled, 1) nucleus seg. filled
    # 2) cell seg. outline 3) nucleus seg. outline
    # We only use channels 0 and 2 to generate a tile
    CELL_FILL_CH = 0
    CELL_OUTLINE_CH = 2

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
    owner = []
    grid_x_len = np.max(data['tile_x'])  # assuming min = 0
    grid_y_len = np.max(data['tile_y'])  # assuming min = 0

    for y in range(grid_y_len+1):  # columns first
        grid_row = []
        owner_row = []
        for x in range(grid_x_len+1):
            i = get_name_index(x, y, im_names)

            my_cells = data[(data['tile_x'] == x) & (data['tile_y'] == y)]
            # same z for entire tile
            if len(my_cells) == 0:
                grid_row.append(np.zeros((3, *ims[i][0, 0].shape)))
                owner_row.append(np.zeros(ims[i][0, 0].shape) - 1)
                continue

            top_z = my_cells.loc[my_cells.index[0]]['z']  # get z of first cell
            tile = (ims[i][top_z, CELL_FILL_CH]).astype(int)  # filled cells
            tile_cp = tile.copy()

            if has_labels:
                # cell_ids = my_cells.loc[my_cells.index]['id']
                cell_rids = my_cells.loc[my_cells.index]['rid']
                cell_ids = my_cells.loc[my_cells.index]['id']

                sort_idx = np.argsort(cell_ids.to_numpy() + 1)
                idx = np.searchsorted(
                    cell_ids.to_numpy(), tile, sorter = sort_idx)
                tile = np.arange(len(cell_ids))[sort_idx][idx]
                rid_tile = cell_rids.to_numpy()[tile]
                rid_tile[tile_cp == 0] = -1 # remove empty pixels

                cell_rids = cell_rids[cell_rids < adata.shape[0]]
                labels = adata.obs['labels'].to_numpy()[cell_rids] + 1

                if labels.size > 0:
                    # since 0 means no cell
                    tile[tile != 0] = labels[tile[tile != 0]]
                else:
                    tile.fill(0)
            else:
                tile[tile != 0] = -1  # if no labels, fill cells with white

            # blackout the cell boundaries
            tile[ims[i][top_z, CELL_OUTLINE_CH] != 0] = 0

            # Color tile and append
            palr, palb, palg = palette_to_rgb(palette, tile.max())
            grid_row.append(np.array([palr[tile], palb[tile], palg[tile]]))
            owner_row.append(rid_tile)

        grid_row = np.concatenate(grid_row, axis=-1)
        owner_row = np.hstack(owner_row)
        grid.append(grid_row)
        owner.append(owner_row)
        logger.info(f'Finished tile row {y}')

    grid = np.hstack(grid)
    grid = np.moveaxis(grid, 0, -1)
    owner = np.vstack(owner)

    if savepath is not None:
        matplotlib.image.imsave(savepath, grid.astype(np.uint8))

    return grid.astype(np.uint8), owner.astype(int)


def generate_10x_spatial(
        path_to_img, path_to_df, path_to_json,
        adata=None, savepath=None, in_tissue=True):
    '''
    in_tissue:
        True: onyl show spots that are in the tissue
        False: show all spots
    '''

    has_labels = False
    if adata is not None:
        if 'labels' in adata.obs:
            has_labels = True

    if has_labels:
        labels = adata.obs['labels']
    else:
        return

    if 'json_dict' in adata.uns:
        dic = adata.uns['json_dict']
    else:
        json_file = open(path_to_json, 'r')
        dic = json.load(json_file)
    scaling_factor = dic['tissue_hires_scalef']  # 0.08250825
    full_d = dic['spot_diameter_fullres']
    r = full_d*scaling_factor/2

    if 'image' in adata.uns:
        small_img = adata.uns['image']
    else:
        small_img = np.array(plt.imread(path_to_img))

    if 'spatial_dict' in adata.uns:
        spatial_dict = adata.uns['spatial_dict']
    else:
        spatial_dict = {}
        points = []
        row_col_dict = {}
        spatial_info = pd.read_csv(
            path_to_df, delimiter=",", header=None).values
        for row in spatial_info:
            if row[1] == 0 and in_tissue:
                continue
            spatial_dict[row[0]] = [row[-2], row[-1]]
            points.append([row[-2], row[-1]])
            row_col_dict[(row[2], row[3])] = row[0]

    owner = np.zeros(small_img.shape[:2]) - 1
    barcodes = list(adata.obs['barcodes'])
    R, G, B = palette_to_rgb()
    for i in range(len(adata.obs)):
        center = np.array(spatial_dict[(barcodes[i])])*scaling_factor
        center = center.astype('int')

        label = labels[i] + 1  # +1 so that dont get white
        color = np.array([R[label], G[label], B[label]])
        circle = draw.disk(center, r)

        small_img[circle[0], circle[1], :] = color
        owner[circle[0], circle[1]] = i

    if savepath is not None:
        plt.imsave(savepath, small_img)

    return small_img, owner.astype(int)
