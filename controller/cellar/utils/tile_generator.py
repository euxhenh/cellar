import os
import tifffile
import pandas as pd
import numpy as np
import matplotlib
import json
import matplotlib.pyplot as plt
import matplotlib.patches as patches

from .exceptions import InvalidArgument
from .colors import palette_to_rgb, PALETTE

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
        small_img = plt.imread(path_to_img)

    if 'spatial_dict' in adata.uns:
        spatial_dict = adata.uns['spatial_dict']
    else:
        spatial_dict = {}
        points = []
        row_col_dict = {}
        spatial_info = pd.read_csv(
            path_to_df, delimiter=",", header=None).values
        for row in spatial_info:
            if row[1] == 0 and in_tissue == True:
                continue
            spatial_dict[row[0]] = [row[-2], row[-1]]
            points.append([row[-2], row[-1]])
            row_col_dict[(row[2], row[3])] = row[0]

    # Display the image
    ax = plt.subplot()
    ax.imshow(small_img)
    label_dict = {}
    barcodes = list(adata.obs['barcodes'])
    for i in range(len(adata.obs)):
        # Create a Rectangle patch
        label_dict[barcodes[i]] = labels[i]
        center = np.array(spatial_dict[(barcodes[i])]) * scaling_factor
        center = center.astype('int')
        color = PALETTE[labels[i]]
        circle = patches.Circle(
            (center[1], center[0]),
            r,
            linewidth=1,
            edgecolor=color,
            facecolor=color,
            alpha=0.6)

        # Add the patch to the Axes
        ax.add_patch(circle)

    fig = ax.get_figure()
    if savepath is not None:
        fig.savefig(savepath)

    # fig2data
    fig.canvas.draw()
    tile = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
    tile = tile.reshape(fig.canvas.get_width_height()[::-1] + (3,))
    # end of fig2data

    return tile
