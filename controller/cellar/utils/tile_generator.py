import os
import tifffile
import pandas as pd
import numpy as np
import matplotlib

from .exceptions import InvalidArgument

PALETTE = [
    "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F",
    "#B3B3B3", "#E5C494", "#9C9DBB", "#E6946B", "#DA8DC4", "#AFCC63",
    "#F2D834", "#E8C785", "#BAB5AE", "#D19C75", "#AC9AAD", "#CD90C5",
    "#B8C173", "#E5D839", "#ECCA77", "#C1B7AA", "#BBA37E", "#BC979D",
    "#C093C6", "#C1B683", "#D8D83E", "#F0CD68", "#C8BAA5", "#A6AB88",
    "#CC958F", "#B396C7", "#CBAB93", "#CCD844", "#F3D05A", "#CFBCA1",
    "#90B291", "#DC9280", "#A699C8", "#D4A0A3", "#BFD849", "#F7D34B",
    "#D6BF9C", "#7BBA9B", "#EC8F71", "#999CC9", "#DD95B3", "#B2D84E",
    "#FBD63D", "#DDC198"
]


def palette_to_rgb(palette=None):
    # Colors when labels are provided
    if palette is None:
        palette = PALETTE
    # add 0 for black and -1 for white
    palette = ["#000000"] + palette + ["#FFFFFF"]
    palr = np.array([int(i[1:3], 16) for i in palette])  # hex string to int
    palb = np.array([int(i[3:5], 16) for i in palette])
    palg = np.array([int(i[5:], 16) for i in palette])
    return palr, palb, palg


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
        print(f'Finished tile row {y}')

    grid = np.moveaxis(np.hstack(grid), 0, -1)

    if savepath is not None:
        matplotlib.image.imsave(savepath, grid.astype(np.uint8))

    return grid.astype(np.uint8)
