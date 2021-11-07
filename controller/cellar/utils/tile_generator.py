import os
from numpy.core.fromnumeric import sort
import tifffile
import pandas as pd
import numpy as np
import matplotlib
import json
import matplotlib.pyplot as plt
import plotly.colors as pc

from .exceptions import InternalError, InvalidArgument, UserError
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


def _read_tiff_images(path_to_tiff):
    """
    Reads images into an list and makes sure each has a z-plane.
    """
    im_names = os.listdir(path_to_tiff)
    # only take .tif files
    im_names = [i for i in im_names if i[-3:] == 'tif']
    ims = [tifffile.imread(os.path.join(path_to_tiff, i)) for i in im_names]

    # In case no z-planes are provided, we add an extra dimension
    for i in range(len(ims)):
        if ims[i].ndim == 3:
            ims[i] = np.expand_dims(ims[i], 0)
        if ims[i].shape != ims[0].shape:
            raise UserError("Found tiles of differing shapes.")
    return ims, im_names


def _local_owner_2_global(local_owner, cell_ids, cell_rids):
    """
    Maps ids in local_owner to rids.
    E.g. if ids = [1, 2, 3], rids = [7, 8, 9]
    and local_owner is [2, 2, 3, 1, 1],
    it will return [8, 8, 9, 7, 7].
    """
    sort_idx = cell_ids.argsort()
    local_owner_idx = sort_idx[
        np.searchsorted(cell_ids, local_owner, sorter=sort_idx)]
    global_owner = cell_rids[local_owner_idx]
    return global_owner


def _get_global_owner_single_image(
        im, data, no_owner=-1, outline=-2, CELL_FILL_CH=0, CELL_OUTLINE_CH=2):
    """
    Given the coordinates of a single image, process it
    and return an ownership tile. All pixels not belonging to a cell
    will be filled with `no_owner`.

    im: array
        A single image of shape (n z_planes, 4, height, width)
         The 4 channels of a tile stand for
            0) cell seg.filled, 1) nucleus seg. filled
            2) cell seg. outline 3) nucleus seg. outline
    data: pd.DataFrame
    no_owner: int
        Default value to fill any pixels which do not belong to a cell.
    outline: int
        Default value to fill any pixels which belong to a cell outline.
    """
    n_zplanes, seg_dim, height, width = im.shape
    if seg_dim != 4:
        raise UserError(
            "Incorrect number of channels in tile. Need to have " +
            "4 channels consisting of cell seg. filled, nucleus seg. " +
            "filled, cell seg. outline, and nucleus seg. outline.")
    # If no cells in this image
    if len(data) == 0:
        logger.info("Empty tile encountered.")
        return np.full((height, width), no_owner)

    # Find best z-plane. All cells in a single image are assumed
    # to have the same best z-plane, so we only look at the first cell.
    top_z = data.loc[data.index[0]]['z']
    if top_z >= n_zplanes:
        raise UserError(
            "Best z-plane is greater than total number of z planes.")

    # Has shape (height, width) and each pixels contains
    # the cell ID that the pixel belongs to.
    # Need to convert cell IDs to cell RIDs to determine global owner.
    local_owner = (im[top_z, CELL_FILL_CH]).astype(int)
    cell_rids = data['rid'].to_numpy().astype(int)
    cell_ids = data['id'].to_numpy().astype(int)
    if cell_ids.min() <= 0:
        raise UserError(
            "Found cell with local ID 0 which is reserved for no-cell.")
    if cell_rids.min() < 0:  # No equality for RIDS
        raise UserError(
            "Found cell with negative RID.")
    # First map each pixel ID to its index in cell_ids
    global_owner = _local_owner_2_global(local_owner, cell_ids, cell_rids)
    global_owner[local_owner == 0] = no_owner  # 0 IDs mean no owner
    # Fill cell boundaries
    outlines = im[top_z, CELL_OUTLINE_CH]
    global_owner[outlines > 0] = outline  # assume outline if positive value

    return global_owner


def generate_tile(
        path_to_tiff, path_to_df, adata=None,
        colors=None, palette=None, savepath=None):
    """
    Given a path to image files, construct a codex tile and return figure.
    Parameters
    __________
    path_to_tiff: string, path to a folder with the Tiff Image files
    path_to_df: string, path to data.csv
    adata: anndata.AnnData object or None
        If present, will store owner under adata.uns['spatial_idx']
    palette: list of strings containing colors
    savepath: path to save the generated tile
    Returns
    _______
    The generated tile as a numpy array (n, m, channels)
    and an ownership array (n, m)
    """
    # Tiff files are now assumed to have shape (z-planes, 4, h, w)
    # The 4 channels stand for
    # 0) cell seg.filled, 1) nucleus seg. filled
    # 2) cell seg. outline 3) nucleus seg. outline
    # We only use channels 0 and 2 to generate a tile
    NO_OWNER = -1
    OUTLINE = -2
    WHITE = 255
    BLACK = 0

    if adata is not None and 'spatial_idx' in adata.uns:
        owner = adata.uns['spatial_idx'].astype(int)
    else:
        ims, im_names = _read_tiff_images(path_to_tiff)  # list of images
        data = pd.read_csv(path_to_df)
        # Begin Grid Construction
        grid = []
        owner = []  # Needed to return which cell a pixel belongs to
        grid_x_len = np.max(data['tile_x']) + 1  # assuming min = 0
        grid_y_len = np.max(data['tile_y']) + 1  # assuming min = 0

        for y in range(grid_y_len):  # columns first
            owner_row = []
            for x in range(grid_x_len):
                # Make sure we are looking at the right image
                i = get_name_index(x, y, im_names)
                my_cells = data[(data['tile_x'] == x) & (data['tile_y'] == y)]
                # Convert the image into a tile of shape (height, width)
                # where each pixel is assigned a RID value.
                # -1 means no owner, and -2 means cell boundary.
                global_owner = _get_global_owner_single_image(
                    ims[i], my_cells, no_owner=NO_OWNER, outline=OUTLINE)
                owner_row.append(global_owner)
            logger.info(f"Finished row {y}.")
            owner.append(owner_row)
        owner = np.block(owner)
        if adata is not None:
            adata.uns['spatial_idx'] = owner

    if colors is None:
        tile = owner.copy()  # Will hold colors
        # Fill all cells with white and set boundaries to black
        tile[owner >= 0] = WHITE
        tile[owner == NO_OWNER] = BLACK
        tile[owner == OUTLINE] = BLACK
        tile = tile.astype(np.uint8)
        # Make it a 3 channel image, channels last
        tile = np.stack([tile, tile, tile], axis=-1)
        if savepath is not None:
            matplotlib.image.imsave(savepath, tile)
        return tile, owner

    if not isinstance(colors, np.ndarray):
        raise InternalError("colors is not a numpy array.")
    if owner.max() >= colors.shape[0]:
        owner[owner >= colors.shape[0]] = NO_OWNER
        logger.warn("Found RIDs greater than the total number of samples. " +
                    "These RIDs will be filtered.")

    tile = owner.copy()
    if colors.dtype == int:
        # Assume these are cluster IDs
        if colors.min() < 0:
            raise InternalError("Negative cluster IDs found.")
        # so assign to each cell a color from the given palette
        # NOTE: palette gets extended at the top with BLACK and
        # end with WHITE
        palr, palb, palg = palette_to_rgb(palette, colors.max())
        # Assign each pixel to its cluster ID
        # Add one not to mess up cells with ID 0 and no owner
        tile[tile >= 0] = colors[tile[tile >= 0]] + 1
        tile[tile < 0] = 0  # Black out all special keys, i.e., outlines
        # 0 will pick up black in the palette, and 1 will pick up white
        tile = np.stack([palr[tile], palb[tile], palg[tile]], axis=-1)
        tile = tile.astype(np.uint8)
        if savepath is not None:
            matplotlib.image.imsave(savepath, tile)
        return tile, owner

    if colors.dtype == float:
        tile = tile.astype(float)
        tile[owner >= 0] = colors[owner[owner >= 0]]
        # Make sure special keys take the smallest value
        min_color = colors.min()
        tile[owner == NO_OWNER] = min_color
        tile[owner == OUTLINE] = min_color

        if savepath is not None:
            matplotlib.image.imsave(savepath, tile, cmap='magma')
        return tile, owner

    raise InternalError("Should not reach this point.")


def generate_10x_spatial(
        path_to_img, path_to_df, path_to_json,
        adata=None, savepath=None, in_tissue=True, palette=None):
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
    R, G, B = palette_to_rgb(palette)
    for i in range(len(adata.obs)):
        center = np.array(spatial_dict[(barcodes[i])])*scaling_factor
        center = center.astype('int')

        label = labels[i] + 1  # +1 so that dont get white
        color = np.array([R[label], G[label], B[label]])
        circle = draw.disk(center, r)

        small_img[circle[0], circle[1], :] = color
        owner[circle[0], circle[1]] = i

    if savepath is not None:
        matplotlib.image.imsave(savepath, small_img)

    return small_img, owner.astype(int)
