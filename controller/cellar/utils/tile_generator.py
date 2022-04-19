import os
import tifffile
import pandas as pd
import numpy as np
import matplotlib
import json
import matplotlib.pyplot as plt

from .misc import _filter_outliers
from .exceptions import InternalError, InvalidArgument, UserError
from .colors import palette_to_rgb, interpolate_grayimage
from skimage import draw
from app import logger


def get_owner_from_coordinates(x, y, NO_OWNER=-1, OUTLINE=-2, cell_radius=10):
    """
    Generates a simple tile given only x and y coordinates.
    """
    x, y = np.array(x), np.array(y)
    xmin, xmax = x.min(), x.max()
    ymin, ymax = y.min(), y.max()
    x = x - xmin
    y = y - ymin
    width = xmax - xmin
    height = ymax - ymin
    n_samples = len(x)

    n_samples_per_row = int((width / height * n_samples) ** (1/2))
    n_samples_per_col = int(n_samples / n_samples_per_row)
    im_height = n_samples_per_col * cell_radius * 2
    im_width = n_samples_per_row * cell_radius * 2
    x_cord = (x * im_width / width).astype(int)
    y_cord = (y * im_height / height).astype(int)

    im = np.full((im_height + cell_radius + 1,
                 im_width + cell_radius + 1), NO_OWNER)
    for i, (x, y) in enumerate(zip(y_cord, x_cord)):
        circle = draw.disk((x, y), max(1, int(cell_radius * 2 / 3)))
        im[circle[0], circle[1]] = i

    return im


def get_name_index(x, y, im_names):
    """
    Given integer coordinates x and y, and a list of image filenames,
    return the index in the list that corresponds to tile (x, y).
    """
    name = f'R001_X00{x+1}_Y00{y+1}.tif'
    if name not in im_names:
        raise InvalidArgument("Invalid image filenames encountered.")
    return im_names.index(name)


def _validate_codex_dataframe(data):
    """
    Checks if data.csv has the right columns.
    id maps to the id in a tile, and rid maps this id to the index in adata
    """
    cols = ['id', 'rid', 'z', 'tile_x', 'tile_y']
    if not set(cols).issubset(data.columns):
        raise UserError("Columns missing from dataframe. All of " +
                        "['id', 'rid', 'z', 'tile_x', 'tile_y'] " +
                        "must be provided.")
    if not (data.dtypes[cols] == int).all():
        raise UserError("Columns were not of integer type.")


def _read_tiff_images(path_to_tiff):
    """
    Reads images into an list and makes sure each has a z-plane.
    """
    if not os.path.isdir(path_to_tiff):
        raise UserError("No images folder found. Make sure all your tif " +
                        "files are stored in a subdirectory named 'images'.")
    im_names = os.listdir(path_to_tiff)
    # only take .tif files
    im_names = [i for i in im_names if i[-3:] == 'tif']
    if len(im_names) == 0:
        raise UserError("No tif images were found.")
    ims = [tifffile.imread(os.path.join(path_to_tiff, i)) for i in im_names]

    # In case no z-planes are provided, we add an extra dimension
    for i in range(len(ims)):
        if ims[i].ndim == 3:
            ims[i] = np.expand_dims(ims[i], 0)
        if ims[i].ndim != 4:
            raise UserError("Encountered incorrect image dimensions.")
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
        A single image of shape (n_z_planes, 4, height, width)
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

    # Has shape (height, width), and each pixel contains
    # the cell ID that the pixel belongs to.
    # Need to convert cell IDs to cell RIDs to determine global owner.
    local_owner = (im[top_z, CELL_FILL_CH]).astype(int)
    cell_rids = data['rid'].to_numpy().astype(int)
    cell_ids = data['id'].to_numpy().astype(int)
    if cell_ids.min() <= 0:
        raise UserError(
            "Found cell with local ID 0 which is reserved for no-cell.")
    if cell_rids.min() < 0:  # No equality for RIDS
        raise UserError("Found cell with negative RID.")
    # First map each pixel ID to its index in cell_ids, then to RID
    global_owner = _local_owner_2_global(local_owner, cell_ids, cell_rids)
    global_owner[local_owner == 0] = no_owner  # 0 IDs mean no owner
    # Fill cell boundaries
    outlines = im[top_z, CELL_OUTLINE_CH]
    global_owner[outlines > 0] = outline  # assume outline if positive value

    return global_owner


def get_owner_from_paths(path_to_tiff, path_to_df, NO_OWNER=-1, OUTLINE=-2):
    ims, im_names = _read_tiff_images(path_to_tiff)  # list of images
    if not os.path.exists(path_to_df):
        raise UserError(
            "No dataframe was found. Please upload " +
            "a 'data.csv' file along with an 'images' subfolder.")
    data = pd.read_csv(path_to_df)
    _validate_codex_dataframe(data)
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
    try:
        owner = np.block(owner)  # Stitch all tiles together
    except ValueError as ve:
        raise UserError("Tile row dimensions do not match. Make sure " +
                        "images form a valid rectangular tile.")

    return owner


def generate_tile(
        path_to_tiff, path_to_df, adata=None,
        colors=None, palette=None, savepath=None, key='spatial_idx'):
    """
    Given a path to image files, construct a codex tile and return figure.
    Parameters
    __________
    path_to_tiff: str
        Path to a folder with the Tiff Image files
    path_to_df: str
        Path to a data.csv containing columns [id, rid, z, tile_x, tile_y]
    adata: anndata.AnnData object
        If present, will store owner under adata.uns[key]
    palette: list of str
        Contains colors in hex format
    savepath: str
        Path to save the generated tile

    Returns
    _______
    The generated tile as a numpy array (n, m, channels),
    and an ownership array (n, m) where every pixel corresponds to the
    cell ID it belongs to.
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

    if adata is not None and key in adata.uns:
        owner = adata.uns[key].astype(int).copy()
        logger.info("Using cached tile.")
    else:
        if 'x' in adata.obs and 'y' in adata.obs:
            logger.info("Reading x and y coordinates from adata.")
            x = adata.obs['x'].to_numpy().astype(float)
            y = adata.obs['y'].to_numpy().astype(float)
            owner = get_owner_from_coordinates(x, y)
        else:
            logger.info("Reading tile information from tiff files.")
            owner = get_owner_from_paths(
                path_to_tiff, path_to_df, NO_OWNER, OUTLINE)
        if adata is not None:
            adata.uns[key] = owner.copy()  # Make sure to copy

    logger.info(f"Found {len(np.unique(owner) - 2)} cells in the tile.")

    if colors is None:
        # Fill all cells with white and set the rest to black
        tile = np.where(owner >= 0, WHITE, BLACK).astype(np.uint8)
        # Make it a 3 channel image, channels last
        tile = np.repeat(tile[..., np.newaxis], 3, axis=-1)
        if savepath is not None:
            matplotlib.image.imsave(savepath, tile)
        return tile, owner

    if not isinstance(colors, np.ndarray):
        raise InternalError("colors is not a numpy array.")
    if owner.max() >= colors.shape[0]:
        owner[owner >= colors.shape[0]] = NO_OWNER
        logger.warn("Found RIDs greater than the total number of samples. " +
                    "These RIDs will be ignored.")

    tile = owner.copy()
    if colors.dtype == int:
        if colors.min() < 0:
            raise InternalError("Negative cluster IDs found.")
        # Assume these are cluster IDs
        # so assign to each cell a color from the given palette
        # NOTE: palette gets extended at the top with BLACK and
        # end with WHITE
        palr, palb, palg = palette_to_rgb(palette, colors.max())
        # Assign each pixel to its cluster ID
        # Add one not to mess up cells with ID 0 and no owner
        tile[tile >= 0] = colors[tile[tile >= 0]] + 1
        tile = tile.clip(min=0)  # Black out all special keys, i.e., outlines
        # 0 will pick up black in the palette which has sufficient colors
        tile = np.stack([palr[tile], palb[tile], palg[tile]], axis=-1)
        tile = tile.astype(np.uint8)
        if savepath is not None:
            matplotlib.image.imsave(savepath, tile)
        return tile, owner

    if colors.dtype == float:
        # This can be protein expression
        tile = tile.astype(float)
        # Filter very high and very low values for better visualization
        colors = _filter_outliers(colors)
        tile[owner >= 0] = colors[owner[owner >= 0]]
        # Make sure special keys take the smallest value - some 10% to
        # distinguish them from cells
        tile[(owner == NO_OWNER) | (owner == OUTLINE)
             ] = colors.min() - np.ptp(colors) * 0.1
        tile = interpolate_grayimage(tile, 'magma')
        if savepath is not None:
            matplotlib.image.imsave(savepath, tile, cmap='magma')
        return tile, owner

    raise InternalError(
        "Error in colors format. Please report this to the developers.")


def _read_verify_10x_json(path_to_json):
    """
    This is a dictionary containing 'tissue_hires_scalef' and
    'spot_diameter_fullres' keys.
    """
    if path_to_json is None or not os.path.exists(path_to_json):
        raise UserError("No json file was provided.")
    try:
        with open(path_to_json, 'r') as jf:
            json_dict = json.load(jf)
    except:
        raise UserError("Could not read json file `scalefactors_json.json`.")

    if any([k not in json_dict
            for k in ['tissue_hires_scalef', 'spot_diameter_fullres']]):
        raise UserError("json_dict is missing one of 'tissue_hires_scalef' " +
                        "or 'spot_diameter_fullres'.")
    return json_dict


def _read_verify_10x_image(path_to_img):
    """
    A tile of shape (height, width, channels).
    """
    if path_to_img is None or not os.path.exists(path_to_img):
        raise UserError("No image file was provided.")
    try:
        image = np.array(plt.imread(path_to_img))
    except:
        raise UserError("Could not read image file `detected_tissue_image.jpg`.")
    return image


def _read_verify_10x_df(path_to_df, in_tissue=True):
    """
    Must be a csv file where the index corresponds to the sample ID
    and contains two columns for x and y coordinates. Assume
    the two coordinates are the last columns.
    """
    if path_to_df is None or not os.path.exists(path_to_df):
        raise UserError("No image file was provided.")
    try:
        df = pd.read_csv(path_to_df, delimiter=",", header=None).values
        spatial_dict = {}
        for row in df:
            if row[1] == 0 and in_tissue:
                continue
            spatial_dict[row[0]] = [row[-2], row[-1]]
    except:
        raise UserError("Could not read dataframe `tissue_positions_list.csv`.")
    return spatial_dict


def generate_10x_spatial(
        path_to_img=None, path_to_df=None, path_to_json=None, colors=None,
        adata=None, savepath=None, in_tissue=False, palette=None):
    '''
    in_tissue:
        True: Only show spots that are in the tissue
        False: Show all spots
    '''
    if 'json_dict' in adata.uns:
        json_dict = adata.uns['json_dict']
    else:
        json_dict = _read_verify_10x_json(path_to_json)
    scaling_factor = json_dict['tissue_hires_scalef']  # 0.08250825
    full_d = json_dict['spot_diameter_fullres']
    radius = full_d * scaling_factor / 2

    if 'image' in adata.uns:
        small_img = adata.uns['image'].copy()
    else:
        small_img = _read_verify_10x_image(path_to_img)

    if 'spatial_dict' in adata.uns:
        spatial_dict = adata.uns['spatial_dict']
    else:
        spatial_dict = _read_verify_10x_df(path_to_df, in_tissue)
    owner = np.zeros(small_img.shape[:2]) - 1
    #barcodes = list(adata.obs['barcodes'])
    barcodes = list(adata.obs.index)
    if colors is None:
        colors = np.full(len(barcodes), 1, dtype=int)
    elif colors.dtype == int:
        palr, palb, palg = palette_to_rgb(palette, colors.max())
        colors += 1
        colors = np.stack([palr[colors], palb[colors], palg[colors]], axis=-1)
    elif colors.dtype == float:
        small_img.fill(0)
        small_img = small_img[..., 0]
        small_img = small_img.astype(float)
    for i in range(len(adata.obs)):
        center = np.array(spatial_dict[(barcodes[i])]) * scaling_factor
        center = center.astype('int')
        circle = draw.disk(center, radius)
        small_img[circle[0], circle[1]] = colors[i]
        owner[circle[0], circle[1]] = i

    if savepath is not None:
        matplotlib.image.imsave(savepath, small_img, cmap='magma')

    return small_img, owner.astype(int)
