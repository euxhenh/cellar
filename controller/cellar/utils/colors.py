import numpy as np
import plotly.colors as pc


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


def get_col_i(pal, i):
    return pal[i % len(pal)]


def palette_to_rgb(
        palette=None, maxval=None, black="#000000", white="#FFFFFF"):
    # Colors when labels are provided
    if palette is None:
        palette = PALETTE
    if maxval is not None:
        while len(palette) - 1 < maxval:
            palette += palette
    # add 0 for black and -1 for white
    palette = [black] + palette + [white]
    palr = np.array([int(i[1:3], 16) for i in palette])  # hex string to int
    palb = np.array([int(i[3:5], 16) for i in palette])
    palg = np.array([int(i[5:], 16) for i in palette])
    return palr, palb, palg


def interpolate_grayimage(im, colorscale='magma'):
    """
    Maps a 2D image to a 3D image using the given colormap.
    We use this manually since plotly is very slow when rendering
    large 2D images using imshow as in the case of CODEX data.
    """
    assert im.ndim == 2
    colorscale = pc.get_colorscale(colorscale)
    colors = pc.colorscale_to_colors(colorscale)
    colors = pc.validate_colors(colors, colortype="tuple")

    # First we bin each value in im so that we know which colors
    # to use for interpolation
    unit = (im - im.min()) / (im.max() - im.min())
    binned = (unit * (len(colors) - 1)) - 0.5
    binned = np.round(binned).astype(int).clip(min=0)
    # empty 3d im
    im3d = np.zeros((*im.shape, 3))
    for i in range(len(colors) - 1):  # Small for loop so it's ok
        r1, g1, b1 = colors[i]
        r2, g2, b2 = colors[i + 1]
        idx = np.where(binned == i)
        minv, maxv = im[idx].min(), im[idx].max()
        intermeds = (im[idx] - minv) / (maxv - minv)
        im3d[idx] = np.stack([
            (r2 - r1) * intermeds + r1,
            (g2 - g1) * intermeds + g1,
            (b2 - b1) * intermeds + b1
        ], axis=-1)
    return im3d
