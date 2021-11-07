import numpy as np


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
