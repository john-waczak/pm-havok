import numpy as np
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib import rcParams
from matplotlib.collections import LineCollection
from cycler import cycler
import warnings



# Set the font files
def matplotlib_defaults(font='montserrat', background_transparent=False):

    fonts_path = "./font-files"
    font_files = font_manager.findSystemFonts(fontpaths=fonts_path)
    for font_file in font_files:
        # print(font_file)
        font_manager.fontManager.addfont(font_file)

    rcParams.update({
                        'figure.dpi': 300,             # default for high quality
                        'savefig.bbox': 'tight',       # tight, standard
                        'savefig.transparent': background_transparent,
                        'font.family': 'sans-serif',
                        'font.sans-serif': font,
                        'mathtext.fontset': 'custom',
                        'mathtext.rm': font,
                        'font.size': 12,
                        'axes.labelsize': 12,          # MDPI uses 10pt font
                        'xtick.labelsize': 10,          # make tick labels slightly smaller
                        'ytick.labelsize': 10,
                        'axes.titlesize': 15,          # make title slightly larger
                        'xtick.major.size': 7,
                        'xtick.major.width': 0.7,
                        'ytick.major.width': 0.7,
                        'ytick.major.size': 7,
                        'xtick.minor.visible': True,
                        'ytick.minor.visible': True,
                        'xtick.minor.size': 3,
                        'ytick.minor.size': 3,
                        'xtick.minor.width': 0.4,
                        'ytick.minor.width': 0.4,
                        'axes.titleweight': 'bold',
                        'axes.prop_cycle': cycler('color', ['#3cd184', '#f97171', '#1e81b0', '#66beb2', '#f99192', '#8ad6cc', '#3d6647', '#000080']),
                        'image.cmap': 'plasma'
                    })    

figsizes = {
    "default" : (4,3),
    "wide" : (7, 3),
}

mints_colors = [
    '#3cd184',
    '#f97171', 
    '#1e81b0', 
    '#66beb2', 
    '#f99192', 
    '#8ad6cc', 
    '#3d6647', 
    '#000080'
]



def colored_line(x, y, c, ax, **lc_kwargs):
    """
    Plot a line with a color specified along the line by a third value.

    It does this by creating a collection of line segments. Each line segment is
    made up of two straight lines each connecting the current (x, y) point to the
    midpoints of the lines connecting the current point with its two neighbors.
    This creates a smooth line with no gaps between the line segments.

    Parameters
    ----------
    x, y : array-like
        The horizontal and vertical coordinates of the data points.
    c : array-like
        The color values, which should be the same size as x and y.
    ax : Axes
        Axis object on which to plot the colored line.
    **lc_kwargs
        Any additional arguments to pass to matplotlib.collections.LineCollection
        constructor. This should not include the array keyword argument because
        that is set to the color argument. If provided, it will be overridden.

    Returns
    -------
    matplotlib.collections.LineCollection
        The generated line collection representing the colored line.
    """
    if "array" in lc_kwargs:
        warnings.warn('The provided "array" keyword argument will be overridden')

    # Default the capstyle to butt so that the line segments smoothly line up
    default_kwargs = {"capstyle": "butt"}
    default_kwargs.update(lc_kwargs)

    # Compute the midpoints of the line segments. Include the first and last points
    # twice so we don't need any special syntax later to handle them.
    x = np.asarray(x)
    y = np.asarray(y)
    x_midpts = np.hstack((x[0], 0.5 * (x[1:] + x[:-1]), x[-1]))
    y_midpts = np.hstack((y[0], 0.5 * (y[1:] + y[:-1]), y[-1]))

    # Determine the start, middle, and end coordinate pair of each line segment.
    # Use the reshape to add an extra dimension so each pair of points is in its
    # own list. Then concatenate them to create:
    # [
    #   [(x1_start, y1_start), (x1_mid, y1_mid), (x1_end, y1_end)],
    #   [(x2_start, y2_start), (x2_mid, y2_mid), (x2_end, y2_end)],
    #   ...
    # ]
    coord_start = np.column_stack((x_midpts[:-1], y_midpts[:-1]))[:, np.newaxis, :]
    coord_mid = np.column_stack((x, y))[:, np.newaxis, :]
    coord_end = np.column_stack((x_midpts[1:], y_midpts[1:]))[:, np.newaxis, :]
    segments = np.concatenate((coord_start, coord_mid, coord_end), axis=1)

    lc = LineCollection(segments, **default_kwargs)
    lc.set_array(c)  # set the colors of each segment

    return ax.add_collection(lc)

