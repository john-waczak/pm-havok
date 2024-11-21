from matplotlib import font_manager
from matplotlib import rcParams
from cycler import cycler

# Set the font files
def matplotlib_defaults(font='montserrat', background_transparent=False):

    fonts_path = "./font-files"
    font_files = font_manager.findSystemFonts(fontpaths=fonts_path)
    for font_file in font_files:
        print(font_file)
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


# matplotlib_defaults(font='montserrat')

# # Create a simple plot with LaTeX text
# import matplotlib.pyplot as plt

# plt.figure()
# plt.plot([0, 1], [0, 1])
# plt.xlabel("PM")
# plt.ylabel(r"$\text{PM}_{2.5}$ $\left( \mu \text{g}\cdot\text{m}^{-3} \right)$")

# # Show the plot
# plt.show()


# matplotlib_defaults(font='palatino')



# import numpy as np
# import matplotlib.pyplot as plt

# x = np.linspace(0, 2*np.pi, 500)
# y = np.sin(x)

# fig, ax = plt.subplots(1, 1)
# ax.plot(x, y)
# ax.set_xlabel("x label")
# ax.set_ylabel("y label")
# ax.set_title("Title")
# ax.grid(which='major', color='#DDDDDD', linewidth=1.5)
# ax.grid(which='minor', color='#DDDDDD', linewidth=0.75)
# ax.minorticks_on()
# plt.show()


