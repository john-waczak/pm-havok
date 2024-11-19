from matplotlib import font_manager
from matplotlib import rcParams
from cycler import cycler

# Set the font files
def matplotlib_defaults(font='montserrat'):

    fonts_path = "./font-files"
    font_files = font_manager.findSystemFonts(fontpaths=fonts_path)
    for font_file in font_files:
        print(font_file)
        font_manager.fontManager.addfont(font_file)

    rcParams.update({
                        'font.family': font,
                        'font.size': 20,
                        'axes.labelsize': 18,
                        'xtick.labelsize': 15,
                        'ytick.labelsize': 15,
                        'xtick.major.size': 7,
                        'xtick.major.width': 1.5,
                        'ytick.major.width': 1.5,
                        'ytick.major.size': 7,
                        'xtick.minor.visible': True,
                        'ytick.minor.visible': True,
                        'xtick.minor.size': 3,
                        'ytick.minor.size': 3,
                        'xtick.minor.width': 0.75,
                        'ytick.minor.width': 0.75,
                        'axes.titleweight': 'bold',
                        'axes.titlesize': 20,
                        'axes.prop_cycle': cycler('color', ['#3cd184', '#f97171', '#1e81b0', '#66beb2', '#f99192', '#8ad6cc', '#3d6647', '#000080']),
                        'image.cmap': 'plasma'
                    })    


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


