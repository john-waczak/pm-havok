import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from matplotlib.colors import LogNorm, LinearSegmentedColormap

import pandas as pd
import os
from tqdm import tqdm # for progress bars on iterators

# matplotlib customizations
from matplotlib_defaults import matplotlib_defaults, figsizes, mints_colors


# handling ips7100 data
# from ips7100 import pm_colors, pm_labels_units, pc_cols, pm_cols, bin_edges, correct_ips7100_for_humidity

# load in our HAVOK functions
from havok import eval_havok_multi, integrate_havok


# set the matplotlib defaults and use font
# for MDPI article
matplotlib_defaults(font='palatino')

# set up plotting paths
outpath = "./output/3__train-best"
if not os.path.exists(outpath):
    os.makedirs(outpath)

# set up data paths
data_path = "./data/pm2_5.csv"
res_path = "./output/2__havok-pm/param_sweep.csv"

# read in data
df = pd.read_csv(data_path) #  , parse_dates=['dateTime']) #, date_format="%Y-&m-%d %H:%M:%S")
df_res = pd.read_csv(res_path) #  , parse_dates=['dateTime']) #, date_format="%Y-&m-%d %H:%M:%S")


df_res = df_res.sort_values(by="rmse_train_mean", ascending=True)
print(df_res.iloc[0])

# Now let's train the model with te winning params
n_embedding = 160
r_model = 3
n_control = 2

# set up dataset
df.dateTime = pd.to_datetime(df['dateTime'], utc=True).dt.tz_convert('US/Central')
df['time_since_start_in_min'] = (df.dateTime - df.dateTime.iloc[0]).dt.total_seconds() / 60
df.dropna(subset=['pm2_5'], inplace=True)

dt = df.dateTime[1] - df.dateTime[0]
is_skip = [False] + [next - cur != dt for cur, next in zip(df.dateTime, df.dateTime[1:])]

group = np.ones(len(is_skip), dtype=int)
val = 0
for i in range(len(df)):
    if is_skip[i]:
        val += 1

    group[i] = val

df['group'] = group
gdf = df.groupby('group')

zs = []
date_times = []
ts = []
n_emb_max = 200

for gp, df in gdf: 
    if len(df) > n_emb_max + 10:
        zs.append(df.pm2_5.values)
        date_times.append(df.dateTime.values)
        ts.append(df.time_since_start_in_min.values)

zs_train = zs[:-1]
ts_train = ts[:-1]
date_times_train = date_times[:-1]

z_test = zs[-1]
t_test = ts[-1]
date_time_test = date_times[-1]

print('Fitting HAVOK model for specified parameters')
zs_x, zs_pred, ts_x, U, s, Vs_out, A, B, fvals = eval_havok_multi(
    zs_train,
    ts_train,
    n_embedding,
    r_model,
    n_control,
)

print('Integrating HAVOK model for test set.')
z_x_test, z_pred_test, t_x_test = integrate_havok(
    z_test,
    t_test,
    n_embedding,
    r_model,
    n_control,
    A, B, U, s
)


# plot the HAVOK operator matrices
fig, axs = plt.subplots(1, 2, gridspec_kw={'width_ratios': [r_model, n_control]}, figsize=figsizes['default'])


axs[0].set_xticks([])
axs[1].set_xticks([])
axs[0].set_yticks([])
axs[1].set_yticks([])

vmin = -0.005
vmax = 0.005
cmap = LinearSegmentedColormap.from_list('list', [mints_colors[2], '#FFFFFF', mints_colors[1]])

im_A = axs[0].imshow(A, cmap=cmap, vmin=vmin, vmax=vmax)
im_B = axs[1].imshow(B, cmap=cmap, vmin=vmin, vmax=vmax)


axs[0].set_title('A')
axs[1].set_title('B')

fig.subplots_adjust(right=0.85, wspace=0.1)
cbar_ax = fig.add_axes([0.875, 0.25, 0.02, 0.5])
# cbar_ax.set_xticks(range(vmin, vmax, 5))
cb = fig.colorbar(im_A, cax=cbar_ax, extend='both')

plt.savefig(os.path.join(outpath, '1__A-B-operator-heatmap.png'))


