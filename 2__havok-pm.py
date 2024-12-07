import numpy as np
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
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
outpath = "./output/2__havok-pm"
if not os.path.exists(outpath):
    os.makedirs(outpath)

# set up data paths
data_path = "./data/pm2_5.csv"

# read in data
df = pd.read_csv(data_path) #  , parse_dates=['dateTime']) #, date_format="%Y-&m-%d %H:%M:%S")

# convert to local timezone so plots make sense
df.dateTime = pd.to_datetime(df['dateTime'], utc=True).dt.tz_convert('US/Central')
df['time_since_start_in_min'] = (df.dateTime - df.dateTime.iloc[0]).dt.total_seconds() / 60

# verify we don't have any missing data
df.dropna(subset=['pm2_5'], inplace=True)

# to start, let's group the data by week and create plots
# that way we can identify any periods with interesting
# spikes or trends

# group by week starting on monday
gdf = df.groupby(pd.Grouper(key='dateTime', freq='W-MON'))

# create path for plots
tspath = os.path.join(outpath, "0__ts-plots")
if not os.path.exists(tspath):
    os.makedirs(tspath)

# create plots for each week
for week, df_i in tqdm(gdf):
    # get formatted week
    wk = week.strftime('%y-%m-%d')    

    fig, ax = plt.subplots(figsize=figsizes['wide'])

    # plot the data
    l,  = ax.plot(df_i.dateTime, df_i.pm2_5, linewidth=1, color=mints_colors[2])

    # turn on major/minor girds
    ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
    ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
    ax.minorticks_on()

    # set limits and ticks
    ax.set_xlim(df_i.dateTime.iloc[0], df_i.dateTime.iloc[-1])
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d/%m/%y"))
    ax.xaxis.set_major_locator(mdates.DayLocator(interval=1))
    ax.xaxis.set_minor_locator(mdates.HourLocator(interval=3))

    ax.set_xlabel("time (local)")
    ax.set_ylabel(r"$\text{PM}_{2.5}$ $\left(\mu \text{g} \cdot \text{m}^{-3} \right)$")
    plt.tight_layout()

    # set axis limits and minor ticks
    plt.savefig(os.path.join(tspath, f"{wk}.png"))
    plt.close()



plt.close()
plt.clf()

# create a polar plot to visualize the the time series
df['angles'] = 2 * np.pi * (df.dateTime.dt.hour * 60 + df.dateTime.dt.minute) / (12*60)

# split into two dataframes for AM/PM
df_am = df.loc[df.dateTime.dt.hour < 12, ['angles', 'pm2_5']]
df_pm = df.loc[df.dateTime.dt.hour >= 12, ['angles', 'pm2_5']]

# df_am.dateTime.dt.hour.unique()
# df_pm.dateTime.dt.hour.unique()

df_am.angles.values.max()
df_pm.angles.values.max()

gdf_am = df_am.groupby('angles', as_index=False)
gdf_pm = df_pm.groupby('angles', as_index=False)




pm_01_am = gdf_am.quantile(0.01)
pm_10_am = gdf_am.quantile(0.10)
pm_25_am = gdf_am.quantile(0.25)
pm_50_am = gdf_am.quantile(0.50)
pm_75_am = gdf_am.quantile(0.75)
pm_90_am = gdf_am.quantile(0.90)
pm_99_am = gdf_am.quantile(0.99)

pm_01_pm = gdf_pm.quantile(0.01)
pm_10_pm = gdf_pm.quantile(0.10)
pm_25_pm = gdf_pm.quantile(0.25)
pm_50_pm = gdf_pm.quantile(0.50)
pm_75_pm = gdf_pm.quantile(0.75)
pm_90_pm = gdf_pm.quantile(0.90)
pm_99_pm = gdf_pm.quantile(0.99)


# https://coolors.co/palettes/popular/5%20colors
q_colors = [
    "#264653",
    "#287271",
    "#2a9d8f",
    "#8ab17d",
    "#e9c46a",
    "#f4a261",
    "#e76f51",
]


plt.close()
plt.clf()

# fig, ax = plt.subplots(figsize=figsizes['default'], subplot_kw={'projection': 'polar'})
fig, axs = plt.subplots(1, 2, figsize=figsizes['wide'], subplot_kw={'projection': 'polar'})
plt.subplots_adjust(hspace=0.25)

l_99, = axs[0].plot(pm_99_am.angles, pm_99_am.pm2_5, linewidth=1.0, color=q_colors[0], zorder=3, label="99%")
l_90, = axs[0].plot(pm_90_am.angles, pm_90_am.pm2_5, linewidth=1.0, color=q_colors[1], zorder=3, label="90%")
l_75, = axs[0].plot(pm_75_am.angles, pm_75_am.pm2_5, linewidth=1.0, color=q_colors[2], zorder=3, label="75%")
l_50, = axs[0].plot(pm_50_am.angles, pm_50_am.pm2_5, linewidth=1.0, color=q_colors[3], zorder=3, label="50%")
l_25, = axs[0].plot(pm_25_am.angles, pm_25_am.pm2_5, linewidth=1.0, color=q_colors[4], zorder=3, label="25%")
l_10, = axs[0].plot(pm_10_am.angles, pm_10_am.pm2_5, linewidth=1.0, color=q_colors[5], zorder=3, label="10%")
l_01, = axs[0].plot(pm_01_am.angles, pm_01_am.pm2_5, linewidth=1.0, color=q_colors[6], zorder=3, label="1%")

axs[1].plot(pm_99_pm.angles, pm_99_pm.pm2_5, linewidth=1.0, color=q_colors[0], zorder=3,)
axs[1].plot(pm_90_pm.angles, pm_90_pm.pm2_5, linewidth=1.0, color=q_colors[1], zorder=3,)
axs[1].plot(pm_75_pm.angles, pm_75_pm.pm2_5, linewidth=1.0, color=q_colors[2], zorder=3,)
axs[1].plot(pm_50_pm.angles, pm_50_pm.pm2_5, linewidth=1.0, color=q_colors[3], zorder=3,)
axs[1].plot(pm_25_pm.angles, pm_25_pm.pm2_5, linewidth=1.0, color=q_colors[4], zorder=3,)
axs[1].plot(pm_10_pm.angles, pm_10_pm.pm2_5, linewidth=1.0, color=q_colors[5], zorder=3,)
axs[1].plot(pm_01_pm.angles, pm_01_pm.pm2_5, linewidth=1.0, color=q_colors[6], zorder=3,)


ticks = ['12', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11']
axs[0].set_xticks([(2 * np.pi * i) / 12 for i in range(12)])
axs[1].set_xticks([(2 * np.pi * i) / 12 for i in range(12)])
axs[0].set_xticklabels(ticks)
axs[1].set_xticklabels(ticks)

axs[0].set_ylim(0, 50)
axs[1].set_ylim(0, 50)
axs[0].set_yticks(range(0, 50, 10))
axs[1].set_yticks(range(0, 50, 10))

axs[0].set_theta_direction(-1)
axs[1].set_theta_direction(-1)
axs[0].set_theta_offset(np.pi/2)
axs[1].set_theta_offset(np.pi/2)

axs[0].tick_params(axis='x', pad=-5)
axs[1].tick_params(axis='x', pad=-5)
axs[0].tick_params(axis='x', pad=-5)
axs[1].tick_params(axis='x', pad=-5)

axs[0].set_xlabel("AM")
axs[1].set_xlabel("PM")

leg = axs[0].legend(handles=[l_01, l_10, l_25, l_50, l_75, l_90, l_99], bbox_to_anchor=(-0.3, 0.5), loc="center left", frameon=False, borderaxespad=0, fontsize=10, handletextpad=0.2, handlelength=0.8)

fig.suptitle(r"${PM}_{2.5}$ Quantiles by Hour $\left( \mu \text{g}\cdot \text{m}^{-3} \right)$", y=1.05)

plt.savefig(os.path.join(outpath, "1__pm-radial-dist.png"))



# do the same but on a 24-hour clock.

df['angles'] = 2 * np.pi * (df.dateTime.dt.hour * 60 + df.dateTime.dt.minute) / (24*60)
gdf_24 = df[['angles', 'pm2_5']].groupby('angles', as_index=False)

pm_01 = gdf_24.quantile(0.01)
pm_10 = gdf_24.quantile(0.10)
pm_25 = gdf_24.quantile(0.25)
pm_50 = gdf_24.quantile(0.50)
pm_75 = gdf_24.quantile(0.75)
pm_90 = gdf_24.quantile(0.90)
pm_99 = gdf_24.quantile(0.99)

fig, ax = plt.subplots(figsize=figsizes['wide'], subplot_kw={'projection': 'polar'})

l_99, = ax.plot(pm_99.angles, pm_99.pm2_5, linewidth=1.0, color=q_colors[0], zorder=3, label="99%")
l_90, = ax.plot(pm_90.angles, pm_90.pm2_5, linewidth=1.0, color=q_colors[1], zorder=3, label="90%")
l_75, = ax.plot(pm_75.angles, pm_75.pm2_5, linewidth=1.0, color=q_colors[2], zorder=3, label="75%")
l_50, = ax.plot(pm_50.angles, pm_50.pm2_5, linewidth=1.0, color=q_colors[3], zorder=3, label="50%")
l_25, = ax.plot(pm_25.angles, pm_25.pm2_5, linewidth=1.0, color=q_colors[4], zorder=3, label="25%")
l_10, = ax.plot(pm_10.angles, pm_10.pm2_5, linewidth=1.0, color=q_colors[5], zorder=3, label="10%")
l_01, = ax.plot(pm_01.angles, pm_01.pm2_5, linewidth=1.0, color=q_colors[6], zorder=3, label="1%")



ticks = ['24'] + [f"{i+1}" for i in range(23)]
ax.set_xticks([(2 * np.pi * i) / 24 for i in range(24)])
ax.set_xticklabels(ticks)

ax.set_ylim(0, 50)
ax.set_yticks(range(0, 50, 10))

ax.set_theta_direction(-1)
ax.set_theta_offset(np.pi/2)

ax.tick_params(axis='x', pad=-5)
ax.tick_params(axis='x', pad=-5)

ax.set_title(r"${PM}_{2.5}$ Quantiles by Hour $\left( \mu \text{g}\cdot \text{m}^{-3} \right)$")

leg = ax.legend(handles=[l_01, l_10, l_25, l_50, l_75, l_90, l_99], bbox_to_anchor=(-0.3, 0.5), loc="center left", frameon=False, borderaxespad=0, fontsize=10, handletextpad=0.2, handlelength=0.8)

plt.savefig(os.path.join(outpath, "1b_pm-radial-dist-24.png"))

plt.clf()
plt.close()





# compute indices for skip and use to group the df
# into groups of continuosu measurements
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


print(gdf.size())


zs = []
date_times = []
ts = []

n_emb_max = 200

for gp, df in gdf: 
    if len(df) > n_emb_max + 10:
        zs.append(df.pm2_5.values)
        date_times.append(df.dateTime.values)
        ts.append(df.time_since_start_in_min.values)


# perform parameter sweep
# n_embeddings = range(10, n_emb_max + 10, 10)
# rs_model = range(3,51)
# ns_control = range(1,26)


n_embeddings = range(30, n_emb_max + 10, 10)
rs_model = range(3,26)
ns_control = range(1,11)


# we must have that r_model + n_control <= n_embedding
valid_params = []
for n_emb in n_embeddings:
    for r_mod in rs_model: 
        for n_con in ns_control: 
            if (n_con < r_mod) and (r_mod + n_con <= n_emb):
                valid_params.append((n_emb, r_mod, n_con))


# split off final chunk for independent testing

zs_train = zs[:-1]
ts_train = ts[:-1]
date_times_train = date_times[:-1]

z_test = zs[-1]
t_test = ts[-1]
date_time_test = date_times[-1]

eval_res = []

i = 0
for params in tqdm(valid_params):
    n_embedding, r_model, n_control = params

    # n_embedding, r_model, n_control = (30, 6, 5)

    zs_x, zs_pred, ts_x, U, s, Vs_out, A, B, fvals = eval_havok_multi(
                                                                      zs_train,
                                                                      ts_train,
                                                                      n_embedding,
                                                                      r_model,
                                                                      n_control
                                                                     )

    # integrate the model for the test set
    z_x_test, z_pred_test, t_x_test = integrate_havok(
                                               z_test,
                                               t_test,
                                               n_embedding,
                                               r_model,
                                               n_control,
                                               A, B, U, s
                                           )


    # compute rmse and mae on each separate time series

    rmse_train = [np.sqrt(np.mean(zs_pred[i] - zs_x[i])**2) for i in range(len(zs_x))]
    mae_train = [np.mean(np.abs(zs_pred[i] - zs_x[i])) for i in range(len(zs_x))]

    rmse_test = np.sqrt(np.mean(z_pred_test - z_x_test)**2)
    mae_test = np.mean(np.abs(z_pred_test - z_x_test))

    res_dict = {}
    res_dict['n_embedding'] = n_embedding
    res_dict['r_model'] = r_model
    res_dict['n_control'] = n_control
    res_dict['rmse_train_mean'] = np.mean(rmse_train)
    res_dict['rmse_train_std'] = np.std(rmse_train)
    res_dict['mae_train_mean'] = np.mean(mae_train)
    res_dict['mae_train_std'] = np.std(mae_train)

    res_dict['rmse_test'] = rmse_test
    res_dict['mae_test'] = mae_test

    eval_res.append(res_dict)

    if i % 100 == 0: 
        # save the dataframe        
        df_res = pd.DataFrame(eval_res)
        df_res = df_res.sort_values(by="rmse_train_mean", ascending=True)
        df_res.to_csv(os.path.join(outpath, "param_sweep.csv"))


    i += 1

df_res = pd.DataFrame(eval_res)
df_res = df_res.sort_values(by="rmse_train_mean", ascending=True)
df_res.to_csv(os.path.join(outpath, "param_sweep.csv"))

print(df_res.head())

