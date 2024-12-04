import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 

import matplotlib.dates as mdates
from matplotlib.colors import LogNorm

from matplotlib_defaults import matplotlib_defaults, figsizes, mints_colors

from ips7100 import pm_colors, pm_labels_units, pc_cols, pm_cols, bin_edges, correct_ips7100_for_humidity

matplotlib_defaults(font='palatino')


# set up plotting paths
figpath = "./figures/0__humidity-correction"
if not os.path.exists(figpath):
    os.makedirs(figpath)

# set up data paths
data_path = "./data"
bam_path = os.path.join(data_path, "bam.csv")
ips_path = os.path.join(data_path, "ips7100-hc.csv")




# read in data
df = pd.read_csv(ips_path) #, parse_dates=['dateTime']) #, date_format="%Y-&m-%d %H:%M:%S")
df_bam = pd.read_csv(bam_path) #, parse_dates=['dateTime'])

# parse dateTime column to dateTime
# df.dateTime = pd.to_datetime(df['dateTime']).dt.tz_localize('UTC')
# df_bam.dateTime = pd.to_datetime(df_bam['dateTime']).dt.tz_localize('UTC')
df.dateTime = pd.to_datetime(df['dateTime']).dt.tz_localize('UTC').dt.tz_convert('US/Central')
df_bam.dateTime = pd.to_datetime(df_bam['dateTime']).dt.tz_localize('UTC').dt.tz_convert('US/Central')


sum(df.dateTime.notnull())

df.describe()
df.head()

type(df.loc[1, 'dateTime'])

df.loc[1, 'dateTime']
df_bam.loc[1, 'dateTime']


# visualize the PM time series data

fig, ax = plt.subplots(figsize=figsizes['wide'])

l_10,  = ax.plot(df.dateTime, df.pm10_0, linewidth=1, label=r"$\text{PM}_{10}$")
l_2_5, = ax.plot(df.dateTime, df.pm2_5, linewidth=1, label=r"$\text{PM}_{2.5}$")
l_1,   = ax.plot(df.dateTime, df.pm1_0, linewidth=1, label=r"$\text{PM}_{1}$")

# turn on major/minor girds
ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

# set axis limits and minor ticks
ax.set_xlim(df.dateTime.iloc[0], df.dateTime.iloc[-1])
ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))
ax.set_ylim(-2, 102)

ax.set_xlabel("time (local)")
ax.set_ylabel(r"$\text{Concentration}$ $\left(\mu \text{g} \cdot \text{m}^{-3} \right)$")

# create legend on top of plot
plt.legend(handles=[l_1, l_2_5, l_10], bbox_to_anchor=(0, 1.015, 1, 0.15), loc="upper center", frameon=False, borderaxespad=0, ncol=3)

plt.tight_layout()

# plt.show()

plt.savefig(os.path.join(figpath, "1__pm-timeseries.png"))



# Visualize particle counts heatmap

# extract PC into a matrix
PC = df.loc[:, pc_cols].values

print(PC.shape)


fig, ax = plt.subplots(figsize=figsizes["wide"])

dt = df.dateTime[1] - df.dateTime[0]
xbins = [d for d in df.dateTime - 0.5*dt]
xbins.append(xbins[-1] + dt)
xbins = pd.to_datetime(xbins)

xbins

ax.xaxis.set_major_formatter(mdates.DateFormatter("%m/%y"))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=1))
ax.xaxis.set_minor_locator(mdates.DayLocator(interval=7))

# pcm = ax.pcolormesh(xbins, bin_edges, np.log10(PC.T))
pcm = ax.pcolormesh(xbins, bin_edges, PC.T, norm=LogNorm(vmin=1, vmax=1e6))

ax.hlines(bin_edges[1:-1], xbins[0], xbins[-1], color='#DDDDDD', linewidth=0.5, linestyle="--")
ax.set_xlabel("time (local)")
ax.set_ylabel(r"Particle Diameter ($\mu$m)")

cb = plt.colorbar(pcm, ax=ax, pad=0.015, extend='both')
# cb.ax.set_ylabel(r"Particle Counts (#$\cdot \text{L}^{-1}$)")
cb.ax.set_ylabel(r"Particle Counts (#$\cdot \text{L}^{-1}$)", size=10)
# cb.ax.tick_params(axis='y', which='major', width=1)
# cb.ax.tick_params(axis='y', which='minor', width=0.75)

cb.ax.tick_params()
cb.ax.minorticks_on()
# ax.tick_params(axis='both', which='major', width=1)
# ax.tick_params(axis='both', which='minor', width=0.75)

plt.tight_layout()
plt.savefig(os.path.join(figpath, "2__heatmap.png"))


# now we can perform the humidity correction
correct_ips7100_for_humidity(df, kappa=0.62, eff=0.35, rh_col='humidity')

df.humidity

# plot the original vs the corrected PM values

fig, axs = plt.subplots(2, 1, height_ratios=[1, 4], sharex=True, figsize=figsizes['wide'])
plt.subplots_adjust(hspace=0.25)

l_orig,  = axs[1].plot(df.dateTime, df.pm2_5, linewidth=0.8, label=r"Original")
l_cor, = axs[1].plot(df.dateTime, df.pm2_5_cor, linewidth=0.8, label=r"Humidity Corrected")

plt.legend(handles=[l_orig, l_cor], bbox_to_anchor=(0, 1.04, 1, 0.01), loc="center right", frameon=False, borderaxespad=0, ncol=3, fontsize=8)

l_hum, = axs[0].plot(df.dateTime, df.humidity, linewidth=0.8, color=mints_colors[3])

# turn on major/minor girds
axs[0].grid(which='major', color='#DDDDDD', linewidth=0.8)
axs[0].grid(which='minor', color='#DDDDDD', linewidth=0.5)
axs[0].minorticks_on()
axs[0].tick_params(axis='x', which='major', length=3)
axs[0].tick_params(axis='x', which='minor', length=2)

axs[1].grid(which='major', color='#DDDDDD', linewidth=0.8)
axs[1].grid(which='minor', color='#DDDDDD', linewidth=0.5)
axs[1].minorticks_on()

# set axis limits and minor ticks
axs[1].set_xlim(df.dateTime.iloc[0], df.dateTime.iloc[-1])
axs[1].xaxis.set_major_formatter(mdates.DateFormatter("%m/%y"))
axs[1].xaxis.set_major_locator(mdates.MonthLocator(interval=1))
axs[1].xaxis.set_minor_locator(mdates.DayLocator(interval=7))
axs[1].set_ylim(-2, 102)
axs[0].set_ylim(-5, 105)

axs[1].set_xlabel("time (local)")
axs[1].set_ylabel(r"$\text{PM}_{2.5}$ $\left(\mu \text{g} \cdot \text{m}^{-3} \right)$")
axs[0].set_ylabel("RH (%)")

plt.savefig(os.path.join(figpath, "3__pm-humidity-corrected.png"))

# plt.show()





# create a scatter diagram comparing the corrected IPS7100 to BAM data

# join the two dataframes

# create new column for bam pm 2.5
#df_pm['bam2_5'] = np.nan

# set the index to datetime
df.set_index(df.dateTime, inplace=True)
df_bam.set_index(df_bam.dateTime, inplace=True)

# merge on the datetime column
df_pm = df.merge(df_bam, left_index=True, right_index=True, how='inner')


# now we can make a scatter diagram


fig, ax = plt.subplots(figsize=figsizes['default'])

ax.grid(which='major', color='#DDDDDD', linewidth=0.8)
ax.grid(which='minor', color='#DDDDDD', linewidth=0.5)
ax.minorticks_on()

s_orig = ax.scatter(df_pm.pm2_5BAM, df_pm.pm2_5, alpha=0.8, s=1, label="Original", zorder=3)
s_hc_p = ax.scatter(df_pm.pm2_5BAM, df_pm.pm2_5HC, alpha=0.8, s=1, label="Prabuddha", zorder=3)
s_hc = ax.scatter(df_pm.pm2_5BAM, df_pm.pm2_5_cor, alpha=0.8, s=1, label="Di Antonio", zorder=3)
ll, = ax.plot([-1000, 1000], [-1000, 1000], color='0.5', linewidth=0.5, label="1:1", zorder=4)

low = min(df_pm.pm2_5BAM.quantile(q=0.01), min(df_pm[['pm2_5', 'pm2_5_cor', 'pm2_5HC']].quantile(q=0.01)), -1)
high = max(df_pm.pm2_5BAM.quantile(q=0.99), max(df_pm[['pm2_5', 'pm2_5_cor', 'pm2_5HC']].quantile(q=0.99)), 1)

low = low - 0.1*abs(low)
high = high + 0.1*abs(high)

# add the legend
leg = ax.legend(handles=[s_orig, s_hc_p, s_hc, ll], bbox_to_anchor=(0, 1.01, 1.0, 0.1), loc="upper center", frameon=False, borderaxespad=0, ncol=4, columnspacing=0.5, fontsize=8, handletextpad=0.0, handlelength=0.8)

ax.set_xlim(low, high)
ax.set_ylim(low, high)

ax.set_xlabel(r"BAM  $\text{PM}_{2.5}$ ($\mu \text{g}\cdot\text{m}^{-3}$)")
ax.set_ylabel(r"IPS7100  $\text{PM}_{2.5}$ ($\mu \text{g}\cdot\text{m}^{-3}$)")
plt.tight_layout()

plt.savefig(os.path.join(figpath, "4__hc-bam-scatter.png"))

# plt.show()


# output the corrected dataframe with just PM 2.5 
df_out = df[['dateTime', 'pm2_5_cor']].copy()
df_out.rename(columns={'pm2_5_cor': "pm2_5"}, inplace=True)
df_out.columns
df_out.to_csv("./data/pm2_5.csv")

