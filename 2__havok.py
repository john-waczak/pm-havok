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
figpath = "./figures/2__havok"
if not os.path.exists(figpath):
    os.makedirs(figpath)

# set up data paths
data_path = "./data/pm2_5.csv"

# read in data
df = pd.read_csv(data_path) #, parse_dates=['dateTime']) #, date_format="%Y-&m-%d %H:%M:%S")

# convert to local timezone so plots make sense
df.dateTime = pd.to_datetime(df['dateTime'], utc=True).dt.tz_convert('US/Central')

# verify we don't have anymore missing data
df.dropna(subset=['pm2_5'], inplace=True)


# compute indices for skip and use to group the df
# into groups of continuosu measurements
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
