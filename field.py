#%%
import matplotlib.pyplot as plt
import pandas as pd
import time
from matplotlib.ticker import FuncFormatter
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os
#%%
def read_data(path):
    files = os.listdir(path)
    data_dict = {}

    for file in files:
        name = file.split('.')[0]
        with open(os.path.join(path, file)) as f:
            df = pd.read_csv(f, sep = ';')
            df = df.dropna()

            df['Time'] = pd.to_timedelta(df['t_base'].str.split().str[1]).astype('timedelta64[s]')
            df = df.set_index('Time')
            # df = df.drop('t_base', axis = 'columns')

        data_dict[name] = df
    
    return data_dict
#%%
def plot_overview(ax, df):
    summed = df[df.keys()[2:]].sum(axis='columns')
    ax.plot(df.index, df['Total ion current'], lw = 1, label = 'Total ion current')
    ax.plot(df.index, summed, lw = 1, label = 'Summed columns')

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon = False, fontsize = 8)
    ax.set(ylabel = 'Intensity', xlabel = 'Time')