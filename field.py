#%%
import matplotlib.pyplot as plt
import pandas as pd
import time
from matplotlib.ticker import FuncFormatter
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os, sys
#%%
def read_data(path, parent_path):
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)

    files = os.listdir(path)
    data_dict = {}

    for file in files:
        name = file.split('.')[0]
        with open(os.path.join(path, file)) as f:
            df = pd.read_csv(f, sep = ';', decimal = ',')
            df = df.dropna()

            df['Time'] = pd.to_timedelta(df['t_base'].str.split().str[1]).astype('timedelta64[s]')
            df = df.set_index('Time')

            for key in df.keys()[1:]:
                df[key] = pd.to_numeric(df[key].str.replace(',', '.'), errors='coerce')
 
            # df = df.drop('t_base', axis = 'columns')

        data_dict[name] = df
    
    return data_dict
#%%
def read_csv(path, parent_path):
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)
    
    files = os.listdir(path)
    data_dict = {}

    for file in files:
        name = file.split('.')[0]
        with open(os.path.join(path, file)) as f:
            df = pd.read_csv(f)
            df = df.dropna()

            df['Time'] = pd.to_timedelta(df['t_base'].str.split().str[1]).astype('timedelta64[s]')
            df = df.set_index('Time')

        data_dict[name] = df
    
    return data_dict

def read_csv_BC(path, parent_path, bc_station):
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)
    
    files = os.listdir(path)
    data_dict = {}

    for file in files:
        if '.csv' in file:
            for st in bc_station:
                if st in file:
                    name = file.split('.')[0]
                    name = name.split('_')[-1]
                    with open(os.path.join(path, file)) as f:
                        df = pd.read_csv(f)

                        new_df = pd.DataFrame()
                        columns = ['Time local (hh:mm:ss)', 'Sample temp (C)', 'Sample RH (%)', 'UV BCc', 'Blue BCc', 'Green BCc', 'Red BCc', 'IR BCc']
                        for col in columns:
                            new_df[col] = df[col]
                        new_df = new_df.dropna()

                    data_dict[name] = new_df
    
    return data_dict 

def plot_test(ax, df):
    df['Sum'] = df[df.keys()[2:]].sum(axis='columns')  

    ax.plot(df.index, df['Sum'], lw = 1, label = 'Summed columns')

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon = False, fontsize = 8)
    ax.set(ylabel = 'Intensity', xlabel = 'Time')

def plot_overview(ax, df, ncol):
    for key in df.keys()[1:]:
        ax.plot(df.index, df[key], lw = 1, label = key)

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon = False, fontsize = 8, ncol = ncol)
    ax.set(ylabel = 'Intensity', xlabel = 'Time')

def plot_PAH_ACSM(ax, df, ncol):
    for key in df.keys()[1:]:
        if 'm' in key:
            ax.plot(df.index, df[key], lw = 1, label = key)

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon = False, fontsize = 8, ncol = ncol)
    ax.set(ylabel = 'Intensity', xlabel = 'Time')
