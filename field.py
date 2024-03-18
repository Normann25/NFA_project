#%%
import matplotlib.pyplot as plt
import pandas as pd
import time
from datetime import datetime
from matplotlib.ticker import FuncFormatter
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import os, sys
#%%
def read_MS_txt(path, parent_path):
    data_dict = {}
    Trafic_before = {}
    Asphalt_peak = {}
    After_AP = {}
    Trafic_after = {}
    dictionaries = [Trafic_before, Asphalt_peak, After_AP, Trafic_after]
    labels = ['Trafic_before', 'Asphalt_peak', 'After_AP', 'Trafic_after']

    
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)
    
    files = os.listdir(path)

    for file in files:
        file_name = file.split('.')[0]
        name = file_name.split('_')[0] + ' ' + file_name.split('_')[1]

        with open(os.path.join(path, file)) as f:
            df = pd.read_table(f, sep = '\t')
            df.columns = ['mass', name]

            if 'TraficBefore' in file:
                Trafic_before[name] = df
            if '_Asphalt' in file:
                Asphalt_peak[name] = df
            if 'AfterAsphalt' in file:
                After_AP[name] = df
            if 'TraficAfter' in file:
                Trafic_after[name] = df

    for i, label in enumerate(dictionaries):
        merged = pd.DataFrame({'mass': []})
        for key in label.keys():
            merged = pd.merge(merged, label[key], on = 'mass', how = 'outer')
        
        for key in merged.keys():
            merged[key] = merged[key].fillna(0)
        
        data_dict[labels[i]] = merged
    
    return data_dict

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

            df['t_base'] = df['t_base'].str.split().str[1] + pd.Timedelta('2 hours')
            # df['t_base'] = (df['t_base'].apply(lambda x: datetime.strptime(x, "%H:%M")) + pd.Timedelta("1 hour")).apply(lambda y: datetime.strftime(y, "%H:%M"))
            df['Time'] = pd.to_timedelta(df['t_base']).astype('timedelta64[s]')     # .str.split().str[1]
            df = df.set_index('Time')

            for key in df.keys()[1:]:
                df[key] = pd.to_numeric(df[key].str.replace(',', '.'), errors='coerce')
 
            # df = df.drop('t_base', axis = 'columns')

        data_dict[name] = df
    
    return data_dict

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

            df['t_base'] = df['t_base'].str.split().str[1] + pd.Timedelta('2 hours')
            # df['t_base'] = (df['t_base'].apply(lambda x: datetime.strptime(x, "%H:%M")) + pd.Timedelta("1 hour")).apply(lambda y: datetime.strftime(y, "%H:%M"))
            df['Time'] = pd.to_timedelta(df['t_base']).astype('timedelta64[s]')     # .str.split().str[1]
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
                        
                        new_df['Time'] = pd.to_timedelta(new_df['Time local (hh:mm:ss)']).astype('timedelta64[s]')  # .str.split().str[1]
                        new_df = new_df.set_index('Time')

                        new_df = new_df.dropna()
                        for key in new_df.keys():
                            if 'BCc' in key:
                                new_df[key] = new_df[key] / 1000

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

def plot_ACSM_BC(ax, df_ACSM, df_BC, acsm_key, n):
    mask = df_ACSM[acsm_key] != 0
    p1, = ax.plot(df_ACSM['Time'][mask], df_ACSM[acsm_key][mask], lw = 1, label = 'Organic carbon', color = 'tab:blue')
    ax2 = ax.twinx()
    p2, = ax2.plot(df_BC.index, df_BC['IR BCc'], lw = 1, label = 'Black carbon', color = 'tab:orange')

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)
    ax2.xaxis.set_major_formatter(formatter)

    ylim = np.array(ax.get_ylim())
    ratio = ylim / np.sum(np.abs(ylim))
    scale = ratio / np.min(np.abs(ratio))
    scale = scale / n
    ax2.set_ylim(np.max(np.abs(ax2.get_ylim())) * scale)

    ax.tick_params(axis = 'y', labelcolor = p1.get_color())
    ax2.tick_params(axis = 'y', labelcolor = p2.get_color())

    ax.legend(frameon = False, fontsize = 8, handles = [p1, p2])

    ax.set_xlabel('Time')
    ax.set_ylabel('Intensity', color = p1.get_color())
    ax2.set_ylabel('Mass concentration / $\mu$g/m$^{3}$', color = p2.get_color())
