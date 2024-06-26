#%%
import matplotlib.pyplot as plt
import pandas as pd
import time
from datetime import datetime
from matplotlib.ticker import FuncFormatter
import numpy as np
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
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

def read_ACMS_txt(path, parent_path, labels_list):
    data_dict = {}

    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)
    
    files = os.listdir(path)
    
    for label in labels_list:
        label_dict = {}
        for file in files:
            if label in file:
                file_name = file.split('.')[0]
                name = file_name.split('_')[0] + ' ' + file_name.split('_')[1]
            
                with open(os.path.join(path, file)) as f:
                    df = pd.read_table(f, sep = '\t')
                    df.columns = ['Time', name]
                    df['Time'] = df['Time'].str.split().str[1] + pd.Timedelta('2 hours')
                    df['Time'] = pd.to_timedelta(df['Time']).astype('timedelta64[s]')     # .str.split().str[1]
                    label_dict[name] = df

                merged = pd.DataFrame({'Time':[]})
                for key in label_dict.keys():
                    merged = pd.merge(merged, label_dict[key], on = 'Time', how = 'outer')
                
                for key in merged.keys():
                    merged[key] = merged[key].fillna(0)
                
                merged.sort_values(by = 'Time', ignore_index = True)

                data_dict[label] = merged

    return data_dict

def read_data(path, parent_path, time_label):
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)

    files = os.listdir(path)
    data_dict = {}

    for file in files:
        file_name = file.split('.')[0]
        name = file_name.split('_')[0] + ' ' + file_name.split('_')[1]
        with open(os.path.join(path, file)) as f:
            df = pd.read_csv(f, sep = ';')
            df = df.dropna()

            df['PAH'] = pd.to_numeric(df['PAH'], errors = 'coerce')

            df[time_label] = df[time_label].str.split().str[1] + pd.Timedelta('2 hours')
            df['Time'] = pd.to_timedelta(df[time_label]).astype('timedelta64[s]')

        data_dict[name] = df
    
    return data_dict

def read_csv(path, parent_path, time_label):
    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)
    
    files = os.listdir(path)
    data_dict = {}

    for file in files:
        file_name = file.split('.')[0]
        name = file_name.split('_')[0] + ' ' + file_name.split('_')[1]
        with open(os.path.join(path, file)) as f:
            df = pd.read_csv(f)
            df = df.dropna()

            df[time_label] = df[time_label].str.split().str[1] + pd.Timedelta('2 hours')
            # df['t_base'] = (df['t_base'].apply(lambda x: datetime.strptime(x, "%H:%M")) + pd.Timedelta("1 hour")).apply(lambda y: datetime.strftime(y, "%H:%M"))
            df['Time'] = pd.to_timedelta(df[time_label]).astype('timedelta64[s]')     # .str.split().str[1]
            # df = df.set_index('Time')

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
                        # new_df = new_df.set_index('Time')

                        new_df = new_df.dropna()
                        for key in new_df.keys():
                            if 'BCc' in key:
                                new_df[key] = new_df[key] / 1000

                    data_dict[name] = new_df
    
    return data_dict 

def plot_overview(ax, df, label, df_keys, ncol):
    for key in df_keys:
        if label in df[key].keys():
            mask = df[key][label] != 0
            ax.plot(df[key]['Time'][mask], df[key][label][mask], lw = 1, label = key)

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.legend(frameon = False, fontsize = 8, ncol = ncol)
    ax.set(ylabel = 'Intensity', xlabel = 'Time')

def plot_PAH_ACSM(ax, data_dict, dict_keys, colors, peak_int):
    for i, dict_key in enumerate(dict_keys):
        df = data_dict[dict_key]
        for j, key in enumerate(df.keys()[1:]):
            ax[i][0].plot(df['Time'][peak_int[i][0]:peak_int[i][1]], df[key][peak_int[i][0]:peak_int[i][1]], lw = 1, color = colors[j])

            if 'm/z' in key:
                ax[i][1].plot(df['Time'][peak_int[i][2]:peak_int[i][3]], df[key][peak_int[i][2]:peak_int[i][3]], lw = 1, color = colors[j])

        formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
        ax[i][0].xaxis.set_major_formatter(formatter)
        ax[i][0].tick_params(axis = 'both', which = 'major', direction = 'out', bottom = True, left = True, labelsize = 8)

        ax[i][0].set_xlabel('Time', fontsize = 8)
        ax[i][0].set_ylabel('PAH$_{est}$ / $\mu$g/m$^{3}$', fontsize = 8)
        ax[i][0].set_title(dict_key, fontsize = 9)

        formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
        ax[i][1].xaxis.set_major_formatter(formatter)
        ax[i][1].tick_params(axis = 'both', which = 'major', direction = 'out', bottom = True, left = True, labelsize = 8)

        ax[i][1].set_xlabel('Time', fontsize = 8)
        ax[i][1].set_ylabel('PAH$_{est}$ / $\mu$g/m$^{3}$', fontsize = 8)
        ax[i][1].set_title(dict_key, fontsize = 9)

        


def plot_PAH_ACSM_wInset(ax, df, key_start, colors, loc, bb2a, height, peak_idx):
    for i, key in enumerate(df.keys()[key_start:]):
        ax.plot(df['Time'], df[key], lw = 1, color = colors[i])

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)

    ax.set(ylabel = 'PAH$_{est}$ conc. / $\mu$g/m$^{3}$', xlabel = 'Time')

    inset_ax = inset_axes(ax,
                        width = "40%", # width = % of parent_bbox
                        height = height, # height : 1 inch
                        loc = loc,
                        bbox_to_anchor = bb2a,
                        bbox_transform = ax.transAxes) # placement in figure

    for i, key in enumerate(df.keys()[key_start+1:]):
        inset_ax.plot(df['Time'][peak_idx[0]:peak_idx[1]], df[key][peak_idx[0]:peak_idx[1]], lw = 1, color = colors[i+1])
    
    inset_ax.set_ylabel(None)
    inset_ax.set_xlabel(None)

    inset_ax.xaxis.set_major_formatter(formatter)

    inset_ax.legend().set_visible(False)

def plot_ACSM_BC(ax, df_ACSM, df_BC, acsm_key, BC_mask, n):
    mask = df_ACSM[acsm_key] != 0
    p1, = ax.plot(df_ACSM['Time'][mask], df_ACSM[acsm_key][mask], lw = 1, label = 'Organic aerosol', color = 'tab:blue')
    ax2 = ax.twinx()
    if 'Vesterbrogade' in acsm_key:
        p2, = ax2.plot(df_BC['Time'][BC_mask], df_BC['IR BCc'][BC_mask], lw = 1, label = 'Black carbon', color = 'k')
        ax2.plot(df_BC['Time'], df_BC['IR BCc'], lw = 1, label = None, color = 'k', ls = '--')
    else:
        p2, = ax2.plot(df_BC['Time'], df_BC['IR BCc'], lw = 1, label = 'Black carbon', color = 'k')

    formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
    ax.xaxis.set_major_formatter(formatter)
    ax.set_xticklabels(ax.get_xticklabels(), size = 8)
    ax2.xaxis.set_major_formatter(formatter)
    ax2.set_xticklabels(ax2.get_xticklabels(), size = 8)

    ylim = np.array(ax.get_ylim())
    ratio = ylim / np.sum(np.abs(ylim))
    scale = ratio / np.min(np.abs(ratio))
    scale = scale / n
    ax2.set_ylim(np.max(np.abs(ax2.get_ylim())) * scale)

    ax.tick_params(axis = 'y', labelcolor = p1.get_color(), labelsize = 8)
    ax2.tick_params(axis = 'y', labelcolor = p2.get_color(), labelsize = 8)

    ax.legend(frameon = False, fontsize = 8, handles = [p1, p2])

    ax.set_xlabel('Time', fontsize = 8)
    ax.set_ylabel('Org / $\mu$g/m$^{3}$', color = p1.get_color(), fontsize = 8)
    ax2.set_ylabel('BC / $\mu$g/m$^{3}$', color = p2.get_color(), fontsize = 8)

def plot_MS(ax, df, key, ttl):
    width = 0.2

    ax.bar(df['mass'], df[key], width)
    
    ax.set(xlabel = 'm/z', ylabel = 'Intensity', title = ttl)

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis = 'both', which = 'major', direction = 'out', bottom = True, left = True, labelsize = 8)
    ax.tick_params(axis = 'both', which = 'minor', direction = 'out', width = 1, length = 2, bottom = True, left = True)
    ax.yaxis.offsetText.set_fontsize(9)

def plot_105_183(ax, acsm_dict, acsm_keys, bc_dict, n, peak_int):
    df1 = acsm_dict['m105']
    df2 = acsm_dict['m183']

    for i, key in enumerate(bc_dict.keys()):
        mask1 = df1[acsm_keys[i]] != 0
        mask2 = df2[acsm_keys[i]] != 0

        x1, y1 = df1['Time'][mask1][peak_int[i][0]:peak_int[i][1]], df1[acsm_keys[i]][mask1][peak_int[i][0]:peak_int[i][1]]
        x2, y2 = df2['Time'][mask2][peak_int[i][0]:peak_int[i][1]], df2[acsm_keys[i]][mask2][peak_int[i][0]:peak_int[i][1]]

        p1, = ax[i].plot(x1, y1, lw = 1, label = 'm/z 105', color = 'blue')
        p2, = ax[i].plot(x2, y2, lw = 1, label = 'm/z 183', color = 'green')

        ax2 = ax[i].twinx()
        
        x3, y3 = bc_dict[key]['Time'][peak_int[i][2]:peak_int[i][3]], bc_dict[key]['IR BCc'][peak_int[i][2]:peak_int[i][3]]
        p3, = ax2.plot(x3, y3, lw = 1, label = 'BC', color = 'k', zorder = -10)

        formatter = FuncFormatter(lambda s, x: time.strftime('%H:%M', time.gmtime(s)))
        ax[i].xaxis.set_major_formatter(formatter)
        ax[i].set_xticklabels(ax[i].get_xticklabels(), size = 8)
        ax2.xaxis.set_major_formatter(formatter)
        ax2.set_xticklabels(ax2.get_xticklabels(), size = 8)

        ylim = np.array(ax[i].get_ylim())
        ratio = ylim / np.sum(np.abs(ylim))
        scale = ratio / np.min(np.abs(ratio))
        scale = scale / n[i]
        ax2.set_ylim(np.max(np.abs(ax2.get_ylim())) * scale)

        ax[i].tick_params(axis = 'y', labelcolor = 'k', size = 8)
        ax2.tick_params(axis = 'y', labelcolor = p3.get_color(), size = 8)

        ax[i].legend(frameon = False, fontsize = 8, handles = [p1, p2, p3])

        ax[i].set_xlabel('Time', fontsize = 8)
        ax[i].set_ylabel('Org / $\mu$g/m$^{3}$', color = 'k', fontsize = 8)
        ax2.set_ylabel('BC / $\mu$g/m$^{3}$', color = p3.get_color(), fontsize = 8)

        ax[i].set_title(acsm_keys[i], fontsize = 9)