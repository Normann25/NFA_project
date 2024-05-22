#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os, sys
#%%
def read_data(path, parent_path):
    new_df = pd.DataFrame({'mass':[]})

    parentPath = os.path.abspath(parent_path)
    if parentPath not in sys.path:
        sys.path.insert(0, parentPath)

    files = os.listdir(path)

    for file in files:
        name = file.split('.')[0]
        with open(os.path.join(path, file)) as f:
            df = pd.read_table(f, sep = '\s+', skiprows = 1, header=None)
            df.columns = ['mass', name]
            new_df = pd.merge(new_df, df, on = 'mass', how = 'outer')

    for key in new_df.keys():
        new_df[key] = new_df[key].fillna(0)

    mask = new_df['mass'] >= 60
    new_df = new_df[mask]

    for key in new_df.keys()[1:]:
        mol_idx = new_df[key].idxmax()
        new_df[key] = (new_df[key] / new_df[key][mol_idx]) * 100
    
    return new_df

def Fracmentation_factor(species_list, data, molecular_weights):
    sums = pd.DataFrame(columns = ['Species', 'Molecular weight', 'Full sum', 'MI', 'FF', 'MI fraction'])

    full_sum = np.zeros(len(species_list))
    MI_frac = np.zeros(len(species_list))
    sum_MI = np.zeros(len(species_list))
    FF = np.zeros(len(species_list))

    for j, key in enumerate(species_list):
        mol_idx = data[key].idxmax()

        main = []
        for i in range(5):
            main.append(data[key][mol_idx-i])
        main.append(data[key][mol_idx+1])

        full_sum[j] += np.sum(data[key])
        MI_frac[j] += data[key][mol_idx] / full_sum[j]

        sum_MI[j] += np.sum(main)
        FF[j] += full_sum[j] / sum_MI[j]
    
        new_row = {'Species': key, 'Molecular weight': molecular_weights[j], 'Full sum': full_sum[j], 'MI': sum_MI[j], 'FF': FF[j], 'MI fraction': MI_frac[j]}
        sums = pd.concat([sums, pd.DataFrame([new_row])], ignore_index=True)

    sums = sums.set_index('Species')

    return sums

def plot_MS(ax, df, key, width):
    ax.bar(df['mass'], df[key], width)
    
    ax.set(xlabel = 'm/z', ylabel = 'Intensity')

    ax.xaxis.set_minor_locator(AutoMinorLocator())
    ax.yaxis.set_minor_locator(AutoMinorLocator())
    ax.tick_params(axis = 'both', which = 'major', direction = 'out', bottom = True, left = True, labelsize = 8)
    ax.tick_params(axis = 'both', which = 'minor', direction = 'out', width = 1, length = 2, bottom = True, left = True)
    ax.yaxis.offsetText.set_fontsize(9)