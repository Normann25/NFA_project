#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
#%%
def read_NIST(path):
    files = os.listdir(path)
    data_dict = {}

    for file in files:
        with open(os.path.join(path, file))as f:
            df = pd.read_table(f, sep = '\s+', skiprows=25, header=None)[:-1]
            df2 = pd.DataFrame(df.stack())
            df3 = pd.DataFrame(df2[0].str.split(',', expand = True).reset_index(drop=True)) #.sort_values(0)
            df3.columns = ['mass', 'intensity']
            data_dict[file] = df3
    
    return data_dict
#%%
def plot_NIST(species_list, data, width, ax, xlim, ylim):
    data_labels = []
    for specie in species_list:
        label = ''.join([specie, '.jdx'])
        data_labels.append(label)
    
    # Find a better, more effiecient way to merge the dataframes
    if len(data_labels) == 2:
        merged = pd.merge(data[data_labels[0]], data[data_labels[1]], on = 'mass')    
    if len(data_labels) == 3:
        merged = pd.merge(data[data_labels[0]], data[data_labels[1]], data[data_labels[2]], on = 'mass')
    if len(data_labels) == 4:
        merged = pd.merge(data[data_labels[0]], data[data_labels[1]], data[data_labels[2]], data[data_labels[3]], on = 'mass')
    if len(data_labels) == 5:
        merged = pd.merge(data[data_labels[0]], data[data_labels[1]], data[data_labels[2]], data[data_labels[3]], data[data_labels[4]], on = 'mass')
    
    for key in merged.keys():
        merged[key] = pd.to_numeric(merged[key])
    
    bottom = np.zeros(len(merged['mass']))

    for i, key in enumerate(merged.keys()[1:]):
        y = merged[key] / 100
        x = merged['mass']

        ax.bar(x, y, width, label = species_list[i], bottom = bottom)
        bottom += y

        ax.legend(frameon = False)
        ax.set(xlabel = 'm/z', ylabel = 'Relative intensity', xlim = xlim, ylim = ylim)
