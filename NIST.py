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
        with open(os.path.join(path, file)) as f:
            df = pd.read_table(f, sep = '\s+', skiprows=25, header=None)[:-1]
            df2 = pd.DataFrame(df.stack())
            df3 = pd.DataFrame(df2[0].str.split(',', expand = True).reset_index(drop=True)) #.sort_values(0)
            df3.columns = ['mass', ('intensity' + file)]

            for key in df3.keys():
                df3[key] = pd.to_numeric(df3[key])

            data_dict[file] = df3
    
    return data_dict
#%%
def merge_NIST(species_list, data):
    data_labels = []
    for specie in species_list:
        label = specie + '.txt'
        data_labels.append(label)
    
    merged = pd.DataFrame({'mass':[]})
    for label in data_labels:
        merged = pd.merge(merged, data[label], on = 'mass', how = 'outer')

    for key in merged.keys():
        merged[key] = merged[key].fillna(0)

    return merged
#%%
def plot_NIST(species_list, data, width, ax, xlim, ylim):
    merged = merge_NIST(species_list, data)

    bottom = np.zeros(len(merged['mass']))

    for i, key in enumerate(merged.keys()[1:]):
        y = merged[key] / 100
        x = merged['mass']

        ax.bar(x, y, width, label = species_list[i], bottom = bottom)
        bottom += y

        ax.legend(frameon = False)
        ax.set(xlabel = 'm/z', ylabel = 'Relative intensity', xlim = xlim, ylim = ylim)
#%%
def sum_columns_NIST(species_list, data):
    merged = merge_NIST(species_list, data)

    full_sum = np.zeros(len(species_list))
    sum_non_MolIon = np.zeros(len(species_list))
    relative_intensity = np.zeros(len(species_list))

    for i, key in enumerate(merged.keys()[1:]):
        full_sum[i] += np.sum(merged[key])
        sum_non_MolIon[i] += full_sum[i] - max(merged[key])
        relative_intensity[i] += sum_non_MolIon[i] / max(merged[key])
    
    return full_sum, sum_non_MolIon, relative_intensity
#%%
def sum_vs_molion(species_full, MWs, data):
    sums = pd.DataFrame(columns = ['Species', 'Molecular weight', 'Full sum', 'Non molecular ion', 'Relative intensity'])

    for i, specie in enumerate(species_full):
        MW = MWs[i]
        full, nonMol, Rel = sum_columns_NIST(specie, data)
        for i, f in enumerate(full):
            new_row = {'Species': specie[i], 'Molecular weight': MW, 'Full sum': full[i], 'Non molecular ion': nonMol[i], 'Relative intensity': Rel[i]}
            sums = pd.concat([sums, pd.DataFrame([new_row])], ignore_index=True)
    
    return sums
#%%
def Fracmentation_factor(species_list, data, molecular_weights):
    sums = pd.DataFrame(columns = ['Species', 'Molecular weight', 'Full sum', 'Molecular ion', 'Fragmentation factor'])

    for i, specie in enumerate(species_list):
        merged = merge_NIST(specie, data)

        full_sum = np.zeros(len(specie))
        sum_MolIon = np.zeros(len(specie))
        FF = np.zeros(len(specie))

        MW = molecular_weights[i]

        for j, key in enumerate(merged.keys()[1:]):
            mol_idx = merged[key].idxmax()
            Main_group = []

            for i in range(len(merged[key][:mol_idx])):
                rel_int = merged[key][mol_idx-i]
                if rel_int != 0:
                    Main_group.append(rel_int)
                if rel_int == 0:
                    break
            for i in range(len(merged[key])):
                if i > mol_idx:
                    Main_group.append(merged[key][i])

            full_sum[j] += np.sum(merged[key])
            sum_MolIon[j] += np.sum(Main_group)
            FF[j] += full_sum[j] / sum_MolIon[j]
        
        for j, f in enumerate(full_sum):
            new_row = {'Species': specie[j], 'Molecular weight': MW, 'Full sum': f, 'Molecular ion': sum_MolIon[j], 'Fragmentation factor': FF[j]}
            sums = pd.concat([sums, pd.DataFrame([new_row])], ignore_index=True)

    sums = sums.set_index('Species')

    return sums