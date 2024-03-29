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
                df3[key] = pd.to_numeric(df3[key], errors = 'coerce')

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
        if 'intensityChrysene.txt' in key:
            y = merged[key] / 10
            x = merged['mass']
        elif 'intensityCoronene.txt' in key:
            y = merged[key] / 10
            x = merged['mass']
        else:
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
    sums = pd.DataFrame(columns = ['Species', 'Molecular weight', 'Full sum', 'MI v1', 'MI v2', 'MI v3', 'FF v1', 'FF v2', 'FF v3', 'MI fraction'])

    for i, specie in enumerate(species_list):
        merged = merge_NIST(specie, data)

        full_sum = np.zeros(len(specie))
        MI_frac = np.zeros(len(specie))
        sum_MI_v1 = np.zeros(len(specie))
        sum_MI_v2 = np.zeros(len(specie))
        sum_MI_v3 = np.zeros(len(specie))
        FF_v1 = np.zeros(len(specie))
        FF_v2 = np.zeros(len(specie))
        FF_v3 = np.zeros(len(specie))

        MW = molecular_weights[i]

        for j, key in enumerate(merged.keys()[1:]):
            mol_idx = merged[key].idxmax()

            Main_group = []

            for i in range(len(merged[key][:mol_idx-1])):
                rel_int = merged[key][mol_idx-i]
                if rel_int != 0:
                    Main_group.append(rel_int)
                if rel_int == 0:
                    break
            for i in range(len(merged[key])):
                if i > mol_idx:
                    Main_group.append(merged[key][i])

            main = merged[key][mol_idx] + merged[key][mol_idx+1] + merged[key][mol_idx-1] + merged[key][mol_idx+2] + merged[key][mol_idx-2]

            main2 = []
            for i in range(5):
                main2.append(merged[key][mol_idx-i])
            main2.append(merged[key][mol_idx+1])

            full_sum[j] += np.sum(merged[key])
            MI_frac[j] += merged[key][mol_idx] / full_sum[j]

            sum_MI_v1[j] += np.sum(Main_group)
            FF_v1[j] += full_sum[j] / sum_MI_v1[j]
            
            sum_MI_v2[j] += main
            FF_v2[j] += full_sum[j] / sum_MI_v2[j]

            sum_MI_v3[j] += np.sum(main2)
            FF_v3[j] += full_sum[j] / sum_MI_v3[j]
        
        for j, f in enumerate(full_sum):
            new_row = {'Species': specie[j], 'Molecular weight': MW, 'Full sum': f, 'MI v1': sum_MI_v1[j], 'MI v2': sum_MI_v2[j], 'MI v3': sum_MI_v3[j], 'FF v1': FF_v1[j], 'FF v2': FF_v2[j], 'FF v3': FF_v3[j], 'MI fraction': MI_frac[j]}
            sums = pd.concat([sums, pd.DataFrame([new_row])], ignore_index=True)

    sums = sums.set_index('Species')

    return sums