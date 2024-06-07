#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from iminuit import Minuit
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os, sys
#%%
def read_cal_spec(path, parent_path):
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

def linear_fit(x, y, a_guess, b_guess):

    Npoints = len(y)

    def fit_func(x, a, b):
        return b + (a * x)

    def least_squares(a, b) :
        y_fit = fit_func(x, a, b)
        squares = np.sum((y - y_fit)**2)
        return squares
    least_squares.errordef = 1.0    # Chi2 definition (for Minuit)

    # Here we let Minuit know, what to minimise, how, and with what starting parameters:   
    minuit = Minuit(least_squares, a = a_guess, b = b_guess)

    # Perform the actual fit:
    minuit.migrad();

    # Extract the fitting parameters:
    a_fit = minuit.values['a']
    b_fit = minuit.values['b']

    Nvar = 2                     # Number of variables 
    Ndof_fit = Npoints - Nvar    # Number of degrees of freedom = Number of data points - Number of variables

    # Get the minimal value obtained for the quantity to be minimised (here the Chi2)
    squares_fit = minuit.fval                          # The chi2 value

    # Calculate R2
    def simple_model(b):
        return b

    def least_squares_simple(b) :
        y_fit = simple_model(b)
        squares = np.sum((y - y_fit)**2)
        return squares
    least_squares_simple.errordef = 1.0    # Chi2 definition (for Minuit)

    # Here we let Minuit know, what to minimise, how, and with what starting parameters:   
    minuit_simple = Minuit(least_squares_simple, b = b_guess)

    # Perform the actual fit:
    minuit_simple.migrad();

    # Get the minimal value obtained for the quantity to be minimised (here the Chi2)
    squares_simple = minuit_simple.fval                          # The chi2 value

    R2 = 1 - (squares_fit / squares_simple)

    # Print the fitted parameters
    print(f"Fit: a={a_fit:6.6f}  b={b_fit:5.3f}  R2={R2:6.6f}")
    
    return a_fit, b_fit, squares_fit, Ndof_fit, R2

def plot_with_LinReg(ax, data_dict, df_keys, x_plot, a_guess, b_guess, lbl, clr, frame, ax_labels):
    a_array = np.zeros(len(a_guess))
    b_array = np.zeros(len(b_guess))
    df_fitted = pd.DataFrame({'MW': x_plot})

    for i, key in enumerate(data_dict.keys()):
        df = data_dict[key]
        x, y = df[df_keys[i][0]], df[df_keys[i][1]]

        lbl_fit = lbl[i] + ' fit'

        a, b, squares, ndof, R2 = linear_fit(x, y, a_guess[i], b_guess[i])
        a_array[i] = a
        b_array[i] = b
        y_fit = a*x_plot + b
        df_fitted[lbl[i]] = y_fit

        ax.plot(x_plot, y_fit, label = lbl_fit, color = clr[i], lw = 1.2)
        ax.scatter(x, y, label = lbl[i], color = clr[i], s = 10)

    ax.legend(fontsize = 8, frameon = frame)
    ax.tick_params(axis = 'both', which = 'major', direction = 'out', bottom = True, left = True, labelsize = 8)
    ax.set(xlabel = ax_labels[0], ylabel = ax_labels[1])

    return a_array, b_array, df_fitted