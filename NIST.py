#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
#%%
path = './NIST_data/'
files = os.listdir(path)
data_dict = {}

for file in files:
    with open(os.path.join(path, file))as f:
        df = pd.read_table(f, sep = '\s+', skiprows=25, header=None)[:-1]
        df2 = pd.DataFrame(df.stack())
        df3 = pd.DataFrame(df2[0].str.split(',', expand = True).reset_index(drop=True)) #.sort_values(0)
        df3.columns = ['mass', 'intensity']
        data_dict[file] = df3

print(data_dict)
# print(data_dict['Fluoranthene.jdx'])
#%%
print(len(data_dict['Pyrene.jdx']['mass']))
print(len(data_dict['Pyrene.jdx']['intensity']))
print(len(pd.to_numeric(data_dict['Pyrene.jdx']['intensity']) / 10))
#%%

bob = pd.merge(data_dict['Pyrene.jdx'], data_dict['Fluoranthene.jdx'], on='mass')
bob['intensity_x'] = pd.to_numeric(bob['intensity_x'])
bob['intensity_y'] = pd.to_numeric(bob['intensity_y'])
bob['pyrene+fluoranthracene'] = bob['intensity_x'] + bob['intensity_y']
bob['mass'] = pd.to_numeric(bob['mass'])
print(bob.head())
fig, ax = plt.subplots(figsize = (6.3, 3))

df = pd.DataFrame()
species = ['Pyrene', 'Fluoranthene']
width = 0.5
bottom = np.zeros(100)

#for specie in species:
#    data = data_dict[''.join([specie, '.jdx'])]
#    y = pd.to_numeric(data['intensity']) / 100
#    x = pd.to_numeric(data['mass'])
#    
#    ax.bar(x, y, width, label = specie) #, bottom = bottom)
    # bottom += y
ax.bar(bob['mass'], bob['intensity_x'] / 100, label='Pyrene')
ax.bar(bob['mass'], bob['intensity_y'] / 100, bottom=bob['intensity_x'] / 100, label='Fluoroanthracene')

ax.legend(frameon = False)
ax.set(xlabel = 'm/z', ylabel = 'Relative intensity', xlim = (190, 210)) #, yscale = 'log') #

fig.tight_layout()
fig.savefig('MW202.png', dpi = 200)
plt.show()
# %%
