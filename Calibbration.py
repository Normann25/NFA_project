#%%
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import os
#%%
def mass_conc(IEs, MWs, IE, I_s, CE_s, mIE, O):
    RIE_s = (IEs / MWs) / (IE / 62)

    C_PAH = np.sum(I_s) / (CE_s * RIE_s * mIE * O)

    return RIE_s, C_PAH