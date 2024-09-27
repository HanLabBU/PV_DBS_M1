# %%
# Import statements
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.transforms import Bbox
from matplotlib.figure import Figure
from matplotlib.colors import LinearSegmentedColormap
import sys
from glob import glob
from pymatreader import read_mat
from scipy import signal
from scipy import stats
from scipy.ndimage import uniform_filter1d

import itertools

#import matlab.engine
# Test out matlab spike detection code
#eng = matlab.engine.start_matlab()
#eng.eval("startup", nargout=0)

# Importing custom python files
import consts
import importlib
importlib.reload(consts)

from open_ephys.analysis import Session

%matplotlib tk

f = os.sep

# %%
# Read in pickle file
interm_data_path =  f + 'home' +f+ 'pierfier' +f+ 'Projects' +f+ 'Pierre Fabris' +f+ 'PV DBS neocortex' +f+ 'Interm_Data' +f+ 'flicker.pkl'
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'Compact' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %% Count number of Neurons in each stim frequency that have an increase in power
combo_df = df.groupby(['mouse_id', 'session_id', 'stim_freq', 'fov_id'])['Power_Delta']\
    .agg(lambda x: 1 if (x ==1).any() else 0).reset_index(name='Power Delta')

filt_combo_df = combo_df[combo_df['Power Delta'].notna()].reset_index()

# Loop through each stimulation frequency and count the neurons that had increase or decrease powers
for stim in filt_combo_df['stim_freq'].unique():
    filt_stim_df = filt_combo_df[filt_combo_df['stim_freq'] == stim]
    inc_pow_count = (filt_stim_df['Power Delta'] == 1).sum()
    dec_pow_count = (filt_stim_df['Power Delta'] == 0).sum()

    print('For ' + str(stim) + ' increase is ' + str(inc_pow_count))
    print('For ' + str(stim) + ' decrease is ' + str(dec_pow_count))


# %%
