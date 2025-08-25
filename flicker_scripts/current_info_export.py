# This file is meant to export the data info for each unique neuron

# %% Import statements
import os
import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import sys
from glob import glob
from pymatreader import read_mat
from scipy import signal
from scipy import stats
from scipy.ndimage import uniform_filter1d
import time
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
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Current' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %% Loop through each stim frequency
current_df = pd.DataFrame()
test_freqs = df['stim_freq'].unique() # TODO change to 
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue
        
        # Stop execution
        #raise Exception("Stopping execution")

        # Add the info parameter to the larger dataframe
        nr_dict = {
            'mouse':fov_df['mouse_id'].unique(),
            'session':fov_df['session_id'].unique(),
            'stim_freq':freq,
            'nr_params':fov_df['stim_param'].unique()
        }
        current_df = pd.concat([current_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')

# Write all of the values to a csv file
current_df.to_csv(savefig_path + 'flicker_current.csv', index=False)
