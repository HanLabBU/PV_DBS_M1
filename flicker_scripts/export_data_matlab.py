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
from scipy.stats import zscore
from scipy.ndimage import uniform_filter1d
import fcwt
import seaborn as sns
import scipy.stats as stats
import scikit_posthocs as posthocs
import itertools
from scipy.io import savemat

%matplotlib inline

f = os.sep

# %% 
# Load dataframe
interm_data_path =  f + 'home' +f+ 'pierfier' +f+ 'Projects' +f+\
      'Pierre Fabris' +f+ 'PV DBS neocortex' +f+ 'Interm_Data' +f

interm_data_pathname =  interm_data_path + 'flicker.pkl'

savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab'\
      +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots'\
          +f+ 'Flicker' +f+ 'Spectra' +f
df = pd.read_pickle(interm_data_pathname)
print(df.columns)

# %%
# Save all trials, spike, flicker time, raster time, data to a matfile 
data_struct = {}

for stim_freq in df['stim_freq'].unique():
    stim_df = df[df['stim_freq'] == stim_freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    non_empty_fovs = [pair for pair in pairings\
                      if not stim_df[(stim_df['mouse_id'] == pair[0]) &\
                        (stim_df['session_id'] == pair[1]) &\
                        (stim_df['fov_id'] == pair[2])].empty ]
    
    # Initialize dictionary fields
    data_struct['f_' + str(stim_freq)] = {}
    all_rawvm_dict = {}
    all_neuron_name = {}
    all_spike_amp_raster_dict = {}
    flicker_raster = []

    # Save a single interp time for each frequency
    data_struct['f_' + str(stim_freq)]['interp_time'] = stim_df['interp_time'].unique()

    # Loop through each neuron
    roi_acum = 0
    for values in non_empty_fovs:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue

        # Save all of the trial raw_vms
        roi_raw_vm = fov_df.pivot(index='trial_id', columns='interp_time', values='detrend_trace').values
        all_rawvm_dict['n_' + str(roi_acum)] = np.transpose(roi_raw_vm)
        
        # Save all of the spike amplitudes
        roi_sp_amp = fov_df.pivot(index='trial_id', columns='interp_time', values='spike_amp_raster').values
        all_spike_amp_raster_dict['n_' + str(roi_acum)] = np.transpose(roi_sp_amp)
                    
        # Add neuron name
        all_neuron_name['n_' + str(roi_acum)] = "_".join([str(e) for e in values])
        
        flicker_raster = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['flicker_raster'].values

        roi_acum += 1

    # Convert arrays to dtype objects for better matlab conversion
    #data_struct['f_' + str(stim_freq)]['raw_vm'] = np.array(\
    #    all_raw_vm, dtype=object)
    
    # (Does not work) Add as a dictionary
    #data_struct['f_' + str(stim_freq)]['raw_vm'] = np.array(list(all_rawvm_dict.values()), dtype=object)

    data_struct['f_' + str(stim_freq)]['raw_vm'] = all_rawvm_dict
    data_struct['f_' + str(stim_freq)]['spike_amp_raster'] = all_spike_amp_raster_dict
    data_struct['f_' + str(stim_freq)]['nr_name'] = all_neuron_name
    data_struct['f_' + str(stim_freq)]['flicker_raster'] = flicker_raster
    
    
# DEBUGGING how to save things
#my_dict = {1: 'apple', 2: 'banana', 3:[1, 2,3]}
#my_dict = list(my_dict.values())

#sample = [np.array([[1, 2, 3]]), np.array([[1], [2], [3]])]
#sample = np.array(sample, dtype=object)
#data_struct['f_140']['sam'] = sample

savemat(interm_data_path + 'v1_flicker.mat', data_struct)