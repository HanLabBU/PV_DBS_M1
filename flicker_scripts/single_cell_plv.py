# %%
# Import statements
import os
import pandas as pd
import numpy as np
import random
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

import matlab.engine
# Test out matlab spike detection code
eng = matlab.engine.start_matlab()
eng.eval("startup", nargout=0)
eng.addpath('..')

# Importing custom python files
import consts
import importlib
importlib.reload(consts)

from open_ephys.analysis import Session

%matplotlib inline

f = os.sep

# Remove warnings from commandline
import warnings
warnings.filterwarnings("ignore")

# %%
# Read in pickle file
interm_data_path =  f + 'home' +f+ 'pierfier' +f+ 'Projects' +f+\
      'Pierre Fabris' +f+ 'PV DBS neocortex' +f+ 'Interm_Data' +f+ 'flicker.pkl'
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab'\
      +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots'\
          +f+ 'Flicker' +f+ 'PLV' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %%
# Calculate the flicker PLV compared to shuffled distribution

num_iter = 500 # TODO change back to 500
wind_dist = 1000 / 1000 # ms
cm = 0.394
samp_freq = 500 # In Hz

# Specify frequency range to probe for each frequency
freqs = np.arange(1, 50+1, 1) # np.array([100, 150]) # TODO change back

# Specify frequency range for filtering of flicker frequency


test_freqs = df['stim_freq'].unique() # [40] #    
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    # Create a list of ways to iterate through each type of neuron
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()
    pairings = list(itertools.product(mouse_list, session_list, fov_list))

    # Pre filter the pairings
    #non_empty_pairings = [
    #(mouse, ses, fov) for mouse, ses, fov in pairings\
    #                if not stim_df[stim_df['mouse_id'] == mouse,
    #                               stim_df['session_id'] == ses,
    #                               stim_df['fov_id'] == fov].empty
    #]

    non_empty_fovs = [pair for pair in pairings\
                      if not stim_df[(stim_df['mouse_id'] == pair[0]) &\
                        (stim_df['session_id'] == pair[1]) &\
                        (stim_df['fov_id'] == pair[2])].empty ]

    # Setup figures
    fig, ax = plt.subplots(len(non_empty_fovs), 3)
    fig.set_size_inches(40.308 * cm, 80.394 * cm)

    #TODO expand this to show the filter signal with the raster overlaid
    # As in like to the left of it so it is easier to see 'how' entrained it looks
    # Loop through unique FOVs
    roi_acum = 0
    for fov_i, values in enumerate(non_empty_fovs):
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue

        # PLVs across the main frequency range
        plvs = np.array([])
        plvs_adj = np.array([])

        # Loop through each trial
        for fr in freqs:
            # Frequency band limits
            low_freq = fr*0.95
            high_freq = fr*1.05
            tr_acum = 0
            for tr_val, trial_df in fov_df.groupby('trial_id'):
                filt_sig = consts.fir_filt(trial_df['interp_subvm'].values,\
                                      low_freq, high_freq, samp_freq, numtaps=30000)
                analytic_sig = signal.hilbert(filt_sig)
                inst_phase = np.angle(analytic_sig)

                #DEBUG trying to determine why the spectra at 100 looks different
                if fr == 100 or fr == 150:
                    fig2 = plt.figure()
                    ax2 = fig2.add_subplot(111)
                    ax2.plot(filt_sig)
                    ax2.set_title('Testing filtered ' + str(fr))
                    plt.show()

                # Save filterd and phase values of signal 
                #TODO switch and append the frequency label and just save all of these columns
                fov_df.loc[fov_df['trial_id'] == tr_val, 'filt_sig' + str(fr)] = filt_sig
                fov_df.loc[fov_df['trial_id'] == tr_val, 'inst_phase' + str(fr)] = inst_phase

                diff_arr = np.diff(trial_df['flicker_raster'].values)
                diff_arr = np.concatenate((np.array([0]), diff_arr), axis=0)
                peak_idxs, _ = signal.find_peaks(diff_arr)

                flick_onset_raster = np.zeros_like(trial_df['flicker_raster'].values)
                flick_onset_raster[peak_idxs] = 1

                fov_df.loc[fov_df['trial_id'] == tr_val, 'flick_onset_raster'] = flick_onset_raster
                tr_acum += 1

                # Plot the filter signal
                if fr == stim_freq:
                    ax[fov_i, 1].plot(trial_df['interp_time'], 0.02*tr_acum + filt_sig)
                    ax[fov_i, 1].plot(trial_df['interp_time'], 0.02*tr_acum + np.max(filt_sig) + 0.04*flick_onset_raster)

            #tr_dict = {
            #    'trial':tr_val,
            #    'timeline':trial_df['']
            #    'filt_sig':filt_sig
            #}        
            #phase_tr_df = pd.concat([phase_tr_df, pd.DataFrame(tr_dict)], ignore_index=True, join='outer')
                        
            # Calculate flicker-Vm PLV, only in the intial flicker onset period
            PLV, PLV2, norm_vecs = consts.event_PLV(\
                fov_df[(fov_df['interp_time'] >= 0) & (fov_df['interp_time'] < 1)],                                    
                'inst_phase' + str(fr), 'flick_onset_raster', 10)

            plvs = np.append(plvs, PLV)
            plvs_adj = np.append(plvs_adj, PLV2)

        # Calculate flicker-Vm PLV, only in the intial flicker onset period
        PLV, PLV2, norm_vecs = consts.event_PLV(\
            fov_df[(fov_df['interp_time'] >= 0) & (fov_df['interp_time'] < 1)],                                    
            'inst_phase8', 'flick_onset_raster', 10)

        #TODO double check that the 8Hz frequency is being probed here, check the new column name from above
        obs_plv = PLV
        obs_plv2 = PLV2

        # Initialize variables for storing the shuffled data
        shuf_plv = np.array([])
        shuf_plv_adj = np.array([])

        # Determine where is the last possible index for randomization

        #Perform the randomization
        temp_fov_df = fov_df.copy(deep=True)
        one_sec_idx_len = fov_df[(fov_df['trial_id'] == fov_df['trial_id'].unique()[0]) &\
                        (fov_df['interp_time'] >= 0) &\
                        (fov_df['interp_time'] < 1)].shape[0]
        
        for i in range(num_iter):
            tr_acum = 0
            for tr_val, trial_df in temp_fov_df.groupby('trial_id'):
                # Use the original and untainted flicker raster
                og_flick_raster = fov_df[fov_df['trial_id'] == tr_val]['flick_onset_raster'].values
                
                onset_idx = np.where(og_flick_raster == 1)[0][0]

                # Pick a time between the original onset and 0 for this current trial
                rand_start_idx = random.randint(0, onset_idx)
                new_flick_raster = np.roll(og_flick_raster, rand_start_idx - onset_idx)

                # Beyond one second from the onset needs to be set to zero so that only this one second PLV is calculated
                new_flick_raster[rand_start_idx + one_sec_idx_len:] = 0
                temp_fov_df.loc[temp_fov_df['trial_id'] == tr_val, 'flick_onset_raster'] = new_flick_raster

                #DEBUG
                #event_idx = np.where(new_flick_raster == 1)[0]
                #plt.plot(trial_df['interp_time'].values, 0.02*tr_acum + trial_df['filt_sig'].values)
                #plt.plot(trial_df['interp_time'].values, 0.02*tr_acum + 0.003*new_flick_raster)
                
                tr_acum += 1
            
            # Calculate the shuffled flicker PLV
            PLV, PLV2, norm_vecs = consts.event_PLV(temp_fov_df,                                    
                'inst_phase8', 'flick_onset_raster', 10)
            
            shuf_plv = np.append(shuf_plv, PLV)
            shuf_plv_adj = np.append(shuf_plv_adj, PLV2)

        # Construct a mask for this neuron in the original dataframe
        cur_nr_mask = (df['stim_freq'] == stim_freq) & \
                (df['mouse_id'] == values[0]) & \
                (df['session_id'] == values[1]) & \
                (df['fov_id'] == values[2])

        # Determine the upper percentile for significance
        perc95 = np.percentile(shuf_plv_adj, 95)
        if obs_plv2 > perc95:
            df.loc[cur_nr_mask, 'plv_etrain'] = 1
            color = 'green'
        else:
            df.loc[cur_nr_mask, 'plv_etrain'] = 0
            color = 'red'

        # Plot the distribution
        ax[fov_i, 0].hist(shuf_plv_adj, bins=1000, alpha=0.7, color='black')
        ax[fov_i, 0].vlines(x=[obs_plv2], ymin=0.5, ymax=10, linestyles='-', colors=color)
        
        # Plot the PLVs across the needed frequency range
        ax[fov_i, 2].plot(freqs, plvs, label='PLVs')

        fig.tight_layout()
        plt.show()

    #raise Exception('Stopping Execution')    
df['plv_etrain'] = df['plv_etrain'].astype('category')

# %% Close all of the windows
plt.close('all')
eng.eval("close all", nargout=0)

# %% Save the updated dataframe to pickle file
df.to_pickle(interm_data_path)

# %% Test PLVs of DBS
num_iter = 500 # TODO change back to 500
wind_dist = 1000 / 1000 # ms
cm = 0.394
samp_freq = 500 # In Hz

# Specify frequency range to probe for each frequency
freqs = np.arange(1, 150+1, 1)

test_freqs = df['stim_freq'].unique() # [40] #      
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    # Create a list of ways to iterate through each type of neuron
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()
    pairings = list(itertools.product(mouse_list, session_list, fov_list))

    # Pre filter the pairings
    #non_empty_pairings = [
    #(mouse, ses, fov) for mouse, ses, fov in pairings\
    #                if not stim_df[stim_df['mouse_id'] == mouse,
    #                               stim_df['session_id'] == ses,
    #                               stim_df['fov_id'] == fov].empty
    #]

    non_empty_fovs = [pair for pair in pairings\
                      if not stim_df[(stim_df['mouse_id'] == pair[0]) &\
                        (stim_df['session_id'] == pair[1]) &\
                        (stim_df['fov_id'] == pair[2])].empty ]

    # Setup figures
    fig, ax = plt.subplots(len(non_empty_fovs), 3, constrained_layout=True)
    fig.set_size_inches(40.308 * cm, 80.394 * cm)

    #TODO expand this to show the filter signal with the raster overlaid
    # As in like to the left of it so it is easier to see 'how' entrained it looks
    # Loop through unique FOVs
    roi_acum = 0
    for fov_i, values in enumerate(non_empty_fovs):
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue

        # PLVs across the main frequency range
        plvs = np.array([])
        plvs_adj = np.array([])

        # Loop through each trial
        for fr in freqs:
            # Frequency band limits
            low_freq = fr*0.05
            high_freq = fr*0.95
            tr_acum = 0
            for tr_val, trial_df in fov_df.groupby('trial_id'):
                filt_sig = consts.fir_filt(trial_df['interp_subvm'].values,\
                                      low_freq, high_freq, samp_freq, numtaps=30000)
                analytic_sig = signal.hilbert(filt_sig)
                inst_phase = np.angle(analytic_sig)

                # Save filterd and phase values of signal 
                #TODO switch and append the frequency label and just save all of these columns
                fov_df.loc[fov_df['trial_id'] == tr_val, 'filt_sig' + str(fr)] = filt_sig
                fov_df.loc[fov_df['trial_id'] == tr_val, 'inst_phase' + str(fr)] = inst_phase

                diff_arr = np.diff(trial_df['stim_raster'].values)
                diff_arr = np.concatenate((np.array([0]), diff_arr), axis=0)
                peak_idxs, _ = signal.find_peaks(diff_arr)

                stim_onset_raster = np.zeros_like(trial_df['stim_raster'].values)
                stim_onset_raster[peak_idxs] = 1

                fov_df.loc[fov_df['trial_id'] == tr_val, 'stim_onset_raster'] = stim_onset_raster
                tr_acum += 1

                # Plot the filter signal
                if fr == stim_freq:
                    ax[fov_i, 1].plot(trial_df['interp_time'], 0.02*tr_acum + filt_sig)
                    ax[fov_i, 1].plot(trial_df['interp_time'], 0.02*tr_acum + np.max(filt_sig) + 0.04*stim_onset_raster)

            #tr_dict = {
            #    'trial':tr_val,
            #    'timeline':trial_df['']
            #    'filt_sig':filt_sig
            #}        
            #phase_tr_df = pd.concat([phase_tr_df, pd.DataFrame(tr_dict)], ignore_index=True, join='outer')
                        
            # Calculate dbs-Vm PLV, only in the intial stim onset period
            PLV, PLV2, norm_vecs = consts.event_PLV(\
                fov_df[(fov_df['interp_time'] >= 1) & (fov_df['interp_time'] < 2)],                                    
                'inst_phase' + str(fr), 'stim_onset_raster', 10)

            plvs = np.append(plvs, PLV)
            plvs_adj = np.append(plvs_adj, PLV2)

        # Calculate dbs-Vm PLV, only in the intial flicker onset period
        PLV, PLV2, norm_vecs = consts.event_PLV(\
            fov_df[(fov_df['interp_time'] >= 1) & (fov_df['interp_time'] < 2)],                                    
            'inst_phase' + str(stim_freq), 'stim_onset_raster', 10)

        #TODO double check that the 8Hz frequency is being probed here, check the new column name from above
        obs_plv = PLV
        obs_plv2 = PLV2

        # Set the title with the PLV value
        ax[fov_i, 0].set_title('PLV is ' + str(obs_plv2))

        ax[fov_i, 2].plot(freqs, plvs_adj)

    plt.suptitle(stim_freq)
    