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

%matplotlib tk

f = os.sep

# %%
# Read in pickle file
interm_data_path =  f + 'home' +f+ 'pierfier' +f+ 'Projects' +f+\
      'Pierre Fabris' +f+ 'PV DBS neocortex' +f+ 'Interm_Data' +f+ 'flicker.pkl'
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab'\
      +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots'\
          +f+ 'Flicker' +f+ 'Pulse_Trig' +f
df =pd.read_pickle(interm_data_path)
print(df.columns)

# %%
# Perform neuronwise visual flicker triggered
extra_trace = 3

test_freqs = df['stim_freq'].unique()
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))

    # Initialize arrays
    all_flick_pre_stim_pulses = np.array([])
    all_flick_dur_stim_pulses = np.array([])
    all_flick_post_stim_pulses = np.array([])

    for values in pairings:#[pairings[0]]: # TODO change back to -> pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue
        
        # Loop through each trial
        for tr_val in fov_df['trial_id'].unique():
            # Single trial dataframe
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            
            # Get the average spike amplitude for normalization
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)

            # If there are no spikes detected, then just set the value to 1
            if np.isnan(avg_sp_amp):
                avg_sp_amp = 1

            # Find the frame index range before, during, and after ESTIM
            start_flicker_idx = np.where(trial_df['flicker_raster'] == 1)[0][0]
            end_flicker_idx = np.where(trial_df['flicker_raster'] == 1)[0][-1]
            start_stim_idx = np.where(trial_df['stim_raster'] == 1)[0][0]
            end_stim_idx = np.where(trial_df['stim_raster'] == 1)[0][-1]
            
            # Grab the initial rise corner of the flicker raster
            diff_arr = np.diff(trial_df['flicker_raster'].values)
            diff_arr = np.concatenate((np.array([0]), diff_arr), axis=0)
            peak_idxs, _ = signal.find_peaks(diff_arr)

            # Determine the windows idx width to store consistent pulse widths
            pulse_width = int(np.ceil(np.mean(np.diff(peak_idxs))))

            # Loop through each pulse start
            for i, pulse_idx in enumerate(peak_idxs):
                            
                cur_pulse = trial_df['interp_subvm'].values[pulse_idx - extra_trace:pulse_idx + pulse_width + extra_trace]/avg_sp_amp
                cur_pulse = cur_pulse.reshape(-1, 1)

                # Onset subtract the pulse
                cur_pulse = cur_pulse - cur_pulse[extra_trace]

                # Check which part of the trial the pulse lies
                if pulse_idx >= start_flicker_idx and pulse_idx < start_stim_idx:
                    try:
                        all_flick_pre_stim_pulses = np.concatenate(\
                            (all_flick_pre_stim_pulses, cur_pulse), axis=1)                                            

                    except Exception as e:
                        all_flick_pre_stim_pulses = cur_pulse
                
                # Check which part of the trial the pulse lies
                if pulse_idx >= start_stim_idx and pulse_idx < end_stim_idx:
                    try:
                        all_flick_dur_stim_pulses = np.concatenate(\
                            (all_flick_dur_stim_pulses, cur_pulse), axis=1)
                    
                    except Exception as e:
                        all_flick_dur_stim_pulses = cur_pulse
                                    
                # Check which part of the trial the pulse lies
                if pulse_idx >= end_stim_idx and pulse_idx < end_flicker_idx:
                    try:
                        all_flick_post_stim_pulses = np.concatenate(\
                            (all_flick_post_stim_pulses, cur_pulse), axis=1)
                    
                    except Exception as e:
                        all_flick_post_stim_pulses = cur_pulse

    # Calculate all of the plotting averages and SEMs
    num_pre_pulses = all_flick_pre_stim_pulses.shape[0]
    num_dur_pulses = all_flick_dur_stim_pulses.shape[0]
    num_post_pulses = all_flick_post_stim_pulses.shape[0]
    avg_pre_flick = np.mean(all_flick_pre_stim_pulses, axis=1)
    avg_dur_flick = np.mean(all_flick_dur_stim_pulses, axis=1)
    avg_post_flick = np.mean(all_flick_post_stim_pulses, axis=1)
    std_pre_flick = np.std(all_flick_pre_stim_pulses, axis=1)
    std_dur_flick = np.std(all_flick_dur_stim_pulses, axis=1)
    std_post_flick = np.std(all_flick_post_stim_pulses, axis=1)
    
    sem_pre_flick = std_pre_flick/np.sqrt(num_pre_pulses)
    sem_dur_flick = std_dur_flick/np.sqrt(num_dur_pulses)
    sem_post_flick = std_post_flick/np.sqrt(num_post_pulses)
    # Plot the average trace pulses
    fig, axs = plt.subplots(3, figsize=(4, 6))
    
    # Calculate the timeline
    timeline = trial_df['interp_time'].values[0:cur_pulse.shape[0]]
    timeline = timeline - timeline[0]
    timeline = timeline - timeline[extra_trace]
    timeline = timeline*1000

    # Parameter to adjust flicker raster height
    flick_height = 0.5

    # Grab a single pulse flicker raster
    flick_pulse_single = trial_df['flicker_raster'].values[peak_idxs[0] - extra_trace:peak_idxs[0] + pulse_width + extra_trace]
    axs[0].plot(timeline, flick_height*flick_pulse_single + np.max(avg_pre_flick))
    axs[0].fill(np.concatenate((timeline, timeline[::-1])), \
        np.concatenate((avg_pre_flick + sem_pre_flick, avg_pre_flick[::-1] - sem_pre_flick[::-1])),
        color='gray')
    axs[0].plot(timeline, avg_pre_flick)
    axs[0].set_title('Pre Stim')
        
    axs[1].plot(timeline, flick_height*flick_pulse_single + np.max(avg_dur_flick))
    axs[1].fill(np.concatenate((timeline, timeline[::-1])), \
        np.concatenate((avg_dur_flick + sem_dur_flick, avg_dur_flick[::-1] - sem_dur_flick[::-1])),
        color='gray')
    axs[1].plot(timeline, avg_dur_flick)
    axs[1].set_title('During Stim')
    axs[2].plot(timeline, flick_height*flick_pulse_single + np.max(avg_post_flick))
    axs[2].fill(np.concatenate((timeline, timeline[::-1])), \
        np.concatenate((avg_post_flick + sem_post_flick, avg_post_flick[::-1] - sem_post_flick[::-1])),
        color='gray')
    axs[2].plot(timeline, avg_post_flick)
    axs[2].set_title('Post Stim')
    plt.suptitle('Flicker Pulse Triggered Avg' + str(freq))
    plt.tight_layout()
    plt.show()
    save_filename = savefig_path + 'Pulse_Triggered_Flicker_' +\
        str(stim_df['stim_freq'].unique()[0]) + 'Hz'
    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')

