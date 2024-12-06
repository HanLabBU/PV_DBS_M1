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
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'Summary' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %%
# Loop through and plot all of the trials for both 40 Hz and 140 Hz
test_freqs = df['stim_freq'].unique()
cm = 0.394
zsco = 1

# Determine whether or not to display names
display_names = 0
test_freqs = df['stim_freq'].unique() # [40] #
for freq in test_freqs:
    freq_df = df[df['stim_freq'] == freq]
    vm_df = freq_df.pivot(index=['mouse_id', 'session_id', 'fov_id', 'trial_id'], columns='interp_time', values='interp_raw_detrend_trace')
    stim_df = freq_df.pivot(index=['mouse_id', 'session_id', 'fov_id', 'trial_id'], columns='interp_time', values='stim_raster')
    flicker_df = freq_df.pivot(index=['mouse_id', 'session_id', 'fov_id', 'trial_id'], columns='interp_time', values='flicker_raster')
    
    # Grab the number of trials for each multi-index
    num_trial_df = vm_df.groupby(level=['mouse_id', 'session_id', 'fov_id'], sort=False).size()

    # Remove fovs that had no trials
    filtered_groups = num_trial_df[num_trial_df > 0].index
    num_trial_df = num_trial_df.loc[filtered_groups]

    #DEBUG
    print(num_trial_df)

    # Store neuron names in a list
    pairings = num_trial_df.index.unique()
    nr_names = pairings.to_list()

    # Create neuron boundaries
    neuron_bound = np.array([0.5])
    for sz in num_trial_df.values:
        print(sz)
        neuron_bound = np.append(neuron_bound, neuron_bound[-1] + sz)

    #print(stim_df.iloc[0])

    #print('Summed mean')
    print(stim_df.mean(axis=0).sum())
    #print('Summed std')
    #print(stim_df.std(axis=0).sum())

    # Setup time variables
    avg_stim_raster = stim_df.mean(axis=0).to_numpy()
    avg_flicker_raster = flicker_df.mean(axis=0).to_numpy()
    timeline = vm_df.columns.values

    t = vm_df.columns.values[avg_stim_raster > 0]
    #plt.figure()
    #plt.plot(t)
    #plt.title(freq)
    block_idx = np.arange(vm_df.shape[0]) + 1 # Create row index array
    
    # Find the onset and offset of each burst period
    gap_start_idx = consts.find_start_idx(stim_df.iloc[0].values)
    gap_end_idx = avg_stim_raster.shape[0] - consts.find_start_idx(stim_df.iloc[0].values[::-1]) - 1
    gap_start_idx = gap_start_idx.astype(int)
    gap_end_idx = gap_end_idx.astype(int)[::-1]
    
    # Setup figure plotting
    plt.rcParams['font.size'] = 14

    fig, ax = plt.subplots(1, 2)
    fig.set_size_inches(40*cm, 35*cm)

    ax[0].set_position(Bbox.from_extents(1, 1, 40*cm, 35*cm))
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].spines['left'].set_visible(False)

    # Determine to zscore
    if zsco == 1:
        vm_df = vm_df.apply(stats.zscore, axis=1)

    # Plot the heatmap
    surf_plot = ax[0].pcolormesh(timeline, block_idx, vm_df.values)

    # Move legend depending on whether or not names are displayed
    if display_names == 1:
        loc = 'right'
    else:
        loc = 'left'

    cbar = plt.colorbar(surf_plot, location=loc, label='Vm zscore: ' + str(bool(zsco)))

    # Label the neuron names mid-way through
    if display_names == 1:
        new_yticks = neuron_bound[1:]
        new_yticks = new_yticks - num_trial_df.values/2
        ax[0].set_yticks(new_yticks)
        ax[0].set_yticklabels(nr_names, rotation=45)
    else:
        new_yticks = neuron_bound[1:]
        new_yticks = new_yticks - num_trial_df.values/2
        ax[0].set_yticks(new_yticks)
        ax[0].set_yticklabels(np.arange(1, new_yticks.shape[0]+1).astype(str))
        
    # Plot the vertical stimulation onset and offset
    t_idx = np.concatenate((gap_start_idx, gap_end_idx))
    t = vm_df.columns.values[t_idx]
    ax[0].vlines(x=t, ymin=1, ymax=np.max(block_idx), linestyles='--', colors=consts.stim_color)

    # Plot the stimulation time bars
    stim_y = np.max(block_idx) + 2
    stim_h = 2
    for start_idx, end_idx in zip(gap_start_idx, gap_end_idx):
        anch_y = stim_y - stim_h/2
        stim_w = vm_df.columns.values[end_idx] - vm_df.columns.values[start_idx]
        rect_patch = patches.Rectangle((vm_df.columns.values[start_idx], anch_y), stim_w, stim_h, facecolor=consts.stim_color)
        ax[0].add_patch(rect_patch)
    ax[0].hlines(y=stim_y, xmin=np.min(timeline), xmax=np.max(timeline), linestyles='-', colors=consts.stim_color)

    # Plot the flicker stimulation wave
    ax[0].plot(timeline, 2*avg_flicker_raster + 5 + vm_df.shape[0], '-b')
       
    # Plot all of the neuron boundaries
    ax[0].hlines(y=neuron_bound, xmin=np.min(timeline), xmax=np.max(timeline), colors='k')
    
    # Rename the Yticks
    #TODO calculate the midpoint of all neuronbounds
    #ax.set_yticks(ytick_loc)
    #ax.set_yticklabels(np.arange(neuron_bound.shape[0]) + 1) 
    
    ax[0].set_xlabel('Time (s)')
    ax[0].set_ylabel('Neuron # (Black outlines indicate all blocks per neuron)')
    ax[0].set_title('Vm')

    # Reduce the X limits to zoom-in and save all trials
    og_xlim = ax[0].get_xlim()
    ax[0].set_xlim(0.9, 2)
    fig.savefig(savefig_path + 'ZoomFlick_Summary_Blocks_Vm_' + str(freq) + '.png', format='png')
    fig.savefig(savefig_path + 'ZoomFlick_Summary_Blocks_Vm_' + str(freq) + '.pdf', format='pdf')
    ax[0].set_xlim(og_xlim)

    # ----------------- Plotting the spike raster -----------------
    spike_df = freq_df.pivot(index=['mouse_id', 'session_id', 'fov_id', 'trial_id'],\
                columns='interp_time', values='interp_spike_raster')

    # Do actual spike plotting using scatter
    #y, x = np.where(spike_df.values)
    #ax[1].scatter(y, x, marker='.')

    # Cleanup the plot
    ax[1].set_position(Bbox.from_extents(1, 1, 11*cm, 9*cm))
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].spines['left'].set_visible(False)


    # Pcolor plotting of spikes
    surf_plot = ax[1].pcolormesh(timeline, block_idx, spike_df.values, cmap='gray_r', vmin=0, vmax=1)

    ax[1].vlines(x=t, ymin=1, ymax=np.max(block_idx), linestyles='--', colors=consts.stim_color)

    # Plot the stimulation time bars
    stim_y = np.max(block_idx) + 2
    stim_h = 2
    for start_idx, end_idx in zip(gap_start_idx, gap_end_idx):
        anch_y = stim_y - stim_h/2
        stim_w = spike_df.columns.values[end_idx] - spike_df.columns.values[start_idx]
        rect_patch = patches.Rectangle((spike_df.columns.values[start_idx], anch_y), stim_w, stim_h, facecolor=consts.stim_color)
        ax[1].add_patch(rect_patch)
    ax[1].hlines(y=stim_y, xmin=np.min(timeline), xmax=np.max(timeline), linestyles='-', colors=consts.stim_color)


    # Plot the flicker stimulation wave
    ax[1].plot(timeline, 2*avg_flicker_raster + 5 + spike_df.shape[0], '-b')

    # Plot all of the neuron boundaries
    ax[1].hlines(y=neuron_bound, xmin=np.min(timeline), xmax=np.max(timeline), colors='k')
    ax[1].set_title('Spike Raster')

    fig.suptitle(str(freq) + 'Hz')    
    #fig.tight_layout()
    fig.subplots_adjust(left=0.2, right=0.85, top=0.85, bottom=0.2, wspace=.03)
    
    fig.savefig(savefig_path + 'Summary_Blocks_Vm_' + str(freq) + ' ' + str(display_names) + '.png', format='png')
    fig.savefig(savefig_path + 'Summary_Blocks_Vm_' + str(freq) + ' ' + str(display_names) + '.pdf', format='pdf')
    

plt.show()

# %%
plt.close('all')


