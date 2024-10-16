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

%matplotlib inline

f = os.sep

# %%
# Read in pickle file
interm_data_path =  f + 'home' +f+ 'pierfier' +f+ 'Projects' +f+ 'Pierre Fabris' +f+ 'PV DBS neocortex' +f+ 'Interm_Data' +f+ 'flicker.pkl'
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'Exemplary' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# Specify default values
cm = 0.394
samp_freq = 500

# %%
# 140 Hz example
exam_df = df[(df['mouse_id'] == '109567_Vb_male') & \
             (df['session_id'] == '20240411') & \
                (df['fov_id'] == 3)]

exam_df = exam_df[exam_df['stim_freq'] == 140]
trial_df = exam_df[exam_df['trial_id'] == exam_df['trial_id'].unique()[2]]

# %%
# 40 Hz example
exam_df = df[(df['mouse_id'] == '109567_Vb_male') & \
             (df['session_id'] == '20240311') & \
                (df['fov_id'] == 2)]

exam_df = exam_df[exam_df['stim_freq'] == 40]
trial_df = exam_df[exam_df['trial_id'] == exam_df['trial_id'].unique()[3]]

# %%
#an example trial from 109567_Vb_male_20240424_fov4_freq140
exam_df = df[(df['mouse_id'] == '109567_Vb_male') & \
             (df['session_id'] == '20240424') & \
                (df['fov_id'] == 4)]

exam_df = exam_df[exam_df['stim_freq'] == 140]
trial_df = exam_df[exam_df['trial_id'] == exam_df['trial_id'].unique()[5]]

# %%
# Another 140 Hz example
exam_df = df[(df['mouse_id'] == '109558_Vb_male') & \
             (df['session_id'] == '20240308') & \
                (df['fov_id'] == 3)]

exam_df = exam_df[exam_df['stim_freq'] == 140]
trial_i = 3
trial_df = exam_df[exam_df['trial_id'] == exam_df['trial_id'].unique()[trial_i]] # Unique 1, 3, 5, 6 are pretty good

# %%
# Plot individual cell trials
low_freq = 8*0.05
high_freq = 8*0.95
#low_freq = 7
#high_freq = 9

fig, axs = plt.subplots(3, sharex=True)
fig.set_size_inches(8.308 * cm, 8.394 * cm)
timeline = trial_df['interp_time'].values
stim_raster = trial_df['stim_raster'].values
spike_raster = trial_df['spike_raster'].values
flicker_start = trial_df[trial_df['flicker_raster'] == 1]['interp_time'].values[0]
timeline = timeline - flicker_start
trace_sbr = trial_df['detrend_trace'].values/trial_df['trace_noise'].values
axs[0].plot(timeline, trace_sbr, 'k', linewidth=0.3)

# Plot the 8 Hz flicker signal
filt_sig8 = consts.fir_filt(trace_sbr,
                            low_freq, high_freq, samp_freq, numtaps=30000)
axs[0].plot(timeline, filt_sig8, 'b', linewidth=0.2)

# Plot the spikes detected in the trial
spike_idx = np.where(spike_raster == 1)[0]
axs[0].plot(timeline[spike_idx], (np.max(trace_sbr) + 1)*np.ones_like(spike_idx), '.k', markersize=2)

# Plotting the stimulation indicators
axs[0].plot(timeline, trial_df['flicker_raster'] + 3 + np.max(trace_sbr), '-b', linewidth=0.3)
axs[0].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + np.max(trace_sbr)), color=consts.pulse_color, linewidth=0.3)
stim_idx = np.where(stim_raster == 1)[0]
axs[0].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + np.max(trace_sbr)), '|',\
         linewidth=20,  color=consts.pulse_color)

# Plot the line bars
#t_bar_length = 500# time in ms
#axs[0].plot([timeline[0], timeline[0] + t_bar_length/1000], [np.min(trace_sbr) - 0.2, np.min(trace_sbr) - 0.2], '-k')
#axs[0].text(timeline[0], np.min(trace_sbr) - 1.1, str(t_bar_length) + ' ms', color='k')

sbr_length = 5
axs[0].plot([timeline[0] + 0.2, timeline[0] + 0.2], [0, sbr_length], '-k')
axs[0].text(timeline[0] + 0.1, 0, str(sbr_length) + ' SBR', rotation=90)

axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].set_yticks([])

# Plot the trials heatpmap
exam_df['trace_sbr'] = exam_df['detrend_trace'] / exam_df['trace_noise']
plot_df = exam_df.pivot(index='trial_id', columns='interp_time', values='trace_sbr')
axs[1].spines['top'].set_visible(False)
axs[1].spines['right'].set_visible(False)
surf_p = axs[1].pcolormesh(timeline, 1 + np.arange(0, plot_df.index.nunique()), plot_df.values)
axs[1].set_ylabel('Trail #')

# Plot the trial average trace with error bar
avg_vm = plot_df.mean(axis=0)
num_trials = plot_df.shape[0]
std_vm = plot_df.std(axis=0)
sem_vm = std_vm /np.sqrt(num_trials)
axs[2].fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
axs[2].plot(timeline, avg_vm, '-k', linewidth=0.3)
axs[2].spines['top'].set_visible(False)
axs[2].spines['right'].set_visible(False)
axs[2].spines['left'].set_visible(False)
axs[2].set_yticks([])
axs[2].set_xlabel('Time from Flicker Stim (S)')

# Plot the SBR bar
sbr_length = 5
axs[2].plot([timeline[0] + 0.2, timeline[0] + 0.2], [0, sbr_length], '-k')
axs[2].text(timeline[0] + 0.1, 0, str(sbr_length) + ' SBR', rotation=90)
plt.xlim(-1, 4)

# Add the colorbar to the side of the heatmap
cbar_ax = fig.add_axes((0.92, 0.35, 0.02, 0.35))
cbar = fig.colorbar(surf_p, label='SBR', cax=cbar_ax)
save_filename = savefig_path + exam_df['mouse_id'].unique()[0] +\
            '_' + exam_df['session_id'].unique()[0] +\
            '_fov' + str(exam_df['fov_id'].unique()[0]) +\
            '_' + str(exam_df['stim_freq'].unique()[0]) +\
            'Hz_trial' + str(trial_i)
plt.savefig(save_filename + '.svg', format='svg')
plt.savefig(save_filename + '.png', format='png')

plt.show()

#TODO how to make the tick marks much larger????

# %%
# Plot the zoomins on single traces

# Specify  the zoom-in time ranges
#flicker_trange = np.array([0.2, 0.7])
#stim_trange = np.array([1.4, 1.9])
flicker_trange = np.array([0, 1])
stim_trange = np.array([1, 2])

# Setup the figure size
fig, axs = plt.subplots(1, 2)
fig.set_size_inches(15 * cm, 2 * cm)

# Grab all of the stimulation time stuff
timeline = trial_df['interp_time'].values
stim_raster = trial_df['stim_raster'].values
spike_raster = trial_df['spike_raster'].values
flicker_start = trial_df[trial_df['flicker_raster'] == 1]['interp_time'].values[0]
timeline = timeline - flicker_start

#-- Flicker period zoom in ---

# Plot the trace with SBR units
trace_sbr = trial_df['detrend_trace'].values/trial_df['trace_noise'].values
axs[0].plot(timeline, trace_sbr, 'k', linewidth=0.3)

filt_sig8 = consts.fir_filt(trace_sbr,
                            low_freq, high_freq, samp_freq, numtaps=30000)
axs[0].plot(timeline, filt_sig8, 'b', linewidth=0.2)

# Plot the spikes detected in the trial
spike_idx = np.where(spike_raster == 1)[0]
axs[0].plot(timeline[spike_idx], (np.max(trace_sbr) + 1)*np.ones_like(spike_idx), '.k', markersize=2)

# Plotting the stimulation indicators
axs[0].plot(timeline, trial_df['flicker_raster'] + 3 + np.max(trace_sbr), '-b', linewidth=0.3)
axs[0].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + np.max(trace_sbr)), color=consts.pulse_color, linewidth=0.3)
stim_idx = np.where(stim_raster == 1)[0]
axs[0].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + np.max(trace_sbr)), '|',\
         linewidth=20,  color=consts.pulse_color)

# Plot the line bars
#t_bar_length = 500# time in ms
#axs[0].plot([timeline[0], timeline[0] + t_bar_length/1000], [np.min(trace_sbr) - 0.2, np.min(trace_sbr) - 0.2], '-k')
#axs[0].text(timeline[0], np.min(trace_sbr) - 1.1, str(t_bar_length) + ' ms', color='k')

sbr_length = 5

# Show the SBR
axs[0].plot([flicker_trange[0] + 0.02, flicker_trange[0] + 0.02], [2, 2 + sbr_length], '-k')
axs[0].text(flicker_trange[0] - 0.01, 2.3, str(sbr_length) + ' SBR', rotation=90)

# Clean up the plot
axs[0].spines['top'].set_visible(False)
axs[0].spines['left'].set_visible(False)
axs[0].spines['right'].set_visible(False)
axs[0].spines['bottom'].set_visible(False)
axs[0].set_yticks([])
axs[0].set_xticks([])

# Find the zoom in
axs[0].set_xlim(flicker_trange[0], flicker_trange[1])

# Show timescale bar
time_length = 50 # in ms
axs[0].plot(.1 + np.array([flicker_trange[0], flicker_trange[0] + time_length/1000]), [-3.5, -3.5], '-k')
axs[0].text(.07 + flicker_trange[0], -4.7, str(time_length) + ' ms')



# -- Plot the estim period zoom in --

# Plot the trace with SBR units
trace_sbr = trial_df['detrend_trace'].values/trial_df['trace_noise'].values
axs[1].plot(timeline, trace_sbr, 'k', linewidth=0.3)
axs[1].plot(timeline, filt_sig8, 'b', linewidth=0.2)

# Plot the spikes detected in the trial
spike_idx = np.where(spike_raster == 1)[0]
axs[1].plot(timeline[spike_idx], (np.max(trace_sbr) + 1)*np.ones_like(spike_idx), '.k', markersize=2)

# Plotting the stimulation indicators
axs[1].plot(timeline, trial_df['flicker_raster'] + 3 + np.max(trace_sbr), '-b', linewidth=0.3)
axs[1].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + np.max(trace_sbr)), color=consts.pulse_color, linewidth=0.3)
stim_idx = np.where(stim_raster == 1)[0]
axs[1].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + np.max(trace_sbr)), '|',\
         linewidth=20,  color=consts.pulse_color)

# Plot the line bars
#t_bar_length = 500# time in ms
#axs[0].plot([timeline[0], timeline[0] + t_bar_length/1000], [np.min(trace_sbr) - 0.2, np.min(trace_sbr) - 0.2], '-k')
#axs[0].text(timeline[0], np.min(trace_sbr) - 1.1, str(t_bar_length) + ' ms', color='k')

sbr_length = 5

# Show the SBR
axs[1].plot([stim_trange[0] + 0.02, stim_trange[0] + 0.02], [2, 2 + sbr_length], '-k')
axs[1].text(stim_trange[0] - 0.01, 2.3, str(sbr_length) + ' SBR', rotation=90)

# Clean up the plot
axs[1].spines['top'].set_visible(False)
axs[1].spines['left'].set_visible(False)
axs[1].spines['right'].set_visible(False)
axs[1].spines['bottom'].set_visible(False)
axs[1].set_yticks([])
axs[1].set_xticks([])

# Find the zoom in
axs[1].set_xlim(stim_trange[0], stim_trange[1])

# Show timescale bar
time_length = 50 # in ms
axs[1].plot(.1 + np.array([stim_trange[0], stim_trange[0] + time_length/1000]), [-3.5, -3.5], '-k')
axs[1].text(.07 + stim_trange[0], -4.7, str(time_length) + ' ms')

# Save the zoom in of exemplary trace
save_filename = savefig_path + exam_df['mouse_id'].unique()[0] +\
            '_' + exam_df['session_id'].unique()[0] +\
            '_fov' + str(exam_df['fov_id'].unique()[0]) +\
            '_' + str(exam_df['stim_freq'].unique()[0]) +\
            'Hz_zoomins_trial' + str(trial_i)

plt.savefig(save_filename + '.svg', format='svg')
plt.savefig(save_filename + '.png', format='png')

# %%
# Plot and save all of the heatmaps and average plots individually
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'All_Exemplary' +f

mouse_list = df['mouse_id'].unique()
session_list = df['session_id'].unique()
fov_list = df['fov_id'].unique()
roi_list = df['roi_id'].unique()
stim_list = df['stim_freq'].unique()

pairings = list(itertools.product(stim_list, mouse_list, session_list, fov_list, roi_list))

# Loop through all unqiue pairings
for values in pairings:
    exam_df = df[(df['stim_freq'] == values[0]) & \
                 (df['mouse_id'] == values[1]) & \
                 (df['session_id'] == values[2]) & \
                 (df['fov_id'] == values[3]) & \
                 (df['roi_id'] == values[4])]

    # Skip if the dataframe is empty
    if exam_df.empty:
        continue
    
    trial_df = exam_df[exam_df['trial_id'] == exam_df['trial_id'].unique()[0]]
    
    # Setup figure
    cm = 0.394
    fig, axs = plt.subplots(3, sharex=True)
    fig.set_size_inches(50 * cm, 50 * cm)
    
    timeline = trial_df['interp_time'].values
    stim_raster = trial_df['stim_raster'].values
    spike_raster = trial_df['spike_raster'].values
    flicker_start = trial_df[trial_df['flicker_raster'] == 1]['interp_time'].values[0]
    
    timeline = timeline - flicker_start
    trace_sbr = trial_df['detrend_trace'].values/trial_df['trace_noise'].values
    axs[0].plot(timeline, trace_sbr, 'k', linewidth=0.3)

    # Plot the spikes detected in the trial
    spike_idx = np.where(spike_raster == 1)[0]
    axs[0].plot(timeline[spike_idx], (np.max(trace_sbr) + 1)*np.ones_like(spike_idx), '.k', markersize=2)

    # Plotting the stimulation indicators
    axs[0].plot(timeline, trial_df['flicker_raster'] + 3 + np.max(trace_sbr), '-b', linewidth=0.3)
    axs[0].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + np.max(trace_sbr)), color=consts.pulse_color, linewidth=0.3)
    stim_idx = np.where(stim_raster == 1)[0]
    axs[0].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + np.max(trace_sbr)), '|',\
             linewidth=20,  color=consts.pulse_color)

    # Plot the line bars
    t_bar_length = 500# time in ms
    #axs[0].plot([timeline[0], timeline[0] + t_bar_length/1000], [np.min(trace_sbr) - 0.2, np.min(trace_sbr) - 0.2], '-k')
    #axs[0].text(timeline[0], np.min(trace_sbr) - 1.1, str(t_bar_length) + ' ms', color='k')

    sbr_length = 5
    axs[0].plot([timeline[0] + 0.2, timeline[0] + 0.2], [0, sbr_length], '-k')
    axs[0].text(timeline[0] + 0.1, 0, str(sbr_length) + ' SBR', rotation=90)

    axs[0].spines['top'].set_visible(False)
    axs[0].spines['left'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].set_yticks([])

        # Plotting the stimulation indicators
    axs[1].plot(timeline, trial_df['flicker_raster'] +3+ exam_df['trial_id'].nunique(), '-b', linewidth=0.3)
    axs[1].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + exam_df['trial_id'].nunique()), color=consts.pulse_color, linewidth=0.3)
    stim_idx = np.where(stim_raster == 1)[0]
    axs[1].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + exam_df['trial_id'].nunique()), '|',\
             linewidth=20,  color=consts.pulse_color)

    # Plot the trials heatpmap
    exam_df['trace_sbr'] = exam_df['detrend_trace'] / exam_df['trace_noise']
    plot_df = exam_df.pivot(index='trial_id', columns='interp_time', values='trace_sbr')
    
    
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    surf_p = axs[1].pcolormesh(timeline, 1 + np.arange(0, plot_df.index.nunique()), plot_df.values)
    axs[1].set_ylabel('Trail #')

    # Plot the trial average trace with error bar
    avg_vm = plot_df.mean(axis=0)
    num_trials = plot_df.shape[0]
    std_vm = plot_df.std(axis=0)
    sem_vm = std_vm /np.sqrt(num_trials)
    axs[2].fill(np.concatenate((timeline, timeline[::-1])), \
                np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
                color='gray')
    axs[2].plot(timeline, avg_vm, '-k', linewidth=0.3)
    axs[2].spines['top'].set_visible(False)
    axs[2].spines['right'].set_visible(False)
    axs[2].spines['left'].set_visible(False)
    axs[2].set_yticks([])
    axs[2].set_xlabel('Time from Flicker Stim (S)')

    # Plotting the stimulation indicators
    axs[2].plot(timeline, trial_df['flicker_raster'] + 3 + np.max(avg_vm), '-b', linewidth=0.3)
    axs[2].plot([timeline[0], timeline[-1]], np.ones((2))* (2 + np.max(avg_vm)), color=consts.pulse_color, linewidth=0.3)
    stim_idx = np.where(stim_raster == 1)[0]
    axs[2].plot(timeline[stim_idx], np.ones_like(stim_idx)*(2 + np.max(avg_vm)), '|',\
             linewidth=20,  color=consts.pulse_color)

    # Plot the SBR bar
    sbr_length = 5
    axs[2].plot([timeline[0] + 0.2, timeline[0] + 0.2], [0, sbr_length], '-k')
    axs[2].text(timeline[0] + 0.1, 0, str(sbr_length) + ' SBR', rotation=90)
    plt.xlim(-1, 4)

    # Add the colorbar to the side of the heatmap
    cbar_ax = fig.add_axes((0.92, 0.35, 0.02, 0.35))
    cbar = fig.colorbar(surf_p, label='SBR', cax=cbar_ax)
    save_filename = savefig_path + exam_df['mouse_id'].unique()[0] +\
                '_' + exam_df['session_id'].unique()[0] +\
                '_fov' + str(exam_df['fov_id'].unique()[0]) +\
                '_roi' + str(exam_df['roi_id'].unique()[0]) +\
                '_' + str(exam_df['stim_freq'].unique()[0]) +\
                'Hz'
    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.pdf', format='pdf')
    plt.savefig(save_filename + '.png', format='png')
    
    #TODO need to fix the zoom in for the subplots, not the colorbar
    axs[1].set_xlim(-0.1, 1.1)
    axs[2].set_xlim(-0.1, 1.1)
    
    plt.savefig(save_filename + '_zoomIn.png', format='png')

    plt.show()

# %%
print(exam_df['trial_id'].nunique())

# %%
plt.close('all')


