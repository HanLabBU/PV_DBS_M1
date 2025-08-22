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
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'Compact' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %% Population Vm average

# Parameters for plot font
# Cosmetic parameters
plt.rcParams['font.size'] = 7
plt.rcParams['pdf.fonttype'] = 42
#plt.rcParams['font.family'] = 'arial'


# Determine which neuron population to plot here
nr_pop = 'all'
#nr_pop = 'inc_power'
#nr_pop = 'dec_power'
zsco = 0

# Loop through each ROI
test_freqs = df['stim_freq'].unique()
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Store all of the neurons trial-averaged vm
    vm_df = pd.DataFrame()

    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue
    
        # Determine which neuron population to analyze
        if nr_pop == 'inc_power' and fov_df['Power_Delta'].unique() == 0:
            continue
        elif nr_pop == 'dec_power' and fov_df['Power_Delta'].unique() == 1:
            continue

        flicker_start = fov_df[fov_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials normalized spike amplitude and set it as a new column
        # in the FOV dataframe
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            fov_df.loc[fov_df['trial_id'] == tr_val, 'interp_norm_vm'] = trial_df['interp_subvm'].values/avg_sp_amp
            
        avg_trial = fov_df.groupby('interp_time')['interp_norm_vm'].mean()

        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_trial.values,
            'interp_time':timeline
        }
        vm_df = pd.concat([vm_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        vm_df['neuron'] = vm_df['neuron'].astype('category')
        roi_acum += 1
    
    

    # Plot the trial averaged Vm
    plot_df = vm_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    
    # Determine to zscore
    if zsco == 1:
        plot_df = plot_df.apply(stats.zscore, axis=1)
    
    avg_vm = plot_df.mean(axis=0)
    std_vm = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_vm = std_vm/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(8.5 * cm, 3.5  * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
    plt.plot(timeline, avg_vm, '-k', linewidth=0.3)
    
    # Plot the flicker protocol
    plt.plot(timeline, 0.07*trial_df['flicker_raster']  + np.max(avg_vm + sem_vm) + .01, '-b')
    
    # Plot the stimulation protocol
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.22 + np.max(avg_vm + sem_vm)),\
              linewidth=0.3, color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.22 + np.max(avg_vm + sem_vm)), '|',\
         linewidth=0.3,  color=consts.pulse_color)
    
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Normalized Vm', labelpad=0)
    plt.title(nr_pop + ' Flicker Average: ' + str(freq) + ' Zscore: ' + str(bool(zsco)))
    
    ax = plt.gca()

    # Remove cluttered ticks
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Set the spacing for y-axis text
    ax.tick_params(axis='y', pad=0)
    ax.set_xlim(np.min(timeline), np.max(timeline))

    plt.tight_layout()

    save_filename = savefig_path + 'Population_Vm_' +\
                    nr_pop + '_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    #plt.savefig(save_filename + '_zsco' + str(zsco) + '.svg', format='svg')
    #plt.savefig(save_filename + '_zsco' + str(zsco) + '.png', format='png')
    plt.savefig(save_filename + '_zsco' + str(zsco) + '.pdf', format='pdf')
    
    plt.show()

# %%
# Plot both the low pass filtered population Vm and high pass filtered Vm
# Loop through each ROI
cutoff_freq = [3.5, 4.5] # hz
samp_freq =  500 # Hz
sos_stop = signal.butter(N=2, Wn=cutoff_freq, btype='bandstop' , fs=samp_freq, output='sos')
sos_pass = signal.butter(N=2, Wn=cutoff_freq, btype='bandpass' , fs=samp_freq, output='sos')

test_freqs = df['stim_freq'].unique()
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Store all of the neurons trial-averaged vm
    vm_df = pd.DataFrame()

    vm_stop_4Hz_df = pd.DataFrame()
    vm_pass_4Hz_df = pd.DataFrame()
    
    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
              
        # Individual FOV bandstop filtered Vm
        fov_bandstop_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Individual FOV highpass filtered Vm
        fov_bandpass_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_bandstop_df.empty:
            continue
        flicker_start = fov_bandstop_df[fov_bandstop_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_bandstop_df[fov_bandstop_df['trial_id'] == fov_bandstop_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials normalized spike amplitude and set it as a new column
        # in the FOV dataframe
        for tr_val in fov_bandstop_df['trial_id'].unique():
            trial_df = fov_bandstop_df[fov_bandstop_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            fov_bandstop_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'] = signal.sosfiltfilt(sos_stop, trial_df['interp_subvm'].values/avg_sp_amp)
            fov_bandpass_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'] = signal.sosfiltfilt(sos_pass, trial_df['interp_subvm'].values/avg_sp_amp)

            timeline = np.arange(np.array(trial_df['interp_subvm'].values/avg_sp_amp).shape[0])

            # DEBUG
            #fig, axs = plt.subplots(2, sharex=True)
            #axs[0].plot(timeline, fov_bandstop_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'], color='r', label='Band Stop Filtered')
            #axs[0].plot(timeline, trial_df['interp_subvm'].values/avg_sp_amp, color='b', label='Original Trace')
            #axs[0].legend()
#
            #axs[0].set_title('Bandstop Filtered Trace Overlay')
#
            
#
            #axs[1].plot(timeline, fov_bandpass_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'], color='r', label='Band Pass Filtered')
            #axs[1].plot(timeline, trial_df['interp_subvm'].values/avg_sp_amp, color='b', label='Original Trace')
            #axs[1].legend()
            #axs[1].set_title('Bandpass Filtered Trace Overlay')
#
            #plt.show()
            
        avg_stop_trial = fov_bandstop_df.groupby('interp_time')['interp_norm_vm'].mean()
        avg_pass_trial = fov_bandpass_df.groupby('interp_time')['interp_norm_vm'].mean()

        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_stop_trial.values,
            'interp_time':timeline
        }

        vm_stop_4Hz_df = pd.concat([vm_stop_4Hz_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        vm_stop_4Hz_df['neuron'] = vm_stop_4Hz_df['neuron'].astype('category')
        
        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_pass_trial.values,
            'interp_time':timeline
        }

        vm_pass_4Hz_df = pd.concat([vm_pass_4Hz_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        vm_pass_4Hz_df['neuron'] = vm_pass_4Hz_df['neuron'].astype('category')
        
        roi_acum += 1
    
    # Plot the trial averaged Vm bandstop
    plot_df = vm_stop_4Hz_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    avg_vm = plot_df.mean(axis=0)
    std_vm = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_vm = std_vm/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(65 * cm, 20 * cm)) #figsize=(25 * cm, 10 * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
    plt.plot(timeline, avg_vm, '-k', linewidth=0.75)
    
    # Plot the stimultion protocol
    plt.plot(timeline, 0.07*trial_df['flicker_raster'] + .10 + np.max(avg_vm + sem_vm), '-b')
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.01 + np.max(avg_vm + sem_vm)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.01 + np.max(avg_vm + sem_vm)), '|',\
         linewidth=20,  color=consts.pulse_color)
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Normalized Vm')
    plt.title('Flicker Average with 4 Hz Removal: ' + str(freq))
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_filename = savefig_path + 'Population_Vm_bandstop_4Hz_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')

    # Plot the trial averaged Vm bandpass
    plot_df = vm_pass_4Hz_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    avg_vm = plot_df.mean(axis=0)
    std_vm = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_vm = std_vm/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(25 * cm, 10 * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
    plt.plot(timeline, avg_vm, '-k', linewidth=0.75)
    
    # Plot the stimultion protocol
    plt.plot(timeline, 0.07*trial_df['flicker_raster'] + .10 + np.max(avg_vm + sem_vm), '-b')
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.01 + np.max(avg_vm + sem_vm)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.01 + np.max(avg_vm + sem_vm)), '|',\
         linewidth=20,  color=consts.pulse_color)
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Normalized Vm')
    plt.title('Flicker Average of just 4 Hz: ' + str(freq))
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_filename = savefig_path + 'Population_Vm_bandpass_4Hz_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')
    
plt.show()

# %%
# Plot the 8 Hz filtered Vm
# Loop through each ROI
cutoff_freq = [7.5, 8.5] # hz
samp_freq =  500 # Hz

sos_stop = signal.butter(N=2, Wn=cutoff_freq, btype='bandstop' , fs=samp_freq, output='sos')
sos_pass = signal.butter(N=2, Wn=cutoff_freq, btype='bandpass' , fs=samp_freq, output='sos')

single_nr_savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker' +f+ 'Vm' +f+ 'Individual' +f


test_freqs = df['stim_freq'].unique()
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Store all of the neurons trial-averaged vm
    vm_df = pd.DataFrame()

    vm_stop_8Hz_df = pd.DataFrame()
    vm_pass_8Hz_df = pd.DataFrame()
    
    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
              
        # Individual FOV bandstop filtered Vm
        fov_bandstop_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Individual FOV bandpass filtered Vm
        fov_bandpass_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_bandstop_df.empty:
            continue
        flicker_start = fov_bandstop_df[fov_bandstop_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_bandstop_df[fov_bandstop_df['trial_id'] == fov_bandstop_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials normalized spike amplitude and set it as a new column
        # in the FOV dataframe
        for tr_val in fov_bandstop_df['trial_id'].unique():
            trial_df = fov_bandstop_df[fov_bandstop_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)

            if np.isnan(avg_sp_amp):
                # Set spike amplitude to 1
                avg_sp_amp = 1

            fov_bandstop_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'] = signal.sosfiltfilt(sos_stop, trial_df['detrend_trace'].values/avg_sp_amp)
            fov_bandpass_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'] = signal.sosfiltfilt(sos_pass, trial_df['detrend_trace'].values/avg_sp_amp)


            # DEBUG
            #fig, axs = plt.subplots(2, sharex=True)
            #timeline = np.arange(np.array(trial_df['interp_subvm'].values/avg_sp_amp).shape[0])
            #axs[0].plot(timeline, fov_bandstop_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'], color='r', label='Band Stop Filtered')
            #axs[0].plot(timeline, trial_df['interp_subvm'].values/avg_sp_amp, color='b', label='Original Trace')
            #axs[0].legend()
#
            #axs[0].set_title('Bandstop Filtered Trace Overlay')
#
            
#
            #axs[1].plot(timeline, fov_bandpass_df.loc[fov_bandstop_df['trial_id'] == tr_val, 'interp_norm_vm'], color='r', label='Band Pass Filtered')
            #axs[1].plot(timeline, trial_df['interp_subvm'].values/avg_sp_amp, color='b', label='Original Trace')
            #axs[1].legend()
            #axs[1].set_title('Bandpass Filtered Trace Overlay')
#
            #plt.show()
            
        avg_stop_trial = fov_bandstop_df.groupby('interp_time')['interp_norm_vm'].mean()
        avg_pass_trial = fov_bandpass_df.groupby('interp_time')['interp_norm_vm'].mean()

        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_stop_trial.values,
            'interp_time':timeline
        }

        vm_stop_8Hz_df = pd.concat([vm_stop_8Hz_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        vm_stop_8Hz_df['neuron'] = vm_stop_8Hz_df['neuron'].astype('category')
        
        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_pass_trial.values,
            'interp_time':timeline
        }

        vm_pass_8Hz_df = pd.concat([vm_pass_8Hz_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        vm_pass_8Hz_df['neuron'] = vm_pass_8Hz_df['neuron'].astype('category')
        
        roi_acum += 1

        #--Plotting individual filtered trials with the final Vm average
        cm = 1/2.54
        fig, axs = plt.subplots(2, 1, figsize=(3*10.5*cm, 3*10*cm))
        fig.patch.set_facecolor('white')
        # Transform matrix to have plotable trials
        trial_vms = np.transpose(fov_bandpass_df.pivot(index='trial_id', columns='interp_time', values='interp_norm_vm').values)

        # This normalization is for plotting purposes
        norm_pl_vms = consts.norm_signals(trial_vms)

        # Loop and plot each column
        for i, tr in enumerate(np.transpose(norm_pl_vms)):
            if np.isnan(tr[0]):
                print('Hello!!')
                raise Exception('Trial does not have the correct shape')

            axs[0].plot(timeline, i + tr, '-k')
        # Plot the flicker protocol
        axs[0].plot(timeline, 1*trial_df['flicker_raster']  + i + 1.5, '-b')
        
        # Plot the electrical stim protocol
        axs[0].plot([timeline[0], timeline[-1]], np.ones((2))* (1.1 + i),\
              linewidth=0.3, color=consts.pulse_color)
        stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
        axs[0].plot(timeline[stim_idx], np.ones_like(stim_idx)*(1.1 + i), '|',\
            linewidth=0.3,  color=consts.pulse_color)
       
        # Plot the trial-average Vm
        axs[1].plot(timeline, avg_pass_trial.values)

        # Plot the flicker protocol
        axs[1].plot(timeline, .05*trial_df['flicker_raster']  + np.max(avg_pass_trial.values), '-b')
        
        # Plot the electrical stim protocol
        axs[1].plot([timeline[0], timeline[-1]], np.ones((2))* (np.max(avg_pass_trial.values)),\
              linewidth=0.3, color=consts.pulse_color)
        stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
        axs[1].plot(timeline[stim_idx], np.ones_like(stim_idx)*(np.max(avg_pass_trial.values)), '|',\
            linewidth=0.3,  color=consts.pulse_color)

        # Remove axis
        axs[0].spines['top'].set_visible(False)
        axs[0].spines['right'].set_visible(False)

        axs[1].spines['top'].set_visible(False)
        axs[1].spines['right'].set_visible(False)

        axs[0].set_title("_".join([str(e) for e in values]) +"_"+ str(freq) + ' Hz ')
        
        save_filename = single_nr_savefig_path + '8Hz_filt_' + "_".join([str(e) for e in values]) +"_"+\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'
        plt.savefig(save_filename + '.png', format='png')
        

    #DEBUG Plot all of the filtered individual
    #raise Exception('After plotting individual neurons stuff')

    
    #--Plot the neuron averaged Vm bandstop
    plot_df = vm_stop_8Hz_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    avg_vm = plot_df.mean(axis=0)
    std_vm = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_vm = std_vm/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(65 * cm, 20 * cm)) #figsize=(25 * cm, 10 * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
    plt.plot(timeline, avg_vm, '-k', linewidth=0.75)
    
    # Plot the flicker protocol
    plt.plot(timeline, 0.07*trial_df['flicker_raster'] + .01 + np.max(avg_vm + sem_vm), '-b')

    # Plot the stimulation protocol
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.10 + np.max(avg_vm + sem_vm)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.10 + np.max(avg_vm + sem_vm)), '|',\
         linewidth=20,  color=consts.pulse_color)
    
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Normalized Vm')
    plt.title('Flicker Average with 8 Hz Removal, stim Freq:' + str(freq) + ' Hz')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_filename = savefig_path + 'Population_Vm_bandstop_8Hz_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')

    #--Plot the neuron averaged Vm bandpass
    plot_df = vm_pass_8Hz_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    avg_vm = plot_df.mean(axis=0)
    std_vm = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_vm = std_vm/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(25 * cm, 10 * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
            color='gray')
    plt.plot(timeline, avg_vm, '-k', linewidth=0.75)
    
    # Plot the flicker protocol
    plt.plot(timeline, 0.07*trial_df['flicker_raster'] + .10 + np.max(avg_vm + sem_vm), '-b')
    
    # Plot the stimulation protocol
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.010 + np.max(avg_vm + sem_vm)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.010 + np.max(avg_vm + sem_vm)), '|',\
         linewidth=20,  color=consts.pulse_color)
    
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Normalized Vm')
    plt.title('Flicker Average of just 8 Hz, stim freq:' + str(freq) + ' Hz')
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_filename = savefig_path + 'Population_Vm_bandpass_8Hz_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')
    
plt.show()

# %%
# Calculate the firing rate for each neuron
spike_rate_win = 50
samp_freq = 500
test_freqs = df['stim_freq'].unique()

# Determine which neuron population to plot here
#nr_pop = 'all'
#nr_pop = 'inc_power'
nr_pop = 'dec_power'

for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Store all of the neurons trial-averaged vm
    fr_df = pd.DataFrame()

    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue

        # Determine which neuron population to analyze
        if nr_pop == 'inc_power' and fov_df['Power_Delta'].unique() == 0:
            continue
        elif nr_pop == 'dec_power' and fov_df['Power_Delta'].unique() == 1:
            continue
        
        flicker_start = fov_df[fov_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials firing rate and set it as a new column
        # in the FOV dataframe
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            
            firing_rate = uniform_filter1d(trial_df['interp_spike_raster'].values*samp_freq, size=25)

            #TODO save the firing rate instead of this Vm stuff
            fov_df.loc[fov_df['trial_id'] == tr_val, 'firing_rate'] = firing_rate

        avg_trial = fov_df.groupby('interp_time')['firing_rate'].mean()

        nr_dict = {
            'neuron':roi_acum,
            'fr_trial_avg':avg_trial.values,
            'interp_time':timeline
        }
        fr_df = pd.concat([fr_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        fr_df['neuron'] = fr_df['neuron'].astype('category')
        roi_acum += 1
    
    # Plot the trial averaged firing rate
    plot_df = fr_df.pivot(index='neuron', columns='interp_time', values='fr_trial_avg')
    avg_fr = plot_df.mean(axis=0)
    std_fr = plot_df.std(axis=0)
    num_neurons = plot_df.shape[0]
    sem_fr = std_fr/np.sqrt(num_neurons)

    cm = 1/2.54
    plt.figure(figsize=(25.5 * cm, 10 * cm))
    plt.fill(np.concatenate((timeline, timeline[::-1])), \
            np.concatenate((avg_fr + sem_fr, avg_fr[::-1] - sem_fr[::-1])),
            color='gray')
    plt.plot(timeline, avg_fr, '-k')
    
    # Plot the stimultion protocol
    plt.plot(timeline, 1*trial_df['flicker_raster'] + 1 + np.max(avg_fr + sem_fr), '-b')
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (.01 + np.max(avg_fr + sem_fr)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(.01 + np.max(avg_fr + sem_fr)), '|',\
         linewidth=20,  color=consts.pulse_color)
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Firing Rate (Hz)')
    plt.title(nr_pop + ' Flicker Average: ' + str(freq))
    
    ax = plt.gca()
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    save_filename = savefig_path + 'Population_FR_' +\
                    nr_pop + '_' +\
                    str(stim_df['stim_freq'].unique()[0]) + 'Hz'

    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')
    
    plt.show()

# %%
plt.close('all')

# %%
# check the padding
raster = np.zeros(150)
raster[(150*np.random.rand(30)).astype(int)] = 1

plt.figure()
plt.plot(raster, label='original')
padded_raster = np.pad(raster, (25, 25), mode='edge')

plt.plot(padded_raster + 0.5, label='padded')
mv_avg = uniform_filter1d(padded_raster, size=25, mode='constant')
plt.plot(mv_avg + 1, label='Firing Rate')
plt.legend()
plt.show()


