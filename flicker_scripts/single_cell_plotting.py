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
savefig_path = f + 'home' +f+ 'pierfier' +f+ 'Dropbox' +f+ 'RKC-HanLab' +f+ 'Pierre PV DBS Project Dropbox' +f+ 'Materials' +f+ 'Plots' +f+ 'Flicker'
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %%
# Find neurons that are significantly depolarized right after the flicker
wind = 100/1000 # in ms

# Determine whether or not to z-score
zsco = 1

test_freqs = df['stim_freq'].unique() # [40]
for freq in test_freqs:
    stim_df = df[df['stim_freq'] == freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    non_empty_fovs = [pair for pair in pairings\
                      if not stim_df[(stim_df['mouse_id'] == pair[0]) &\
                        (stim_df['session_id'] == pair[1]) &\
                        (stim_df['fov_id'] == pair[2])].empty ]
    
    # Initialize the empty dataframes
    act_vm_df = pd.DataFrame()
    sup_vm_df = pd.DataFrame()
    non_mod_vm_df = pd.DataFrame()
    
    roi_acum = 0
    
    # Loop through each neuron
    for values in non_empty_fovs:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                         (stim_df['fov_id'] == values[2])]
        
        # Time for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values

        # A mask for the neuron for the giant dataframe
        nr_mask = (df['stim_freq'] == freq) & \
                    (df['mouse_id'] == values[0]) & \
                    (df['session_id'] == values[1]) & \
                    (df['fov_id'] == values[2])

        # Calcualte the trial average Vm
        avg_trial = fov_df.groupby('interp_time')['interp_subvm'].mean()
        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_trial.values,
            'interp_time':timeline
        }

        # Get the Vm during flicker that is within the window right after the onset
        flicker_onset_mask = (fov_df['interp_time'] > 0) & (fov_df['interp_time'] <= 0 + wind)
        baseline_mask = (fov_df['interp_time'] > 0 - wind) & (fov_df['interp_time'] <= 0)

        # Calculate the average baseline Vm
        avg_base_vm_df = fov_df[baseline_mask].groupby('trial_id')['interp_subvm'].mean()
        avg_flicker_vm_df = fov_df[flicker_onset_mask].groupby('trial_id')['interp_subvm'].mean()

        stat, p_value = stats.wilcoxon(avg_base_vm_df.values, avg_flicker_vm_df.values)

        # Label neuron as significantly depolarized
        if p_value < 0.05 and np.mean(avg_base_vm_df.values) < np.mean(avg_flicker_vm_df.values):
            df.loc[nr_mask, 'Vm_flick_mod'] = 1
            act_vm_df = pd.concat([act_vm_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        elif p_value < 0.05 and np.mean(avg_base_vm_df.values) > np.mean(avg_flicker_vm_df.values):
            df.loc[nr_mask, 'Vm_flick_mod'] = -1
            sup_vm_df = pd.concat([sup_vm_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        else:
            df.loc[nr_mask, 'Vm_flick_mod'] = 0
            non_mod_vm_df = pd.concat([non_mod_vm_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        
        roi_acum += 1

    # Pivot all of the dataframes
    try:
        act_vm_df = act_vm_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    except:
        pass
    try:
        non_mod_vm_df = non_mod_vm_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    except:
        pass
    try:
        sup_vm_df = sup_vm_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    except:
        pass

    # Zscore each population
    if zsco == 1:
        act_vm_df = act_vm_df.apply(stats.zscore, axis=1)
        non_mod_vm_df = non_mod_vm_df.apply(stats.zscore, axis=1)
        sup_vm_df = sup_vm_df.apply(stats.zscore, axis=1)

    try:
        # Sort the activated neurons
        flick_time = (act_vm_df.columns > 0) & (act_vm_df.columns < 1)
        act_vm_df['flick_avg'] = act_vm_df.loc[:, flick_time].mean(axis=1)
        act_vm_df_sorted = act_vm_df.sort_values(by='flick_avg', ascending=True)
        act_vm_df_sorted = act_vm_df_sorted.drop(columns=['flick_avg'])
    except:
        act_vm_df_sorted = pd.DataFrame()

    try:    
        # Sort the non-modulated neurons
        flick_time = (non_mod_vm_df.columns > 0) & (non_mod_vm_df.columns < 1)
        non_mod_vm_df['flick_avg'] = non_mod_vm_df.loc[:, flick_time].mean(axis=1)
        non_mod_vm_df_sorted = non_mod_vm_df.sort_values(by='flick_avg', ascending=True)
        non_mod_vm_df_sorted = non_mod_vm_df_sorted.drop(columns=['flick_avg'])
    except:
        non_mod_vm_df_sorted = pd.DataFrame()

    try:
        # Sort the suppressed neurons
        flick_time = (sup_vm_df.columns > 0) & (sup_vm_df.columns < 1)
        sup_vm_df['flick_avg'] = sup_vm_df.loc[:, flick_time].mean(axis=1)
        sup_vm_df_sorted = sup_vm_df.sort_values(by='flick_avg', ascending=True)
        sup_vm_df_sorted = sup_vm_df_sorted.drop(columns=['flick_avg'])
    except:
        sup_vm_df_sorted = pd.DataFrame()

    # actine all of the neurons into a printable dataframe
    comb_vm_df = pd.concat([act_vm_df_sorted, non_mod_vm_df_sorted, sup_vm_df_sorted], ignore_index=True, join='outer')
    #comb_vm_df = act_vm_df_sorted

    # Plot all of the Vms
    cm = 1/2.54
    plt.figure(plt.figure(figsize=(20.5 * cm, 15 * cm)))
    surf_h = plt.pcolormesh(timeline, np.arange(1, comb_vm_df.shape[0] + 1), comb_vm_df.values)
    plt.colorbar(surf_h, label='SBR zscore: ' + str(zsco))
    plt.xlabel('Time (S)')
    plt.ylabel('Neuron')
    plt.vlines(x=[1, 2], ymin=0.5, ymax=comb_vm_df.shape[0] + 0.5, linestyles='--', colors=consts.stim_color)
    plt.title('Neuron trial-averaged Firing Rate ' + str(freq) + ' Separated')

    # Separate each neuron population based on sizes
    num_act = 0
    num_non = 0
    num_sup = 0
    if not act_vm_df.empty:
        num_act = act_vm_df.shape[0]
    if not non_mod_vm_df.empty:
        num_non = non_mod_vm_df.shape[0]
    if not sup_vm_df.empty:
        num_sup = sup_vm_df.shape[0]
    
    plt.hlines(y=[0.5 + num_act, 
                  0.5 + num_act + num_non,
                  0.5 + num_act + num_non + num_sup],
                xmin=-1, xmax=np.max(timeline), colors='k')

    # Set Xlimits based on timeline min max
    plt.xlim(np.min(timeline)), np.max(timeline)

    # Save figure
    save_filename = savefig_path +f+ 'Neuronwise' +f+ 'Trial_Avg_Vm_' + str(freq) + 'Hz_Separated'
    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')
# %% 
# Perform single cell analysis Vm

# Flag to determine whether to zscore or not
zsco = 1

#nr_pop = 'all'
#nr_pop = 'separate'

#nr_pop = 'inc_power'
#nr_pop = 'dec_power'

# Population with/out 8 Hz entrainment change from base to flicker
#nr_pop = '8entrain'
#nr_pop = 'non_8entrain'

# Population with/out plv entrainment
#nr_pop = 'plv_etrain'
#nr_pop = 'plv_nonetrain'



# Loop through each stim frequency
test_freqs = df['stim_freq'].unique() # [40]
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
        elif nr_pop == '8entrain' and fov_df['8entrain'].unique() == 0:
            continue
        elif nr_pop == 'non_8entrain' and fov_df['8entrain'].unique() == 1:
            continue
        elif nr_pop == 'plv_etrain' and fov_df['plv_etrain'].unique() == 0:
            continue
        elif nr_pop == 'plv_nonetrain' and fov_df['plv_etrain'].unique() == 1:
            continue

        flicker_start = fov_df[fov_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start

        # Calculate the trial average Vm for each neuron
        for tr_val, trial_df in fov_df.groupby('trial_id'):
            # Get the average spike height
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            if np.isnan(avg_sp_amp):
                avg_sp_amp = 1
            
            fov_df.loc[fov_df['trial_id'] == tr_val, 'interp_norm_vm'] = trial_df['interp_subvm'].values/avg_sp_amp

        # Calculate average and add to the total population
        avg_trial = fov_df.groupby('interp_time')['interp_norm_vm'].mean()
        nr_dict = {
            'neuron':roi_acum,
            'vm_trial_avg':avg_trial.values,
            'interp_time':timeline
        }
        vm_df = pd.concat([vm_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        roi_acum += 1
    
    # Shorten the memory usage for neuron tag
    vm_df['neuron'] = vm_df['neuron'].astype('category')

    # Pivot and sort the Vms by the 8Hz entrainment value
    vm_df = vm_df.pivot(index='neuron', columns='interp_time', values = 'vm_trial_avg')

    # Z-score the Vms for each neuron
    if zsco == 1:
        vm_df = vm_df.apply(stats.zscore, axis=1)

    # Average the flicker component and sort by that value
    flick_time = (vm_df.columns > 0) & (vm_df.columns < 1)
    vm_df['flick_avg'] = vm_df.loc[:, flick_time].mean(axis=1)
    vm_df_sorted = vm_df.sort_values(by='flick_avg', ascending=True)
    vm_df_sorted = vm_df_sorted.drop(columns=['flick_avg'])

    # Plot all of the Vms
    cm = 1/2.54
    plt.figure(plt.figure(figsize=(20.5 * cm, 15 * cm)))
    surf_h = plt.pcolormesh(timeline, np.arange(1, vm_df_sorted.shape[0] + 1), vm_df_sorted.values)
    plt.colorbar(surf_h, label='SBR zscore: ' + str(zsco))
    plt.xlabel('Time (S)')
    plt.ylabel('Neuron')
    plt.vlines(x=[1, 2], ymin=0.5, ymax=vm_df_sorted.shape[0] + 0.5, linestyles='-', colors=consts.stim_color)
    plt.title('Neuron trial-averaged Vm ' + str(freq) + ' ' + nr_pop)

    # Save figure
    save_filename = savefig_path +f+ 'Neuronwise' +f+ 'Trial_Avg_Vm_' + str(freq) + 'Hz_' + nr_pop
    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')

# %% 
# Perform single cell analysis Firing Rate

# Flag to determine whether to zscore or not
zsco = 1

samp_freq = 500

nr_pop = 'all'
#nr_pop = 'inc_power'
#nr_pop = 'dec_power'

# Loop through each stim frequency
test_freqs = df['stim_freq'].unique() # TODO change to 
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

        # Calculate the trial average Vm for each neuron
        for tr_val, trial_df in fov_df.groupby('trial_id'):
            firing_rate = uniform_filter1d(trial_df['interp_spike_raster'].values*samp_freq, size=25)           
            fov_df.loc[fov_df['trial_id'] == tr_val, 'firing_rate'] = firing_rate

        # Calculate average and add to 
        avg_trial = fov_df.groupby('interp_time')['firing_rate'].mean()
        nr_dict = {
            'neuron':roi_acum,
            'fr_trial_avg':avg_trial.values,
            'interp_time':timeline
        }
        fr_df = pd.concat([fr_df, pd.DataFrame(nr_dict)], ignore_index=True, join='outer')
        roi_acum += 1
    
    # Shorten the memory usage for neuron tag
    fr_df['neuron'] = fr_df['neuron'].astype('category')

    # Pivot and sort the Vms by the 8Hz entrainment value
    fr_df = fr_df.pivot(index='neuron', columns='interp_time', values = 'fr_trial_avg')

    # Z-score the Vms for each neuron
    if zsco == 1:
        fr_df = fr_df.apply(stats.zscore, axis=1)

    # Average the flicker component and sort by that value
    flick_time = (fr_df.columns > 0) & (fr_df.columns < 1)
    fr_df['flick_avg'] = fr_df.loc[:, flick_time].mean(axis=1)
    fr_df_sorted = fr_df.sort_values(by='flick_avg', ascending=True)
    fr_df_sorted = fr_df_sorted.drop(columns=['flick_avg'])

    # Plot all of the Vms
    cm = 1/2.54
    plt.figure(plt.figure(figsize=(20.5 * cm, 15 * cm)))
    surf_h = plt.pcolormesh(timeline, np.arange(1, fr_df_sorted.shape[0] + 1), fr_df_sorted.values)
    plt.colorbar(surf_h, label='SBR zscore: ' + str(zsco))
    plt.xlabel('Time (S)')
    plt.ylabel('Neuron')
    plt.vlines(x=[1, 2], ymin=0.5, ymax=fr_df_sorted.shape[0] + 0.5, linestyles='-', colors=consts.stim_color)
    plt.title('Neuron trial-averaged Firing Rate ' + str(freq) + ' ' + nr_pop)

    # Save figure
    save_filename = savefig_path +f+ 'Neuronwise' +f+ 'Trial_Avg_FR_' + str(freq) + 'Hz_' + nr_pop
    plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')

# %% Close all of the figures
plt.close('all')