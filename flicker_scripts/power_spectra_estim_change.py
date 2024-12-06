# This script is meant to perform the regression of the change in 8 Hz power from onset to ESTIM period
# to increase in stimulation frequency power

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
import statsmodels.api as sm

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
          +f+ 'Flicker' +f+ 'Spectra' +f
df = pd.read_pickle(interm_data_path)
print(df.columns)

# %%
# Functions to help with violin plotting

# Returns a string representation for 
def get_violin_opts():
    matlab_exp = "opts.ShowMedian = true;\
            opts.ShowMean = false;\
            opts.MedianColor = [1 1 1];\
            opts.MarkerSize = 5;\
            opts.MedianMarkerSize = 10;\
            opts.BoxWidth = 0.1;\
            opts.BoxColor = [0 0 0];\
            opts.ViolinAlpha = {[0.1], [0.1]};"
    return matlab_exp

# Using the matlab violin plot stuff
#def specify_violin_color(violin_i, color):
#    pass

# %%
# Population average power spectra

# Reset stats file
open(savefig_path + 'stats.txt', 'w').close()

# Neuron population with 8 hz change between flicker and DBS 
nr_pop = 'all'
#nr_pop = 'inc_power'
#nr_pop = 'dec_power'

# Population with/out 8 Hz entrainment change from base to flicker
#nr_pop = '8entrain'
#nr_pop = 'non_8entrain'

# Population with/out plv entrainment
#nr_pop = 'plv_etrain'
#nr_pop = 'plv_nonetrain'

samp_freq = 500
freq_limit = [1, 200]
freq_nums = 3*(freq_limit[1] - freq_limit[0])
wavelet = 'morl'

# Has the spectra for all neurons for all conditions into 1 dataframe
all_spectra_df = pd.DataFrame()

test_freqs = df['stim_freq'].unique() # [40] #  [140] #    
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    non_empty_fovs = [pair for pair in pairings\
                      if not stim_df[(stim_df['mouse_id'] == pair[0]) &\
                        (stim_df['session_id'] == pair[1]) &\
                        (stim_df['fov_id'] == pair[2])].empty ]

    # Store all of the neurons trial-averaged vm
    vm_df = pd.DataFrame()

    power_spec_df = pd.DataFrame()

    # Loop through unique FOVs
    roi_acum = 0
    for values in non_empty_fovs:
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
        
        #Calculate each trials normalized spike amplitude Vm and set it as a new column
        # in the FOV dataframe
        nr_power_spec_df = pd.DataFrame()
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            norm_subvm = trial_df['detrend_trace'].values/avg_sp_amp # or 'sub_vm' Switching between the full Vm and the Sub Vm
            fov_df.loc[fov_df['trial_id'] == tr_val, 'norm_vm'] = norm_subvm

            # Calculate the trial power spectra
            calc_freqs, coeffs = fcwt.cwt(norm_subvm, samp_freq, freq_limit[0], freq_limit[1], freq_nums)
            power_vals = np.abs(coeffs)

            # Linearize power spectra to add into dataframe
            linear_power = power_vals.ravel()
            t_mesh, freq_mesh = np.meshgrid(timeline, calc_freqs)
            t_lin = t_mesh.ravel()
            freq_lin = freq_mesh.ravel()

            tr_dict = {
                'trial':tr_val,
                'power':linear_power,
                'timeline':t_lin,
                'freq_spec':freq_lin
            }

            nr_power_spec_df = pd.concat([nr_power_spec_df, pd.DataFrame(tr_dict)], ignore_index=True, join='outer')
            nr_power_spec_df['trial'] = nr_power_spec_df['trial'].astype('category')


        # Calculate the trial-averaged Vm
        avg_tr_power_df = nr_power_spec_df.groupby(['timeline', 'freq_spec'])['power'].mean()
        avg_tr_power_df = avg_tr_power_df.unstack(level='timeline')
        avg_tr_power_df.sort_index(ascending=False, inplace=True)
        nr_norm_vals_df = avg_tr_power_df

        #--Testing different normalizations
        #Z-score across frequencies 

        #nr_norm_vals_df = nr_norm_vals_df.apply(lambda col: zscore(col, ddof=0), axis=0)

        #Z-score across time
        #nr_norm_vals_df = nr_norm_vals_df.apply(lambda row: zscore(row, ddof=0), axis=1)
        
        # (x - B)/(A + B) normalization where 
        # Set the time limits for power normalization
        #base_limit = [0, 1] # Flicker only baseline period
        #stim_limit = [1, 2] # Estim period

        base_limit = [-1, 0] # Base only period 
        stim_limit = [0, 1] # Flicker only period

        # Calculate the flicker baseline average across frequencies
        base_flicker_time_mask = (nr_norm_vals_df.columns > base_limit[0]) & (nr_norm_vals_df.columns < base_limit[1])
        base_flicker_pow_df = nr_norm_vals_df[nr_norm_vals_df.columns[base_flicker_time_mask]]
        base_flicker_pow_avg = base_flicker_pow_df.mean(axis=1)

        # Calculate the stimulation period average across frequencies
        stim_time_mask = (nr_norm_vals_df.columns > stim_limit[0]) & (nr_norm_vals_df.columns < stim_limit[1])
        stim_pow_df = nr_norm_vals_df[nr_norm_vals_df.columns[stim_time_mask]]
        stim_pow_avg = stim_pow_df.mean(axis=1)

        # Perform the final formulaic normalization
        summed_periods = base_flicker_pow_avg + stim_pow_avg
        nr_norm_vals_df = (nr_norm_vals_df.sub(base_flicker_pow_avg, axis=0))\
                        .div(summed_periods, axis=0)

        # Setup the dataframe to be concatenated for the neuron population
        nr_norm_vals_df = nr_norm_vals_df.stack()
        nr_norm_vals_df.name = 'power'
        nr_norm_vals_df = nr_norm_vals_df.reset_index()
        
        #Originally was just an incremental number
        #nr_norm_vals_df['neuron'] = roi_acum
        #DEBUG using a neuron name
        nr_norm_vals_df['neuron'] = str(values[0]) + '_' + str(values[1]) + '_fov' + str(values[2]) + '_stim' + str(stim_freq)

        power_spec_df = pd.concat([power_spec_df, nr_norm_vals_df], ignore_index=True, join='outer')
        power_spec_df['neuron'] = power_spec_df['neuron'].astype('category')
        roi_acum += 1
    
    # Add all of the power spectra to the total spectra dataframe
    power_spec_temp_df = power_spec_df.copy()
    power_spec_temp_df['stim_freq'] = stim_freq
    all_spectra_df = pd.concat([all_spectra_df, power_spec_temp_df],\
                            ignore_index=True, join='outer')


    # Print out the number of neurons for each frequency
    print(stim_freq)
    print('Num neurons ' + str(roi_acum + 1))

    # Get the 2D average of the heatmap
    avg_power_spec = power_spec_df.groupby(['timeline', 'freq_spec'])['power'].mean()

    # Create a 2D dataframe for    
    heatmap_df = avg_power_spec.unstack(level='timeline')
            
    # Chop starting and end points to no longer display edge effects
    t_trim = 50 # 1 #  
    display_timeline = timeline[t_trim:-t_trim]
    t_trimmed_idx = (heatmap_df.columns >= display_timeline[0]) & (heatmap_df.columns <= display_timeline[-1])
    heatmap_df = heatmap_df[heatmap_df.columns[t_trimmed_idx]]
    
    cm = 1/2.54
    plt.figure(figsize=(10.5 * 10* cm, 7 * 10* cm))
    plt.rcParams['font.size'] = 40
    surf_p = plt.pcolormesh(display_timeline, heatmap_df.index.values, heatmap_df.values, cmap='jet') #, vmin=-1, vmax=1
        
    # Plot the stimultion protocol
    plt.plot(display_timeline,
             2*trial_df[t_trimmed_idx]['flicker_raster'] + 13 + np.max(heatmap_df.index.values), '-b')
    
    plt.plot([display_timeline[0], display_timeline[-1]], np.ones((2))* (10 + np.max(heatmap_df.index.values)), color=consts.pulse_color)
    stim_idx = np.where(trial_df[t_trimmed_idx]['stim_raster'].values == 1)[0]
    plt.plot(display_timeline[stim_idx], np.ones_like(stim_idx)*(10 + np.max(heatmap_df.index.values)), '|',\
         linewidth=20,  color=consts.pulse_color)
    
    # Remove the spines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)

    #Labels
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title(nr_pop + ' Spectra ' + str(stim_freq))

    # Add colorbar axis
    fig = plt.gcf()
    cbar_ax = fig.add_axes((0.92, 0.35, 0.02, 0.35))
    cbar = fig.colorbar(surf_p, label='Normalized Power Change', cax=cbar_ax)

    # Set the default colorbar scale similar to main figures
    surf_p.set_clim(-0.2, 0.8)

    # If color bar is too high, cap it
    #c_limits = surf_p.get_clim()
    #if c_limits[1] > 1.5:
    #    surf_p.set_clim(c_limits[0], 1.5)

    save_filename = savefig_path + 'Population_Flicker_Stim_Lesion_' +\
            nr_pop + '_' +\
            str(stim_df['stim_freq'].unique()[0]) + 'Hz_MaxFreq' + str(freq_limit[1])
    
    # Spectra SVGs are incredibly hard to edit in illustrator
    #plt.savefig(save_filename + '.svg', format='svg')
    plt.savefig(save_filename + '.png', format='png')
    plt.show()

    ## Plot the trial averaged Vm
    #plot_df = vm_df.pivot(index='neuron', columns='interp_time', values='vm_trial_avg')
    #avg_vm = plot_df.mean(axis=0)
    #std_vm = plot_df.std(axis=0)
    #num_neurons = plot_df.shape[0]
    #sem_vm = std_vm/np.sqrt(num_neurons)
    #plt.figure(figsize=(13, 4))
    #plt.fill(np.concatenate((timeline, timeline[::-1])), \
    #        np.concatenate((avg_vm + sem_vm, avg_vm[::-1] - sem_vm[::-1])),
    #        color='gray')
    #plt.plot(timeline, avg_vm, '-k')
    #
    #plt.xlabel('Time from flicker onset (s)')
    #plt.ylabel('Normalized Vm')
    #plt.title('Flicker Average: ' + str(freq))
    #plt.show()

    # Make violin plotsW

    # TODO create timeline masks for the below df and just use that as a way to better reference all of these 
    # also, add a violin plot for the baseline, flicker, stim+flicker, flicker post-stim, then do statistics between each of these periods
    mask_base_pre = (power_spec_df['timeline'] < 0)
    mask_flicker_pre = (power_spec_df['timeline'] >= 0) & (power_spec_df['timeline'] < 1)
    mask_stim = (power_spec_df['timeline'] >= 1) | (power_spec_df['timeline'] < 2)
    mask_flicker_post = (power_spec_df['timeline'] >= 2) & (power_spec_df['timeline'] < 3)
    
    mask_8hz_spec = (power_spec_df['freq_spec'] > 7.5) & (power_spec_df['freq_spec'] < 8.5)

    # Take the average for each neuron across time
    base_pre_power = power_spec_df[mask_base_pre & mask_8hz_spec].groupby(['neuron']).mean()
    flicker_pre_power = power_spec_df[mask_flicker_pre & mask_8hz_spec].groupby(['neuron']).mean()
    stim_power = power_spec_df[mask_stim & mask_8hz_spec].groupby(['neuron']).mean()
    flicker_post_power = power_spec_df[mask_flicker_post & mask_8hz_spec].groupby(['neuron']).mean()

    #flicker_only_period_data = flicker_only_period_data.groupby(['neuron']).mean()
    #dbs_flicker_data = dbs_flicker_data.groupby(['neuron']).mean()

    # (Deprecated, hopefully I do not need this) Determine if the neuron has an increase or decrease in flicker entrainment
    #for nr_name, row in base_period_data.iterrows():
    #    flicker_point = flicker_only_period_data.loc[nr_name, 'power']
    #    dbs_flicker_point = dbs_flicker_data.loc[nr_name, 'power']
#
    #    # Parse out neuron id fields
    #    vals = nr_name.split('_')
    #    cur_mouse_id = "_".join(vals[0:3])
    #    cur_session_id = vals[3]
    #    cur_fov_id = int(vals[4].replace("fov", ""))
    #    cur_stim_freq = int(vals[5].replace('stim', ''))
#
    #    mask = (df['stim_freq'] == cur_stim_freq) & \
    #            (df['mouse_id'] == cur_mouse_id) & \
    #            (df['session_id'] == cur_session_id) & \
    #            (df['fov_id'] == cur_fov_id)
#
    #    # Determine if neurons have a higher or lower entrainment 8Hz during stimulation
    #    # when compared to flicker stimulation during baseline
    #    
    #    if dbs_flicker_point > flicker_point:            
    #        df.loc[mask, 'Power_Delta'] = 1
    #        
    #    elif dbs_flicker_point < flicker_point:
    #        df.loc[mask, 'Power_Delta'] = 0

    data = np.vstack((base_pre_power['power'].values, flicker_pre_power['power'].values,
             stim_power['power'].values, flicker_post_power['power'].values))
    labels = ['Base', 'Flicker Onset', 'ESTIM', 'Flicker Post']


# %% Regression of 8 Hz power change between flicker onset and stim to
# stim frequency change

# Calculate the 8 Hz power change for each neuron
mask_8hz_spec = (all_spectra_df['freq_spec'] > 7.5) & (all_spectra_df['freq_spec'] < 8.5)
mask_flicker_pre = (all_spectra_df['timeline'] >= 0) & (all_spectra_df['timeline'] < 1)
mask_stim = (all_spectra_df['timeline'] >= 1) & (all_spectra_df['timeline'] < 2)

flicker_pre_8hz_df = all_spectra_df[mask_8hz_spec & mask_flicker_pre].groupby(['neuron', 'stim_freq'])['power'].mean()
stim_8hz_df = all_spectra_df[mask_8hz_spec & mask_stim].groupby(['neuron', 'stim_freq'])['power'].mean()

pow_diff_df = pd.DataFrame()
pow_diff_df['8Hz_pow_diff'] = stim_8hz_df - flicker_pre_8hz_df 

# Calculate the power difference at stimulation frequency
for stim_freq in all_spectra_df['stim_freq'].unique():
    mask_stim_freq = (all_spectra_df['freq_spec'] > 0.95*stim_freq) &\
        (all_spectra_df['freq_spec'] < 1.05*stim_freq)
    
    flicker_pre_stim_df = all_spectra_df[mask_stim_freq &\
        mask_flicker_pre].groupby(['neuron', 'stim_freq'])['power'].mean()
    
    stim_pow_df = all_spectra_df[mask_stim_freq &\
        mask_stim].groupby(['neuron', 'stim_freq'])['power'].mean()

    pow_diff_df['stim_pow_diff'] = stim_pow_df - flicker_pre_stim_df

pow_diff_df = pow_diff_df.reset_index(level='stim_freq')

# Plot the regression plot
plt.figure(figsize=(10, 10))
plt.rcParams['font.size'] = 20

#TODO perform the regression for each group
for stim_freq in all_spectra_df['stim_freq'].unique():
    stim_df = pow_diff_df[pow_diff_df['stim_freq'] == stim_freq]
    line1, = plt.plot(stim_df['8Hz_pow_diff'].values,\
             stim_df['stim_pow_diff'].values, '.' , label=str(stim_freq),\
             markersize=10)
    
    # Perform regression and plot line
    X = sm.add_constant(stim_df['8Hz_pow_diff'].values)

    model = sm.OLS(stim_df['stim_pow_diff'].values, X)
    results = model.fit()
    y_pred = results.predict(X)

    prev_color = line1.get_color()    
    # Plot the regression line
    plt.plot(stim_df['8Hz_pow_diff'].values, y_pred,\
             label='p= ' + str(round(results.pvalues[1], 2))\
                 + ' r= ' + str(round(results.params[1], 2)),\
                color=prev_color)

plt.xlabel('8Hz Power Difference')
plt.ylabel('Stimulation Frequency Power Difference')

plt.legend(loc='best')
plt.show()
