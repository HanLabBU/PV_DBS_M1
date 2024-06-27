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
          +f+ 'Flicker' +f+ 'Spectra' +f
df =pd.read_pickle(interm_data_path)
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
samp_freq = 500
freq_limit = [1, 60]
freq_nums = 3*(freq_limit[1] - freq_limit[0])
wavelet = 'morl'

test_freqs = df['stim_freq'].unique()
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Store all of the neurons trial-averaged vm
    vm_df = pd.DataFrame()

    power_spec_df = pd.DataFrame()

    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue
        flicker_start = fov_df[fov_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials normalized spike amplitude and set it as a new column
        # in the FOV dataframe
        nr_power_spec_df = pd.DataFrame()
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            norm_interp_subvm = trial_df['interp_subvm'].values/avg_sp_amp
            fov_df.loc[fov_df['trial_id'] == tr_val, 'interp_norm_vm'] = norm_interp_subvm

            # Calculate the trial power spectra
            calc_freqs, coeffs = fcwt.cwt(norm_interp_subvm, samp_freq, freq_limit[0], freq_limit[1], freq_nums)
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


        #TODO change this for calculating each neuron's average and then the final population average
        avg_tr_power_df = nr_power_spec_df.groupby(['timeline', 'freq_spec'])['power'].mean()
        avg_tr_power_df = avg_tr_power_df.unstack(level='timeline')
        avg_tr_power_df.sort_index(ascending=False, inplace=True)
        nr_norm_vals_df = avg_tr_power_df

        #--Testing different normalizations
        #Z-score across frequencies 

        #nr_norm_vals_df = nr_norm_vals_df.apply(lambda col: zscore(col, ddof=0), axis=0)

        #Z-score across time
        #nr_norm_vals_df = nr_norm_vals_df.apply(lambda row: zscore(row, ddof=0), axis=1)
        
        #TODO lmao good luck doing this in a data frame
        # (x - B)/(A + B) normalization where 
        # B is the 1 sec flicker period
        # A is the 1 sec electrical stimulation
        #print(power_vals.shape)
        #flicker_base = power_vals[:, np.where((timeline > 0) & (timeline <=1))[0]]
        #stim_period = power_vals[:, np.where((timeline > 1) & (timeline <=2))[0]]
        #avg_flicker_base = np.mean(flicker_base, axis=1).reshape((-1, 1))
        #avg_stim_period = np.mean(stim_period, axis=1).reshape((-1, 1))
        #norm_vals = (power_vals - avg_flicker_base)/(avg_flicker_base + avg_stim_period)

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
    
    # Print out the number of neurons for each frequency
    print(stim_freq)
    print('Num neurons ' + str(roi_acum + 1))

    # Get the 2D average of the heatmap
    avg_power_spec = power_spec_df.groupby(['timeline', 'freq_spec'])['power'].mean()
    #TODO need to figure out how to pivot back with the indices in the above Series
    heatmap_df = avg_power_spec.unstack(level='timeline')
            
    cm = 1/2.54
    plt.figure(figsize=(10.5 * 10* cm, 7 * 10* cm))
    surf_p = plt.pcolormesh(timeline, heatmap_df.index.values, heatmap_df.values) #, vmin=-1, vmax=1
        
    # Plot the stimultion protocol
    plt.plot(timeline, 2*trial_df['flicker_raster'] + 13 + np.max(heatmap_df.index.values), '-b')
    plt.plot([timeline[0], timeline[-1]], np.ones((2))* (10 + np.max(heatmap_df.index.values)), color=consts.pulse_color)
    stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
    plt.plot(timeline[stim_idx], np.ones_like(stim_idx)*(10 + np.max(heatmap_df.index.values)), '|',\
         linewidth=20,  color=consts.pulse_color)
    
    # Remove the spines
    ax = plt.gca()
    ax.spines['top'].set_visible(False)

    #Labels
    plt.xlabel('Time from flicker onset (s)')
    plt.ylabel('Frequency (Hz)')
    plt.title(stim_freq)

    # Add colorbar axis
    fig = plt.gcf()
    cbar_ax = fig.add_axes((0.92, 0.35, 0.02, 0.35))
    cbar = fig.colorbar(surf_p, label='Power', cax=cbar_ax)

    save_filename = savefig_path + 'Population_Flicker_Stim_Lesion_' +\
            str(stim_df['stim_freq'].unique()[0]) + 'Hz'
    
    plt.savefig(save_filename + '.svg', format='svg')
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

    #TODO need to change this to have values that are neuron wise for individual points and average across the periods
    base_period_df = power_spec_df[(power_spec_df['timeline'] < 0) | (power_spec_df['timeline'] > 3)]
    flicker_only_period_df = power_spec_df[((power_spec_df['timeline'] >= 0) & (power_spec_df['timeline'] < 1) ) |\
                                      ((power_spec_df['timeline'] >= 2) & (power_spec_df['timeline'] < 3))]
    dbs_flicker_period_df = power_spec_df[(power_spec_df['timeline'] >= 1) & (power_spec_df['timeline'] < 2)]

    # Get the stim frequency power bands
    # TODO need to look at each of the power bands
    base_period_data = base_period_df[(base_period_df['freq_spec'] > 7.5) &\
                                      (base_period_df['freq_spec'] < 8.5) ]
    flicker_only_period_data = flicker_only_period_df[(flicker_only_period_df['freq_spec'] > 7.5) &\
                                                      (flicker_only_period_df['freq_spec'] < 8.5)]
    dbs_flicker_data = dbs_flicker_period_df[(dbs_flicker_period_df['freq_spec'] > 7.5) &\
                                             (dbs_flicker_period_df['freq_spec'] < 8.5)]

    # Take the average for each neuron across time
    base_period_data = base_period_data.groupby(['neuron']).mean()
    flicker_only_period_data = flicker_only_period_data.groupby(['neuron']).mean()
    dbs_flicker_data = dbs_flicker_data.groupby(['neuron']).mean()

    data = [base_period_data['power'].values, flicker_only_period_data['power'].values, dbs_flicker_data['power'].values]
    labels = ['Base', 'Flicker', 'ESTIM+Flicker']

    # Perform statistical tests to determine if the 8 Hz is actually reduced
    statistic, p_val = stats.kruskal(data[0], data[1], data[2])

    print('Kruskal Wallis Test: ' + str(stim_freq) + 'Hz')
    print('H stat: ' + str(statistic))
    print('P-value: ' + str(p_val))
    
    # Perform post hoc comaparisons test
    p_values = posthocs.posthoc_dunn(data)
    print(p_values)
    print('1-Base, 2-flicker, 3-flicker+ESTIM')

    # Compare flicker period and flicker+DBS
    statistic, p_val = stats.mannwhitneyu(data[1], data[2])
    print('Flicker vs. DBS Comparison Test: ' + str(stim_freq) + 'Hz')
    print('H stat: ' + str(statistic))
    print('P-value: ' + str(p_val))

    #-- Violin plotting using seaborn --
    
    #(Deprecated) -- using matplotlib for violins
    #plt.figure()
    #plt.rcParams.update({'font.size': 12})
    #sns.violinplot(data, inner='box', alpha=0.5)
#
    ## Fix the cosmetics of the plot
    #ax = plt.gca()
    #ax.spines['top'].set_visible(False)
    #ax.spines['right'].set_visible(False)
    #ax.set_xticklabels(labels)
    #ax.set_ylabel('8 Hz Vm Power (A.U.)')
    #
    ##plt.legend([str(p_val)])
    #plt.title(str(stim_freq) + ' Hz ESTIM')
#
    #save_filename = savefig_path + '8Hz_Vm_Power_violin_' +\
    #        str(stim_df['stim_freq'].unique()[0]) + 'Hz'
    #
    #plt.savefig(save_filename + '.svg', format='svg')
    #plt.savefig(save_filename + '.png', format='png')
    
    # Plotting using matlab's engine and violinplot functions
    eng.feval('figure', 'Renderer', 'Painters', nargout=0)
    eng.eval(get_violin_opts(), nargout=0)

    data = base_period_data['power'].values
    data = np.concatenate((data, flicker_only_period_data['power'].values))
    data = np.concatenate((data, dbs_flicker_data['power'].values))

    eng.workspace['data'] = data

    # Concatenate columnwise to plot between points
    data_col = base_period_data['power'].values.reshape(-1, 1)
    data_col = np.concatenate((data_col, flicker_only_period_data['power'].values.reshape(-1, 1)), axis=1)
    data_col = np.concatenate((data_col, dbs_flicker_data['power'].values.reshape(-1, 1)), axis=1)

    eng.workspace['data_col'] = data_col

    labels = ['Base' for _ in range(len(base_period_data['power'].values))]
    labels.extend(['Flicker' for _ in range(len(flicker_only_period_data['power'].values))])
    labels.extend(['ESTIM+Flicker' for _ in range(len(dbs_flicker_data['power'].values))])
    eng.workspace['labels'] = labels

    # Plotting the violins in matlab
    
    # Store matlab commands in this variable
    matlab_exp = ""

    # Plotting the violin
    matlab_exp += "v = violinplot(data, labels, 'GroupOrder',"
    matlab_exp += "{'Base', 'Flicker', 'ESTIM+Flicker'}, opts);\n"
    matlab_exp += "v(1).ViolinColor = {Multi_func.base_color};\n"
    matlab_exp += "v(2).ViolinColor = {Multi_func.flick_color};\n"
    matlab_exp += "v(3).ViolinColor = {Multi_func.stim_color};\n"
    
    # Clean the axis
    matlab_exp += "Multi_func.set_default_axis(gca);\n"
    
    # Plot the lines between each points of the violin plot
    matlab_exp += "hold on;\n"
    matlab_exp += "plot([1, 2, 3], data_col, '-k');\n"

    # Adding stim frequency to title
    eng.workspace['freq'] = stim_freq
    matlab_exp += "title([num2str(freq) ' Hz' ], 'Interpreter', 'none');\n"
    eng.workspace['savefig_path'] = savefig_path
    matlab_exp += "saveas(gcf, [savefig_path '8Hz_Vm_Power_violin_' num2str(freq) 'Hz.png']);\n"
    matlab_exp += "saveas(gcf, [savefig_path '8Hz_Vm_Power_violin_' num2str(freq) 'Hz.svg']);\n"

    eng.eval(matlab_exp, nargout=0)
    
    # Plot each neuron's violin vales
    #for nr_name, row in base_period_data.iterrows():
    #    base_point = base_period_data.loc[nr_name, 'power']
    #    flicker_point = flicker_only_period_data.loc[nr_name, 'power']
    #    dbs_flicker_point = dbs_flicker_data.loc[nr_name, 'power']
#
    #    # Create a figure
    #    plt.figure()
    #    plt.plot(np.array([1, 2, 3]), np.array([base_point, flicker_point, dbs_flicker_point]))
    #    plt.title(nr_name)

    plt.show()


# %%
plt.close('all')
eng.eval("close all", nargout=0)

# %%
# Plot the power spectra for each neuron
samp_freq = 500
freq_limit = [1, 60]
freq_nums = 3*(freq_limit[1] - freq_limit[0])
wavelet = 'morl'

test_freqs = df['stim_freq'].unique()
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        # Skip if dataframe is empty()
        if fov_df.empty:
            continue
        flicker_start = fov_df[fov_df['flicker_raster'] == 1]['interp_time'].values[0]
        # Timeline for all
        timeline = fov_df[fov_df['trial_id'] == fov_df['trial_id'].unique()[0]]['interp_time'].values - flicker_start
        
        #Calculate each trials normalized spike amplitude and set it as a new column
        # in the FOV dataframe
        tr_power_spec_df = pd.DataFrame()
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)
            norm_interp_subvm = trial_df['interp_subvm'].values/avg_sp_amp
            fov_df.loc[fov_df['trial_id'] == tr_val, 'interp_norm_vm'] = norm_interp_subvm

            # Calculate the trial power spectra
            calc_freqs, coeffs = fcwt.cwt(norm_interp_subvm, samp_freq, freq_limit[0], freq_limit[1], freq_nums)
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

            tr_power_spec_df = pd.concat([tr_power_spec_df, pd.DataFrame(tr_dict)], ignore_index=True, join='outer')
            tr_power_spec_df['trial'] = tr_power_spec_df['trial'].astype('category')

        avg_tr_power_df = tr_power_spec_df.groupby(['timeline', 'freq_spec'])['power'].mean()
        nr_2d_power_df = avg_tr_power_df.unstack(level='timeline')
        nr_2d_power_df.sort_index(ascending=False, inplace=True)
        
        # Calculate the trial-averaged sub Vm 
        avg_trial_vm = fov_df.groupby('interp_time')['interp_norm_vm'].mean()

        # Setup figure
        cm = 1/2.54
        fig, axs = plt.subplots(3, 1, figsize=(3*10.5*cm, 3*10*cm))
        timeline

        # Plot the neuron power spectra
        axs[0].pcolormesh(timeline, calc_freqs, nr_2d_power_df.values)
        
        # Plot the stimultion protocol
        axs[0].plot(timeline, 2*trial_df['flicker_raster'] + 13 + np.max(calc_freqs), '-b')
        axs[0].plot([timeline[0], timeline[-1]], np.ones((2))* (10 + np.max(calc_freqs)), color=consts.pulse_color)
        stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
        axs[0].plot(timeline[stim_idx], np.ones_like(stim_idx)*(10 + np.max(calc_freqs)), '|',\
             linewidth=20,  color=consts.pulse_color)

        # Plot the neuron sub Vm
        axs[1].plot(timeline, avg_trial_vm.values)
        axs[1].plot(timeline, 2*trial_df['flicker_raster'] + 1 + avg_trial_vm.max(), '-b')
        axs[1].plot([timeline[0], timeline[-1]], np.ones((2))* (1 + avg_trial_vm.max()), color=consts.pulse_color)
        stim_idx = np.where(trial_df['stim_raster'].values == 1)[0]
        axs[1].plot(timeline[stim_idx], np.ones_like(stim_idx)*(1 + avg_trial_vm.max()), '|',\
             linewidth=20,  color=consts.pulse_color)

        # Plot the line for each period
        #Old and deprecated way of doing things
        #hz8_idx = (avg_tr_power_df.index < 8.5) & (avg_tr_power_df.index > 7.5)
        #base_period_data = avg_tr_power_df.loc[hz8_idx, (avg_tr_power_df.columns.values < 0).astype(bool) | (avg_tr_power_df.columns.values > 3).astype(bool)]
        #flicker_only_period_data = avg_tr_power_df.loc[hz8_idx, \
        #    ((avg_tr_power_df.columns.values >=0).astype(bool) & (avg_tr_power_df.columns.values < 1).astype(bool) )|\
        #    ((avg_tr_power_df.columns.values >= 2).astype(bool) & (avg_tr_power_df.columns.values < 3).astype(bool))]
        #dbs_flicker_data= avg_tr_power_df.loc[hz8_idx, (avg_tr_power_df.columns.values >= 1).astype(bool) | (avg_tr_power_df.columns.values < 2).astype(bool)]
        #
        #spec_8hz_data = np.array([np.mean(base_period_data), np.mean(flicker_only_period_data),\
        #                          np.mean(dbs_flicker_data)])
        
        # Filter data frame from column values
        avg_tr_power_df = avg_tr_power_df.reset_index()
        base_period_data = avg_tr_power_df[\
            ((avg_tr_power_df['timeline'] < 0) | (avg_tr_power_df['timeline'] > 3)) &\
            ((avg_tr_power_df['freq_spec'] > 7.5) & (avg_tr_power_df['freq_spec'] < 8.5)) ]

        flicker_only_period_data = avg_tr_power_df[\
            (((avg_tr_power_df['timeline'] >= 0) & (avg_tr_power_df['timeline'] < 1)) | \
            ((avg_tr_power_df['timeline'] >= 2) & (avg_tr_power_df['timeline'] < 3)))
            & ((avg_tr_power_df['freq_spec'] > 7.5) & (avg_tr_power_df['freq_spec'] < 8.5)) ]
        
        dbs_flicker_period_data = avg_tr_power_df[\
            ((avg_tr_power_df['timeline'] >= 1) & (avg_tr_power_df['timeline'] < 2)) &\
            ((avg_tr_power_df['freq_spec'] > 7.5) & (avg_tr_power_df['freq_spec'] < 8.5)) ]
        
        # Consolidate all of the 8 Hz data
        spec_8hz_data = np.array([base_period_data['power'].mean(), \
                    flicker_only_period_data['power'].mean(), dbs_flicker_period_data['power'].mean()])

        axs[2].plot(np.array([1, 2, 3]), spec_8hz_data, '-k')
        axs[2].set_xticks([1, 2, 3])
        axs[2].set_xticklabels(['base', 'flicker', 'dbs+flicker'])

        # Add title to identify neuron
        fig.suptitle("".join(str(values)) + str(stim_freq) + ' Hz')
        # Save the individual plots
        
        # Set the individual figure name
        save_filename = savefig_path + 'Individual' +f + str(values[0]) +'_'+ str(values[1]) \
                            +'_fov'+ str(values[2]) +'_'+\
                            str(stim_freq) + 'Hz_Whole_Period_Average'

        # Save figure
        plt.savefig(save_filename + '.svg', format='svg')
        plt.savefig(save_filename + '.png', format='png')

# %%
plt.close('all')
eng.eval("close all", nargout=0)

# %%