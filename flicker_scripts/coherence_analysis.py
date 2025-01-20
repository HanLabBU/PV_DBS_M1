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
import decimal
from glob import glob
from pymatreader import read_mat
from scipy import signal

from scipy import stats
from scipy.stats import zscore
from scipy.ndimage import uniform_filter1d
import fcwt
from pycwt import wavelet as wv
import seaborn as sns
import scipy.stats as stats
import scikit_posthocs as posthocs
import itertools

import matlab.engine
# Test out matlab code
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
          +f+ 'Flicker' +f+ 'Coherence' +f
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
# Calculate the wavelet coherence by using pycwt
samp_freq = 500
freq_limit = [1, 40]
freq_nums = 3*(freq_limit[1] - freq_limit[0])

# Reset stats file
open(savefig_path + 'Individual' +f+ 'single_cell_coh_stats.txt', 'w').close()

# Cosmetic parameters
plt.rcParams['pdf.fonttype'] = 42

#  Loop through each frequency
test_freqs = df['stim_freq'].unique() # [40] #   
for stim_freq in test_freqs:
    stim_df = df[df['stim_freq'] == stim_freq]
    
    mouse_list = stim_df['mouse_id'].unique()
    session_list = stim_df['session_id'].unique()
    fov_list = stim_df['fov_id'].unique()

    pairings = list(itertools.product(mouse_list, session_list, fov_list))
    
    # Setup dataframe for population coherence
    pop_coh_df = pd.DataFrame()

    # Loop through unique FOVs
    roi_acum = 0
    for values in pairings:
        fov_df = stim_df[(stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])]
        
        cur_nr_mask = (df['stim_freq'] == stim_freq) & \
                (df['mouse_id'] == values[0]) & \
                (df['session_id'] == values[1]) & \
                (df['fov_id'] == values[2])
        
        cur_nr_id = "_".join([str(e) for e in values])

        # Skip if dataframe is empty()
        if fov_df.empty:
            continue

        # Calculate the stim protocol signal

        nr_coh_df = pd.DataFrame()
        # Calculate spike-amplitude normalized Vm
        i = 0
        for tr_val in fov_df['trial_id'].unique():
            trial_df = fov_df[fov_df['trial_id'] == tr_val]
            avg_sp_amp = np.nanmean(trial_df['spike_amp_raster'].values)

            # If no spikes were detected, just use a value of 1
            if np.isnan(avg_sp_amp):
                avg_sp_amp = 1

            norm_vm = trial_df['detrend_trace'].values/avg_sp_amp # or 'sub_vm' Switching between the full Vm and the Sub Vm
            fov_df.loc[fov_df['trial_id'] == tr_val, 'norm_vm'] = norm_vm

            # Flicker protocol signal
            flicker_sig = trial_df['flicker_raster'].values

            # Manual method -- Get the wavelet transform raw values
            #tr_calc_freqs, tr_coeffs = fcwt.cwt(norm_vm, samp_freq, freq_limit[0], freq_limit[1], freq_nums)
            #fl_calc_freqs, fl_coeffs = fcwt.cwt(flicker_sig, samp_freq, freq_limit[0], freq_limit[1], freq_nums)
            #w_coh = consts.wav_coherence(tr_coeffs, fl_coeffs)

            #TODO need to do more reading on this package's documentation
            # Infor about the wavelet parameters
            # s0: small scale refers to higher frequency sampling
            # J: for max scale to s0 * 2**(J * dj)
            # dj: spacing between scales
            wct_res = wv.wct(norm_vm, flicker_sig, 1/samp_freq, s0=1/freq_limit[1], dj=1/30, sig=False)
            freqs = wct_res[3] # The frequencies used
            wav_coh = wct_res[0] # The coherences across time

            # Linear the power
            linear_coh = wav_coh.ravel()
            t_mesh, freq_mesh = np.meshgrid(timeline, freqs)
            t_lin = t_mesh.ravel()
            freq_lin = freq_mesh.ravel()

            tr_dict = {
                'trial':tr_val,
                'power':linear_coh,
                'timeline':t_lin,
                'freq_spec':freq_lin
            }

            nr_coh_df = pd.concat([nr_coh_df, pd.DataFrame(tr_dict)], ignore_index=True, join='outer')
            nr_coh_df['trial'] = nr_coh_df['trial'].astype('category')
                
            # Calculate the timeline
            timeline = fov_df[fov_df['trial_id'] == tr_val]['interp_time'].values
            
            ##--Plot the individual trial coherence with flicker
            #cm = 1/2.54
            #fig, axs = plt.subplots(3, 1, figsize=(3*10.5*cm, 3*10*cm))
            ##Plot the stim protocol
            #axs[0].plot(timeline, np.max(fl_calc_freqs) + flicker_sig, '-b')
            ## Plot the coherence spectrum
            #surf_p = axs[0].pcolormesh(timeline, freqs, wav_coh, cmap='jet')
            #axs[0].set_title('Wavelet Coherence pycwt')
            #surf_p = axs[1].pcolormesh(timeline, tr_calc_freqs, np.abs(tr_coeffs), cmap='jet')
            #axs[1].set_title('Trial CWT')
            #surf_p = axs[2].pcolormesh(timeline, fl_calc_freqs, np.abs(fl_coeffs), cmap='jet')
            #axs[2].set_title('Flicker CWT')

            
            ##--Plot the MATLAB wavelet
            #eng.workspace['norm_vm'] = norm_vm
            #eng.workspace['flicker_sig'] = flicker_sig
            #eng.workspace['fs'] = samp_freq
            #eng.workspace['timeline'] = timeline
            #eng.workspace['freq_limits'] = freq_limit
            #            
            ## Calculate wavelet transform in MATALB
            #matlab_exp = ""
            ## Need to convert the frequency limits
            #matlab_exp += "freq_limits = cell2mat(freq_limits);\n"
            #matlab_exp += "[wcoh, wcs, f] = wcoherence(norm_vm, flicker_sig, fs,"
            #matlab_exp +=       " FrequencyLimits=freq_limits);\n"
#
            ## Start the figure
            ##eng.feval('figure', 'Renderer', 'Painters', nargout=0)
#
            ## Plot the wavelet coherence
            ##matlab_exp += "disp(f)\n"
            #matlab_exp += "surface(timeline, f, wcoh, "
            #matlab_exp += "'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none')\n"
#
            ##eng.eval(matlab_exp, nargout=0)
            
            # Wait for a few iterations
            i += 1
            if i == 2:
                continue
                raise Exception('One Trial')

        # Calculate trial-averaged power
        avg_tr_coh_df = nr_coh_df.groupby(['timeline', 'freq_spec'])['power'].mean()
        avg_tr_coh_df.name = 'power' # Specify the mean value to power before unstacking
        avg_tr_coh_df = avg_tr_coh_df.unstack(level='timeline')
        avg_tr_coh_df.sort_index(ascending=False, inplace=True)

        # Setup dataframe for concatenation to neuron population
        avg_tr_coh_df = avg_tr_coh_df.stack()
        avg_tr_coh_df.name = 'power' # Specify the mean value to power before unstacking
        avg_tr_coh_df = avg_tr_coh_df.reset_index()

        avg_tr_coh_df['neuron'] = cur_nr_id
        avg_tr_coh_df['neuron'] = avg_tr_coh_df['neuron'].astype('category')

        ##--Plot neuron average coherence
        #TODO plot the individual spectral power with the coherences as well
        #cm = 1/2.54
        #plt.figure(figsize=(24 * cm, 12 * cm))
        #
        ## Setup neuron coherence heatmap
        #nr_heatmap_df = avg_tr_coh_df.copy()
        #nr_heatmap_df = nr_heatmap_df.pivot(index='freq_spec', columns='timeline', values='power')
        #
        #surf_p = plt.pcolormesh(timeline, nr_heatmap_df.index.values, nr_heatmap_df.values, cmap='jet')
        #
        ## Remove the spines
        #ax = plt.gca()
        #ax.spines['top'].set_visible(False)
        #
        ##Labels
        #plt.xlabel('Time from flicker onset (s)')
        #plt.ylabel('Frequency (Hz)')
        #plt.title('Coherence ' + avg_tr_coh_df['neuron'].unique()[0] +"_"+ str(stim_freq))
        #
        ## Save neuron coherence
        #save_filename = savefig_path + 'Individual' +f+ 'Coh_'+ avg_tr_coh_df['neuron'].unique()[0] +"_"+ str(stim_freq)
        #plt.savefig(save_filename +'.png', format='png')

        # Append neuron average to population dataframe
        pop_coh_df = pd.concat([pop_coh_df, avg_tr_coh_df], ignore_index=True, join='outer')
        pop_coh_df['neuron'] = pop_coh_df['neuron'].astype('category')

        #--Perform multi-comparison of the power between the stim periods
        mask_base_pre = (nr_coh_df['timeline'] < 0)
        mask_flicker_pre = (nr_coh_df['timeline'] >= 0) & (nr_coh_df['timeline'] < 1)
        mask_stim = (nr_coh_df['timeline'] >= 1) & (nr_coh_df['timeline'] < 2)
        mask_flicker_post = (nr_coh_df['timeline'] >= 2) & (nr_coh_df['timeline'] < 3)

        mask_8Hz_spec = (nr_coh_df['freq_spec'] > 7.5) & (nr_coh_df['freq_spec'] < 8.5)

        nr_base_pre_power = nr_coh_df[mask_base_pre & mask_8Hz_spec].groupby(['trial']).mean()
        nr_flicker_pre_power = nr_coh_df[mask_flicker_pre & mask_8Hz_spec].groupby(['trial']).mean()
        nr_stim_power = nr_coh_df[mask_stim & mask_8Hz_spec].groupby(['trial']).mean()
        nr_flicker_post_power = nr_coh_df[mask_flicker_post & mask_8Hz_spec].groupby(['trial']).mean()

        #DEBUG trying to find the nan values
        if np.any(np.isnan(nr_base_pre_power)):
            print(nr_base_pre_power)
            print(nr_coh_df[nr_coh_df['trial'] == 6])

            raise Exception('Found the NaN neuron')

        # Compile into data structures for comparisons
        data = np.vstack((nr_base_pre_power['power'].values, nr_flicker_pre_power['power'].values,
                          nr_stim_power['power'].values, nr_flicker_post_power['power'].values))
        labels = ['Base', 'Flicker_Onset', 'ESTIM', 'Flicker Post']

        # Write to large stats data
        #TODO maybe also save all of the powers for each period in arrays to be plotted later?
        with open(savefig_path + 'Individual' +f+ 'single_cell_coh_stats.txt', 'a') as stats_file:
            # Stop executions if NaNs are found
            if np.isnan(statistic):
                raise Exception('Found nan')

            stats_file.write("------"+ cur_nr_id +'_'+ str(stim_freq) +"------\n")
            # Wilcoxon sign-rank for flicker_onset vs stim periods
            statistic, p_val = stats.wilcoxon(data[1], data[2])
            stats_file.write('Sign-Rank Test\n')
            stats_file.write('Statistic: ' + str(statistic) + '\n')
            stats_file.write('P-value: ' + str(p_val) + '\n')
            
            # Determine the modulation of the stimulation to flicker coherence
            if p_val < 0.05 and np.mean(data[2] - data[1]) > 0:
                mod = 1
            elif p_val < 0.05 and np.mean(data[2] - data[1]) < 0:
                mod = -1
            else:
                mod = 0

            # TODO the stim_df may not be necessary
            # Set the mod factor in the stimulation dataframe
            cur_nr_mask = (stim_df['mouse_id'] == values[0]) & \
                         (stim_df['session_id'] == values[1]) & \
                            (stim_df['fov_id'] == values[2])

            # Set the flicker modulation
            stim_df.loc[cur_nr_mask, 'stim_mod'] = mod

            # Set stim_mod to the coherence dataframe
            pop_coh_df.loc[pop_coh_df['neuron'] == cur_nr_id, 'stim_mod'] = mod

            # Kruskal wallis for all periods
            statistic, p_val = stats.kruskal(data[0], data[1], data[2], data[3])
            stats_file.write('KW for periods\n')
            stats_file.write('H stat: %s\n' % statistic)
            stats_file.write('P-value: %s\n' % p_val)

            # Perform post hoc comaparisons test
            p_values = posthocs.posthoc_dunn(data)
            stats_file.write('Post-hoc Dunn Test\n')
            stats_file.write('|  Base    | flick_s |  ESTIM |flick_p |\n')
            stats_file.write(str(p_values))

            stats_file.write('\n\n')

        #DEBUG
        #raise Exception('Save neuron plot')

        
    # Calculate the population average coherence and create 2D heatmap matrix
    avg_pop_coh_df = pop_coh_df.groupby(['timeline', 'freq_spec'])['power'].mean()
    avg_imprv_coh_df = pop_coh_df[pop_coh_df['stim_mod'] == 1].groupby(['timeline', 'freq_spec'])['power'].mean()
    avg_wors_coh_df = pop_coh_df[pop_coh_df['stim_mod'] == -1].groupby(['timeline', 'freq_spec'])['power'].mean()
    avg_unch_coh_df = pop_coh_df[pop_coh_df['stim_mod'] == 0].groupby(['timeline', 'freq_spec'])['power'].mean()

    heatmap_df = avg_pop_coh_df.unstack(level='timeline')

    #--Plot the 2D matrix heatmap, all neurons, and the individual populations
    cm = 1/2.54
    fig, axs = plt.subplots(4, 1, figsize=(20 * cm, 25 * cm))
    surf_p = axs[0].pcolormesh(timeline, heatmap_df.index.values, heatmap_df.values, cmap='jet')

    # Remove the spines
    axs[0].spines['top'].set_visible(False)

    #Labels
    axs[0].set_xlabel('Time from flicker onset (s)')
    axs[0].set_ylabel('Frequency (Hz)')
    axs[0].set_title('Coherence all neurons' + str(stim_freq))

    # Set the spacing for y-axis text
    axs[0].tick_params(axis='y', pad=0)
    num_nrs = pop_coh_df['neuron'].unique().shape[0]
    axs[0].legend(title='Num nrs: ' + str(num_nrs))
    # Set the position of the main figure in relation to the colorbar
    #axs[0].set_position((0.15, 0.1, 0.62, 0.82))

    # Check that there are improved coherence neurons
    if not avg_imprv_coh_df.empty:
        heatmap_df = avg_imprv_coh_df.unstack(level='timeline')
        surf_p = axs[1].pcolormesh(timeline, heatmap_df.index.values, heatmap_df.values, cmap='jet')
        # Remove the spines
        axs[1].spines['top'].set_visible(False)

        #Labels
        axs[1].set_xlabel('Time from flicker onset (s)')
        axs[1].set_ylabel('Frequency (Hz)')
        axs[1].set_title('Coherence of improved neurons ' + str(stim_freq))

        # Set the spacing for y-axis text
        axs[1].tick_params(axis='y', pad=0)
        num_imprv_nrs = pop_coh_df[pop_coh_df['stim_mod'] == 1]['neuron'].unique().shape[0]
        axs[1].legend(title='Num nrs: ' + str(num_imprv_nrs))

    # Set the position of the main figure in relation to the colorbar
    #axs[1].set_position((0.15, 0.1, 0.62, 0.82))

    #TODO may need to figure this one out Add colorbar axis
    #fig = plt.gcf()
    #cbar_ax = fig.add_axes((0.80, 0.35, 0.02, 0.35))
    #cbar = fig.colorbar(surf_p, label='Normalized Power Change', cax=cbar_ax)

    # Set the default colorbar scale similar to main figures
    #surf_p.set_clim(-0.2, 0.6)
    if not avg_wors_coh_df.empty:
        heatmap_df = avg_wors_coh_df.unstack(level='timeline')
        surf_p = axs[2].pcolormesh(timeline, heatmap_df.index.values, heatmap_df.values, cmap='jet')
        # Remove the spines
        axs[2].spines['top'].set_visible(False)

        #Labels
        axs[2].set_xlabel('Time from flicker onset (s)')
        axs[2].set_ylabel('Frequency (Hz)')
        axs[2].set_title('Coherence of worsened neurons ' + str(stim_freq))

        # Set the spacing for y-axis text
        axs[2].tick_params(axis='y', pad=0)
        num_wors_nrs = pop_coh_df[pop_coh_df['stim_mod'] == -1]['neuron'].unique().shape[0]
        axs[2].legend(title='Num nrs: ' + str(num_wors_nrs))

    if not avg_unch_coh_df.empty:
        heatmap_df = avg_unch_coh_df.unstack(level='timeline')
        surf_p = axs[3].pcolormesh(timeline, heatmap_df.index.values, heatmap_df.values, cmap='jet')
        # Remove the spines
        axs[3].spines['top'].set_visible(False)

        #Labels
        axs[3].set_xlabel('Time from flicker onset (s)')
        axs[3].set_ylabel('Frequency (Hz)')
        axs[3].set_title('Coherence of unchanged neurons ' + str(stim_freq))

        # Set the spacing for y-axis text
        axs[3].tick_params(axis='y', pad=0)
        num_unch_nrs = pop_coh_df[pop_coh_df['stim_mod'] == 0]['neuron'].unique().shape[0]
        axs[3].legend(title='Num nrs: ' + str(num_unch_nrs))

    fig.tight_layout()

    # Save neuron coherence
    save_filename = savefig_path + 'Coherence_' + str(stim_freq) + "Hz"
    plt.savefig(save_filename +'.png', format='png')

    # Consolidate the powers for each period and keep as a population
    mask_base_pre = (pop_coh_df['timeline'] < 0)
    mask_flicker_pre = (pop_coh_df['timeline'] >= 0) & (pop_coh_df['timeline'] < 1)
    mask_stim = (pop_coh_df['timeline'] >= 1) & (pop_coh_df['timeline'] < 2)
    mask_flicker_post = (pop_coh_df['timeline'] >= 2) & (pop_coh_df['timeline'] < 3)

    mask_8Hz_spec = (pop_coh_df['freq_spec'] > 7.5) & (pop_coh_df['freq_spec'] < 8.5)

    #pop_base_pre_power = pop_coh_df[mask_base_pre & mask_8Hz_spec].groupby(['neuron']).mean()
    pop_flicker_pre_power = pop_coh_df[mask_flicker_pre & mask_8Hz_spec].groupby(['neuron']).mean()
    pop_stim_power = pop_coh_df[mask_stim & mask_8Hz_spec].groupby(['neuron']).mean()
    pop_flicker_post_power = pop_coh_df[mask_flicker_post & mask_8Hz_spec].groupby(['neuron']).mean()

    # Setup variables for MATLAB violin plots

    # Plot each violin
    data = pop_flicker_pre_power['power'].values
    data = np.concatenate((data, pop_stim_power['power'].values))
    data = np.concatenate((data, pop_flicker_post_power['power'].values))

    # To plot lines across each violin
    data_col = pop_flicker_pre_power['power'].values.reshape(-1, 1)
    data_col = np.concatenate((data_col, pop_stim_power['power'].values.reshape(-1, 1)), axis=1)
    data_col = np.concatenate((data_col, pop_flicker_post_power['power'].values.reshape(-1, 1)), axis=1)
    
    labels = ['Flicker Onset' for _ in range(len(pop_flicker_pre_power['power'].values))]
    labels.extend(['STIM' for _ in range(len(pop_stim_power['power'].values))])
    labels.extend(['Flicker Post' for _ in range(len(pop_flicker_post_power['power'].values))])
    eng.workspace['labels'] = labels

    # Transfer the data variable
    #data = nr_base_pre_power['power'].values.reshape(-1, 1)
    #data = np.concatenate((data, nr_flicker_pre_power['power'].values.reshape(-1, 1)), axis=1)
    #data = np.concatenate((data, nr_stim_power['power'].values.reshape(-1, 1)), axis=1)

    # Plot the violins in MATLAB
    eng.feval('figure', 'Renderer', 'Painters', nargout=0)
    eng.eval(get_violin_opts(), nargout=0)

    eng.workspace['data'] = data
    eng.workspace['data_col'] = data_col
    eng.workspace['labels'] = labels

    matlab_exp = ""

    # Plotting the violin
    matlab_exp += "v = violinplot(data, labels, 'GroupOrder',"
    matlab_exp += "{'Flicker Onset', 'STIM', 'Flicker Post'}, opts);\n"
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
    matlab_exp += "title([num2str(freq) ' Hz Coherence'], 'Interpreter', 'none');\n"
    matlab_exp += "ylabel('Coherence (A.U.)', 'Interpreter', 'none');\n"
    matlab_exp += "set(gca, 'FontSize', 14);\n"
    eng.workspace['savefig_path'] = savefig_path
    matlab_exp += "saveas(gcf, [savefig_path 'Flicker_8Hz_coherence_violin_' num2str(freq) 'Hz.png']);\n"
    matlab_exp += "saveas(gcf, [savefig_path 'Flicker_8Hz_coherence_violin_' num2str(freq) 'Hz.pdf']);\n"

    eng.eval(matlab_exp, nargout=0)


#TODO might also need to compare to MATLAB's wavelet coherence method

# %%
plt.close('all')
eng.eval("close all", nargout=0)