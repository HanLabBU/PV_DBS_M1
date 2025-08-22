import os
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from matplotlib.colors import LinearSegmentedColormap

# File meant to store all static variables

f = os.path.sep
save_path = f + 'home' + f + 'pierfier' + f + 'Dropbox' + f + 'RKC-HanLab' + f + 'Pierre CaMK DBS Project Dropbox' + f + 'Materials' + f + 'Plots' + f 

# Gonna make a gray color
stim_color = np.array([20, 33, 61])/255

# Color for the pulse bars
pulse_color = np.array([170, 176, 97])/255

# Red purplse blue colormap
red_purple_blue_cmap = LinearSegmentedColormap.from_list("rpb_cmap",\
        [(0.0, tuple(val / 255 for val in (0, 119, 182))),
         (0.5, tuple(val / 255 for val in (230, 230, 250))),
         (1.0, tuple(val / 255 for val in (214, 40, 40)))\
            ])

# Find the starts of each of the burst periods
def find_start_idx(raster):
    gap_idx = np.array([])
    stim_idx = np.where(raster > 0)[0]
    diff_idx = np.diff(stim_idx)
    
    # Get the minimum sized gaps
    min_dist = np.min(diff_idx)

    gap_idx = np.append(gap_idx, stim_idx[0])

    # Find the large gaps
    for index, value in enumerate(diff_idx):
        # Needs to be higher than a tolerance value
        #print(np.abs(min_dist - value))
        if np.abs(min_dist - value) > 200: # Gaps that are larger than 0.25 seconds given a sampling of 10,000
            gap_idx = np.append(gap_idx, stim_idx[index + 1])

    return gap_idx

# Normalize signals across columns
def norm_signals(sig_mat):
    return (sig_mat - np.min(sig_mat, axis=0))/(np.max(sig_mat, axis=0) - np.min(sig_mat, axis=0))

# Calculate the wavelet coherence from two wavelet info
def wav_coherence(wt1, wt2):
    S_xy = np.conj(wt1) * wt2
    S_xx = np.abs(wt1)**2
    S_yy = np.abs(wt2)**2

    print(np.abs(S_xy)**2)

    #print((np.abs(S_xy)**2)[10][10])
    #print(S_yy[10][10])
    #print(S_xx[10][10])
    #print((S_xx*S_yy)[10][10])
    #print('\n')

    return np.abs(S_xy)**2 / (S_xx * S_yy)
    


def get_ephys_rise_indices(ephys_signal):
    # Check whether to flip the trigger signal
    upper_volt = np.max(ephys_signal)
    lower_volt = np.min(ephys_signal)

    volt_range = abs(upper_volt - lower_volt)

    last_points = np.mean(ephys_signal[4990:5500])
    diff_low = abs(last_points - lower_volt)
    diff_high = abs(last_points - upper_volt)

    if diff_low > diff_high:
        ephys_signal = -1*ephys_signal;

    # Find the rising signal
    first_der = np.concatenate((np.array([0]), np.diff(ephys_signal)))
    #print(first_der.shape)

    thres = 0.03*(np.max(first_der))

    #DEBUG
    #plt.figure()
    ##plt.plot(ephys_signal)
    #plt.plot(first_der, '-g')
    

    first_der[first_der < thres] = 0
    peak_idx, _ = signal.find_peaks(first_der)
    
    #plt.plot(peak_idx, first_der[peak_idx], 'r|')
    #plt.hlines(y=thres, xmin=0, xmax=6000000, colors='black', linewidth=2)
    #plt.title('Derivative')
    #plt.show()


    #print(peak_idx.shape)

    #DEBUG
    #plt.figure()
    ##plt.plot(ephys_signal)
    #plt.plot(peak_idx, first_der[peak_idx], 'r|')
    #plt.plot(first_der, '-g')
    #plt.title('Derivative')
    #plt.show()

    return peak_idx

# Filter signal using FIR filter
def fir_filt(data, low_freq, high_freq, fs, numtaps=400, mode='same'):
    h = signal.firwin(numtaps, [low_freq, high_freq], pass_zero=False, fs=fs)
    return signal.fftconvolve(data, h, mode=mode)

# Calculate PLVs from the dataframe while specifying the column labels
def event_PLV(fov_df, phase_lbl, event_lbl, exc_criteria, timeshift=0): # Is it worth using timeshift here
    indi_phase_vecs = np.exp(1j*fov_df[fov_df[event_lbl] == 1][phase_lbl].values)


    # Number of events
    num_events = fov_df[event_lbl].sum()
    PLV = (1/num_events) * (np.abs(np.sum(indi_phase_vecs)))
    PLV2 = (1/(num_events - 1))*((PLV**2)*num_events - 1)

    #DEBUG
    #print('Num events: ' + str(num_events))
    #print('PLV: ' + str(PLV))
    #print('PLV2: ' + str(PLV2))
    #print('size of df ' + str(fov_df.shape))
    #print('event label: ' + event_lbl)
    return PLV, PLV2, indi_phase_vecs