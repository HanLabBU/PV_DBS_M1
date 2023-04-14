clear all;
close all;
clc;

f = filesep;
local_rootpath = '~/Projects/';

data_path = [local_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'Stim Recordings' f];
[files, path] = uigetfile([data_path '*.mat'], 'Multiselect', 'on');
savepath = path;

% Use fixed framerate
Fs = 1/(1.2*10^-3);

% Read in traces
data = load([path files]);

trial_nums = unique(data.result.trial_vec);

all_subVm = [];
all_rasters = [];
for trial=trial_nums
    trial
    cur_trace = data.result.traces(data.result.trial_vec == trial);    
    spike_info = spike_detect_SNR_sim3(cur_trace, 3.75, 4, 7);

    % Store the subthreshold Vm
    [baseline, coeff] = Multi_func.exp_fit_Fx(spike_info.trace_ws', round(Fs));
    %cur_subVm = spike_info.trace_ws - baseline;
    cur_subVm = cur_trace';
    all_subVm = horzcat_pad(all_subVm, cur_subVm');

    % Store the spike raster
    cur_raster = spike_info.roaster;
    all_rasters = horzcat_pad(all_rasters, cur_raster');
end

% Show the average sub Threshold
figure;
timeline = (1:size(all_subVm, 1))./Fs - 1;
plot(timeline, mean(all_subVm, 2, 'omitnan'));
title('Plotting average subthreshold Vm');

% Plot the spike rasters
figure;
timeline = (1:size(all_subVm, 1))./Fs - 1;
for i=1:size(all_rasters, 2)
    cur_color = [rand, rand, rand];
    plot(timeline(all_rasters(:, i) == 1), i*5, '.', 'Color', cur_color);
    hold on;
end
title('Plotting spike raster');

%TODO Show the subthreshold spectra

