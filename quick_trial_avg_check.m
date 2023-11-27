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
    cur_subVm = spike_info.trace_ws - baseline;
    %cur_subVm = cur_trace';
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
savefig(gcf, [path files(1:end - 5) '_subVmAvg.fig']);

% Plot the spike rasters
figure;
timeline = (1:size(all_subVm, 1))./Fs - 1;
for i=1:size(all_rasters, 2)
    %cur_color = [rand, rand, rand];
    cur_color = [0 0 0];
    spike_idx = find(all_rasters(:, i) == 1);

    if length(spike_idx) == 0
        continue;
    end

    %plot(timeline(spike_idx), repmat(i*5, length(spike_idx), 1), '.', 'Color', cur_color);
    plot(timeline(spike_idx), i*5, '.', 'Color', cur_color);
    hold on;
end
title('Plotting spike raster');
savefig(gcf, [path files(1:end - 5) '_spike_raster.fig']);

% Plot all of the traces
figure;
surface(timeline, 1:size(all_subVm, 2), all_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
title('All Traces Heatmap');
savefig(gcf, [path files(1:end - 5) '_trials_heatmap.fig']);

%TODO Show the subthreshold spectra
Fs = 828;
[wt, f] = get_power_spec(mean(all_subVm, 2, 'omitnan'), Fs);
figure;
surface(timeline, ... 
        f, ...
        abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
title('Power spectra of average signal');
savefig(gcf, [path files(1:end - 5) '_SpecAvg.fig']);

% Calculate cwt for input signal and 
function [wt, f] = get_power_spec(signal, samp_freq)
    freqLimits = [0 150];
    fb = cwtfilterbank(SignalLength=length(signal),...
                       SamplingFrequency=samp_freq,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(signal, FilterBank=fb);
end
