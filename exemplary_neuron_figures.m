% Housekeeping clear all stuff
clc;
clear all;
close all;
f = filesep;

% Other sample traces to try
%617100_M1_rec20211111_FOV3_140 trial #4 (has delta)
% 31556noeartag_M1_rec20221206_FOV1_140 (hyperpolarization) trial 3 or 9

% Specify which photobleaching detrending to use
% (2) for photobleaching only the baseline
% (1) for the photobleaching taking into account the stimulation bump
% (0) for the basic full trace expontential fit subtraction detrending

sophis_bleachdetrend = 1;

% Option for excluding or including first 200ms
exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;

% Maingear office computer
local_root_path = '/home/pierfier/Projects/';
server_root_path = '~/handata_server/';
% Local windows machine
%local_root_path = 'Z:\Local_Data\';
%server_root_path = 'X:\';
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
pv_data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Path to save the figures
%savefig_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f 'Exemplary' f];
% handata3 server path
%savefig_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'Plots' f];
savefig_path = Multi_func.save_plot();


% Check if the figure path exists
if ~exist(savefig_path)
    disp('Figure path not found');
    return;
end

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
%Load the data
load(save_all_data_file);
f_region = 'r_M1';
data_bystim = region_data.(f_region);

%% Info for a 140Hz neuron regular polarized
neuron = '617100_M1_rec20211110_FOV4_140_60';
trial_num = 3;
nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
neuron_data = data_bystim.f_140;

%% Info for a 40Hz example without delta
neuron = '50373_M1_rec20230801_FOV2_40_200';
trial_num = 8;
nr_idx = find(contains(data_bystim.f_40.neuron_name, neuron));
neuron_data = data_bystim.f_40;

%% Info for a 40Hz example with delta
neuron = '617100_M1_rec20211110_FOV3_40_60';
trial_num = 4;
nr_idx = find(contains(data_bystim.f_40.neuron_name, neuron));
neuron_data = data_bystim.f_40;

%% Infor for a 140Hz neuron hyperpolarized
neuron = '31556noeartag_M1_rec20221206_FOV1_140';
trial_num = 1;
nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
neuron_data = data_bystim.f_140;

%% Code to plot a heatmap with raster plot
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
tiledlayout(3, 1, 'TileSpacing', 'loose', 'Padding', 'loose', 'Units', 'centimeters', 'InnerPosition', [4 5 7 9]);
timeline = neuron_data.trace_timestamps(:, nr_idx);

%-- Plot example trial
nexttile;
raw_tr = neuron_data.all_trial_rawVm{nr_idx}(:, trial_num);
tr_noise = neuron_data.all_trial_trace_noise{nr_idx}(trial_num);
plot(timeline, raw_tr./tr_noise, 'k');
hold on;

spike_idx = neuron_data.all_trial_spikeidx{nr_idx}(:, trial_num);
%spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
%spike_idx = spike_idx - front_frame_drop + 1;
%spike_idx = spike_idx + 1;
spike_idx(isnan(spike_idx)) = [];
plot(timeline(spike_idx), ones(size(spike_idx)).*(max(raw_tr(spike_idx)./tr_noise) + 1), '.r');
hold on;

% Plot the stim pulses
stim_idx = neuron_data.stim_timestamps(:, nr_idx);
plot(stim_idx, ones(size(stim_idx))*(max(raw_tr(spike_idx)./tr_noise) + 3), '|k');
hold on;
plot(timeline, ones(size(timeline))*(max(raw_tr(spike_idx)./tr_noise) + 3), '-k');
hold on;
plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)

xlim([min(timeline) max(timeline)]);
ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);

%-- Plot the fluorescence heatmap
nexttile;
vm_map = [];
spidx_map = [];
% Loop through each trial of given neuron
for i = 1:size(neuron_data.all_trial_rawVm{nr_idx}, 2)
    vm_map = [vm_map; neuron_data.all_trial_rawVm{nr_idx}(:, i)'];
    cur_spikeidx = neuron_data.all_trial_spikeidx{nr_idx}(:, i);
    spidx_map = [spidx_map; cur_spikeidx'];
end
imagesc('XData', timeline, 'YData', 1:size(vm_map, 1), 'CData', vm_map);
colorbar;
hold on;

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(vm_map, 1)]);

ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

%-- Plot the population average trace
nexttile;

all_traces = neuron_data.all_trial_rawVm{nr_idx};
snr_traces = all_traces./neuron_data.all_trial_trace_noise{nr_idx};

% Calculate the average and SEM
avg_trace = nanmean(snr_traces, 2);
std_trace = nanstd(snr_traces, 0, 2);
num_trials = size(snr_traces, 2);
sem_trace = std_trace./sqrt(num_trials);

% Plot neuron's trial SEM
fill_h = fill([timeline; flip(timeline)], [avg_trace + sem_trace; flipud(avg_trace - sem_trace)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;

% Plot trial-averaged trace
plot(timeline, avg_trace, 'k');

hold on;

% Plotting the SNR scale
plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)

xlim([min(timeline) max(timeline)]);
ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);


% Save figure
saveas(gcf, [savefig_path 'Exemplary' f neuron '_heatmap.png']);
saveas(gcf, [savefig_path 'Exemplary' f neuron '_heatmap.pdf']);


%% Code to plot individual trials with raster plot and average
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
tiledlayout(3, 1, 'TileSpacing', 'loose', 'Padding', 'loose', 'Units', 'centimeters', 'InnerPosition', [4 5 8 15]);
timeline = neuron_data.trace_timestamps(:, nr_idx);


% Plot trial-averaged fluorescence
nexttile;
all_traces = neuron_data.all_trial_rawVm{nr_idx};
snr_traces = all_traces./neuron_data.all_trial_trace_noise{nr_idx};

% Calculate the average and SEM
avg_trace = nanmean(snr_traces, 2);
std_trace = nanstd(snr_traces, 0, 2);
num_trials = size(snr_traces, 2);
sem_trace = std_trace./sqrt(num_trials);

% Plot neuron's trial SEM
fill_h = fill([timeline; flip(timeline)], [avg_trace + sem_trace; flipud(avg_trace - sem_trace)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;

% Plot trial-averaged trace
plot(timeline, avg_trace, 'k');

hold on;

% Plot the stim pulses
stim_idx = neuron_data.stim_timestamps(:, nr_idx);
plot(stim_idx, ones(size(stim_idx))*(max(avg_trace + sem_trace) + 1), '|k');
hold on;
plot(timeline, ones(size(timeline))*(max(avg_trace + sem_trace) + 1), '-k');
hold on;

% Plotting the SNR scale
plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)

xlim([min(timeline) max(timeline)]);
ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);
% Plot each individual trial trace
nexttile;
trace_map = [];
spidx_map = [];
% Loop through each trial of given neuron
for i = 1:size(neuron_data.all_trial_rawVm{nr_idx}, 2)
    tr_trace = neuron_data.all_trial_rawVm{nr_idx}(:, i)';
    norm_trace = (tr_trace - min(tr_trace)) ./(max(tr_trace) - min(tr_trace));
    trace_map = [trace_map; norm_trace + (i - 0.5)];
    cur_spikeidx = neuron_data.all_trial_spikeidx{nr_idx}(:, i);
    spidx_map = [spidx_map; cur_spikeidx'];
end
plot(timeline, trace_map, 'k');
hold on;
% loop through and plot the spikes detected
for i =1:size(spidx_map, 1)
    temp_spidx = spidx_map(i, :);
    temp_spidx(isnan(temp_spidx)) = [];
    plot(timeline(temp_spidx), trace_map(i, temp_spidx), '.r');
    hold on;
end

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(trace_map, 1) + 0.5]);

ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

% Save figure
saveas(gcf, [savefig_path 'Exemplary' f neuron '_tr_traces_average.png']);
saveas(gcf, [savefig_path 'Exemplary' f neuron '_tr_traces_average.pdf']);



%% Code to plot a raster plot
% Plot the raster of the current neuron

raster_map = NaN(size(spidx_map));
raster_map(find(~isnan(spidx_map))) = 1;

raster_map = raster_map.*(1:size(spidx_map, 1))';

% Clean up the nans
spidx_map(isnan(spidx_map)) = [];
raster_map(isnan(raster_map)) = [];

% Plot all of the raster points
plot(timeline(spidx_map), raster_map, '.k');

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(trace_map, 1)]);
ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

%% Get exemplary trace at 140 for V1
example_matfile = [pv_data_path '611284_V1_rec20210827_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 5;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['611284_V1_rec20210827_FOV1_140_60_.mat tr 5' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 15;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);



title('Exemplary V1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 140Hz Power Spectra');


%% Get 2nd exemplary trace at 140 for V1
example_matfile = [pv_data_path '611284_V1_rec20221108_FOV2_140_50_.mat'];
data = load(example_matfile);
trial_idx = 2;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['611284_V1_rec20221108_FOV2_140_50_.mat tr 2' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 24;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary V1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 140Hz Power Spectra');

%% Get exemplary M1 trace at 140

% Parameters for plot zoom in
base_zoom = [-600 -400];
stim_zoom = [10 210];
offset_zoom = [1400 1600];

example_matfile = [pv_data_path '617100_M1_rec20211111_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 8;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure

figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
% Show file name
sgtitle(['617100_M1_rec20211111_FOV1_140_60_.mat tr 8'], 'Interpreter', 'none');

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 20;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 30, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -9;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-2, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);

% Plot boxes for each zoom in part
%baseline
dim = [];
dim(1) = base_idx(1);
dim(2) = -7;
dim(3) = range(base_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.base_color, 'LineStyle', '--');
hold on;

%stimulation
dim = [];
dim(1) = stim_ped_idx(1);
dim(2) = -7;
dim(3) = range(stim_ped_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.stim_color, 'LineStyle', '--');
hold on;

%offset
dim = [];
dim(1) = offset_idx(1);
dim(2) = -7;
dim(3) = range(offset_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.post_color, 'LineStyle', '--');

set(gca, 'Units', 'centimeters', 'Position', [5 20 10 3.00], 'PositionConstraint', 'innerposition');

title('Exemplary M1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.pdf']);


% Baseline zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
base_spike_idx = intersect(base_idx, spike_idx);

plot(detrend_traces(base_idx)./trace_noise, 'Color', Multi_func.base_color);
hold on;
plot(base_spike_idx - base_idx(1) + 1, detrend_traces(base_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy, ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 18]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Baseline Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_BaseZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_BaseZoomIn.pdf']);

% Stimulation zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
stim_spike_idx = intersect(stim_ped_idx, spike_idx);

plot(detrend_traces(stim_ped_idx)./trace_noise, 'Color', Multi_func.stim_color);
hold on;
plot(stim_spike_idx - stim_ped_idx(1) + 1, detrend_traces(stim_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Stimulation Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_StimZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_StimZoomIn.pdf']);


% Ofset zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
offset_spike_idx = intersect(offset_idx, spike_idx);

plot(detrend_traces(offset_idx)./trace_noise, 'Color', Multi_func.post_color);
hold on;
plot(offset_spike_idx - offset_idx(1) + 1, detrend_traces(offset_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Offset Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_OffsetZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_OffsetZoomIn.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 140Hz Power Spectra');

%------

%% Get 2nd exemplary M1 trace at 140
%example_matfile = [pv_data_path '617100_M1_rec20211110_FOV5_140_60_.mat'];
%trial_idx = 6;
example_matfile = [pv_data_path '617100_M1_rec20211110_FOV4_140_60_.mat'];
trial_idx = 3;
data = load(example_matfile);

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;


% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['617100_M1_rec20211110_FOV5_140_60_.mat tr 6'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 20;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -9;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary M1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.pdf']);

%% Get exemplary trace at 40 for V1
example_matfile = [pv_data_path '23072_V1_rec20220217_FOV3_40_220_.mat'];
data = load(example_matfile);
trial_idx = 3;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['23072_V1_rec20220217_FOV3_40_220_.mat tr 3'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 15;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -10;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('Exemplary V1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.pdf']);


%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 40Hz Power Spectra');


%% Get 2nd exemplary trace at 40 for V1
example_matfile = [pv_data_path '23072_V1_rec20220223_FOV2_40_250_.mat'];
data = load(example_matfile);
trial_idx = 8;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1

end 

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['23072_V1_rec20220223_FOV2_40_250_.mat tr 8'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 18;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -10;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('Exemplary V1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.pdf']);


%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 40Hz Power Spectra');

%% Get exemplary M1 trace at 40

% Parameters for plot zoom in
base_zoom = [-400 -200];
stim_zoom = [10 210];
offset_zoom = [1100 1300];

example_matfile = [pv_data_path '617100_M1_rec20211110_FOV3_40_60_.mat'];
data = load(example_matfile);
trial_idx = 9;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end
detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_idx = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_idx*sam_freq);
stim_start = data.raw.trial{trial_idx}.raw_stimulation_time(1);

% Calculate the idx time for base, stim, and offset
base_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > base_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < base_zoom(2)./1000);
stim_ped_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > stim_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < stim_zoom(2)./1000);
offset_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > offset_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < offset_zoom(2)./1000);

% Generate figure
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1000]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['617100_M1_rec20211110_FOV3_40_60_.mat tr 9' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 22;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);

% Plot boxes for each zoom in part
%baseline
dim = [];
dim(1) = base_idx(1);
dim(2) = -7;
dim(3) = range(base_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.base_color, 'LineStyle', '--');
hold on;

%stimulation
dim = [];
dim(1) = stim_ped_idx(1);
dim(2) = -7;
dim(3) = range(stim_ped_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.stim_color, 'LineStyle', '--');
hold on;

%offset
dim = [];
dim(1) = offset_idx(1);
dim(2) = -7;
dim(3) = range(offset_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.post_color, 'LineStyle', '--');

set(gca, 'Units', 'centimeters', 'Position', [5 20 9 3.00], 'PositionConstraint', 'innerposition');

title('Exemplary M1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.pdf']);


% Baseline zoom in of M1 40Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
base_spike_idx = intersect(base_idx, spike_idx);

plot(detrend_traces(base_idx)./trace_noise, 'Color', Multi_func.base_color);
hold on;
plot(base_spike_idx - base_idx(1) + 1, detrend_traces(base_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy, ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 18]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Baseline Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_BaseZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_BaseZoomIn.pdf']);

% Stimulation zoom in of M1 40Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
stim_spike_idx = intersect(stim_ped_idx, spike_idx);

plot(detrend_traces(stim_ped_idx)./trace_noise, 'Color', Multi_func.stim_color);
hold on;
plot(stim_spike_idx - stim_ped_idx(1) + 1, detrend_traces(stim_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Stimulation Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_StimZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_StimZoomIn.pdf']);

% Ofset zoom in of M1 140Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
offset_spike_idx = intersect(offset_idx, spike_idx);

plot(detrend_traces(offset_idx)./trace_noise, 'Color', Multi_func.post_color);
hold on;
plot(offset_spike_idx - offset_idx(1) + 1, detrend_traces(offset_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Offset Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_OffsetZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_OffsetZoomIn.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 40Hz Power Spectra');


%% 2nd example M1 trace at 40
example_matfile = [pv_data_path '31556eartag_M1_rec20221214_FOV1_40_80_.mat'];
data = load(example_matfile);
trial_idx = 4;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces,...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces,...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces);
end
detrend_traces = data.raw.trial{trial_idx}.raw_traces - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_idx = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_idx*sam_freq);

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['31556eartag_M1_rec20221214_FOV1_40_80_.mat tr 4'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 18;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary M1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 40Hz Power Spectra');
