% Housekeeping clear all stuff
clc;
clear all;
close all;
f = filesep;

% Specify which photobleaching detrending to use
% True (1) for the photobleaching taking into account the stimulation bump
% False (0) for the basic full trace expontential fit subtraction detrending
sophis_bleachdetrend = 1;

% Maingear office computer
local_root_path = '/home/pierfier/Projects/';
server_root_path = '~/handata_server/';
% Local linux machine
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
pv_data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Path to save the figures
%savefig_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f 'Exemplary' f];
savefig_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'Plots' f];

%% Get exemplary trace at 140 for V1
example_matfile = [pv_data_path '611284_V1_rec20210827_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 5;

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.eps'], 'epsc');
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

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.eps'], 'epsc');
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
example_matfile = [pv_data_path '617100_M1_rec20211111_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 8;

if sophis_bleachdetrend
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
title('Exemplary M1 140 Hz trace');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.pdf']);

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
example_matfile = [pv_data_path '617100_M1_rec20211110_FOV5_140_60_.mat'];
data = load(example_matfile);
trial_idx = 6;

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.pdf']);

%% Get exemplary trace at 40 for V1
example_matfile = [pv_data_path '23072_V1_rec20220217_FOV3_40_220_.mat'];
data = load(example_matfile);
trial_idx = 3;

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.eps'], 'epsc');
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

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.eps'], 'epsc');
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
example_matfile = [pv_data_path '617100_M1_rec20211110_FOV3_40_60_.mat'];
data = load(example_matfile);
trial_idx = 9;

if sophis_bleachdetrend
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
title('Exemplary M1 40 Hz trace');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.pdf']);

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

if sophis_bleachdetrend
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
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.eps'], 'epsc');
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
