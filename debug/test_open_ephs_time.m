clc;
clear all;
close all;

f = filesep;

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
addpath('..');

pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

front_frame_drop = 15 + round((828*.200));
back_frame_drop = 2496;

% Calculate population average of all M1
%ses = dir([pv_data_path '*617100*_140*']);
%ses = dir([pv_data_path '*617100*_40*']);
%ses = dir([pv_data_path '*_V1_*_140*']);
ses = dir([pv_data_path '*_V1_*_40*']);
matfiles = {ses.name};
pop_avg = [];
pre_pop_avg = [];
pre_frame_time = [];
time = [5:2500]*1.2;
all_stim_time = [];
all_frame_time = [];
%figure;
n = 1;
for i=1:length(matfiles)
    matfile = matfiles{i};
    data = load([pv_data_path matfile]); 

    ri = strsplit(matfile, '_');
    trial_idxs = find(~cellfun(@isempty, data.align.trial));
    try
        trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI1']);
    catch
        trial_ignr_list = [];
    end

    trial_idxs = setdiff(trial_idxs, trial_ignr_list);

    if length(trial_idxs) <= 2
        continue;
    end

    cur_fov_trace = [];
    pre_cur_fov_trace = [];
    cur_fov_stim_time = [];
    cur_fov_frame_time = [];
    pre_cur_fov_frame_time = [];
    matfile
    for j = trial_idxs

        align_trial = data.align.trial{j};
        raw_trial = data.raw.trial{j};
        pre_trace = align_trial.spike_info375.trace_ws(1, :);
        [baseline, coeff] = Multi_func.exp_fit_Fx(pre_trace', round(align_trial.camera_framerate));
        detrend_subVm = pre_trace - baseline;
        pre_cur_fov_trace(:, end + 1) = detrend_subVm';

        cur_trace = align_trial.spike_info375.trace_ws(1, front_frame_drop:back_frame_drop);
        
        [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace', round(align_trial.camera_framerate));
        detrend_subVm = cur_trace - baseline;
        cur_fov_trace(:, end + 1) = detrend_subVm';

        stim_start = raw_trial.raw_stimulation_time(1);
        cur_fov_stim_time(:, end + 1) = raw_trial.raw_stimulation_time(1:str2num(ri{5}))' - stim_start;
        cur_fov_frame_time(:, end + 1) = align_trial.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
        pre_cur_fov_frame_time(:, end + 1) = align_trial.camera_frame_time - stim_start;

    end
    cur_trace = mean(cur_fov_trace, 2);
    pre_trace = mean(pre_cur_fov_trace, 2);
    %min_trace = min(cur_trace);
    %max_trace = max(cur_trace);
    %norm_trace = (cur_trace - min_trace)./(max_trace - min_trace);

    pop_avg(:, end + 1) = cur_trace;
    pre_pop_avg(:, end + 1) = pre_trace;
    all_stim_time(:, end + 1) = mean(cur_fov_stim_time, 2, 'omitnan');
    all_frame_time(:, end + 1) = mean(cur_fov_frame_time, 2, 'omitnan');
    pre_frame_time(:, end + 1) = mean(pre_cur_fov_frame_time, 2);
    figure;
    % Plot all of the camera frames
    for j = 1:length(data.align.trial)
        n = j;
        %n = n + 1;
        align_trial = data.align.trial{j};
        
        raw_trial = data.raw.trial{j};
        %stim_start = raw_trial.raw_stimulation_time(1)*1000;
        camera_start = raw_trial.raw_camera_start_time*1000;
        plot(align_trial.camera_frame_time*1000 - camera_start, repmat(n, 1, length(align_trial.camera_frame_time)), '-k');
        hold on;
        
        norm_trace = norm(align_trial.detrend_traces);
        plot(align_trial.camera_frame_time*1000 - camera_start, norm_trace + n, '-g');  
        hold on;
        plot(raw_trial.raw_stimulation_time*1000 - camera_start, repmat(n+ 0.5, 1, length(raw_trial.raw_stimulation_time)), '|m');
        hold on;
    end
    
    xline(1000);
end
plot(time, repmat(0, 1, length(time)), '-r');
xlim([800 1050]);

close all;

figure;
plot(mean(all_frame_time, 2, 'omitnan'), mean(pop_avg, 2))
hold on;
xline(mean(all_stim_time, 2, 'omitnan') ); %+ 10*10^(-3)
hold on;
xline(0, 'r')
%xlim([800 1050]);
title('Population Average');

% Plotting average before chopping
figure;
plot(mean(pre_frame_time, 2, 'omitnan'), mean(pre_pop_avg, 2))
hold on;
xline(mean(all_stim_time, 2, 'omitnan') ); %+ 10*10^(-3)
hold on;
xline(0, 'r')
%xlim([800 1050]);
title('Pre chopping Population Average');

% SUGGESTION: Use cameraframes to plot the traces
return;

figure;
for i =1:size(pop_avg, 2)
    plot(time, pop_avg(:, i) + i - 1, 'k');
    hold on;
end
xline(1000);




return;
% Read in all V1 neurons


ses = dir([pv_data_path '*_V1_*140*.mat']);
V1_matfiles = {ses.name};

all_v1_traces = [];

% Loop through each matfile
for i = 1:length(V1_matfiles)
    matfile = V1_matfiles{i};
    
    data = load([pv_data_path matfile]);

    % Raw timestamps
    %figure;
    %% Loop through each trial
    %for j = 1:length(data.raw.trial)
    %    raw_trial = data.raw.trial{j};
    %    camera_start = raw_trial.raw_camera_start_time;

    %    % Plot Camera timestamps
    %    plot(raw_trial.raw_camera_frame_time - camera_start, repmat(j, 1, length(raw_trial.raw_camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    %% Aligned timestamps
    %figure;
    %for j = 1:length(data.align.trial)
    %    raw_trial = data.raw.trial{j};
    %    align_trial = data.align.trial{j};
    %    camera_start = align_trial.camera_start_time;
    %    
    %    plot(align_trial.camera_frame_time - camera_start, repmat(j, 1, length(align_trial.camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    % Plot traces with stimulation 
    %figure;
    cur_fov_trace = [];
    for j = 1:length(data.align.trial)
        align_trial = data.align.trial{j};
        time = [5:2500]*1.2;
        min_trace = min(align_trial.detrend_traces);
        max_trace = max(align_trial.detrend_traces);
        norm_trace = (align_trial.detrend_traces - min_trace)./(max_trace - min_trace);
        %plot(time, norm_trace + j);
        %hold on;

        % collect for this neuron
        cur_fov_trace(:, end + 1) = align_trial.detrend_traces;
    end
    all_v1_traces(:, end + 1) = mean(cur_fov_trace, 2, 'omitnan');
    %xline(1000);
end

% Plot the population V1 Vm
figure;
time = [5:2500]*1.2;
plot(time, mean(all_v1_traces, 2, 'omitnan'))
hold on;
xline(1000);
title('Pop average 140');

% Test the traces just from calculating the average exposure rate for each frame

ses = dir([pv_data_path '*_V1_*40*.mat']);
V1_matfiles = {ses.name};

all_v1_traces = [];

% Loop through each matfile
for i = 1:length(V1_matfiles)
    matfile = V1_matfiles{i};
    
    data = load([pv_data_path matfile]);


    % Raw timestamps
    %figure;
    %% Loop through each trial
    %for j = 1:length(data.raw.trial)
    %    raw_trial = data.raw.trial{j};
    %    camera_start = raw_trial.raw_camera_start_time;

    %    % Plot Camera timestamps
    %    plot(raw_trial.raw_camera_frame_time - camera_start, repmat(j, 1, length(raw_trial.raw_camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    %% Aligned timestamps
    %figure;
    %for j = 1:length(data.align.trial)
    %    raw_trial = data.raw.trial{j};
    %    align_trial = data.align.trial{j};
    %    camera_start = align_trial.camera_start_time;
    %    
    %    plot(align_trial.camera_frame_time - camera_start, repmat(j, 1, length(align_trial.camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    % Plot traces with stimulation 
    %figure;
    cur_fov_trace = [];
    for j = 1:length(data.align.trial)
        align_trial = data.align.trial{j};
        time = [5:2500]*1.2;
        min_trace = min(align_trial.detrend_traces);
        max_trace = max(align_trial.detrend_traces);
        norm_trace = (align_trial.detrend_traces - min_trace)./(max_trace - min_trace);
        %plot(time, norm_trace + j);
        %hold on;

        % collect for this neuron
        cur_fov_trace(:, end + 1) = align_trial.detrend_traces;
    end
    all_v1_traces(:, end + 1) = mean(cur_fov_trace, 2, 'omitnan');
    %xline(1000);
end

% Plot the population V1 Vm
figure;
time = [5:2500]*1.2;
plot(time, mean(all_v1_traces, 2, 'omitnan'))
hold on;
xline(1000);
title('Pop average 40');

% Read in all 617100 neurons

pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
ses = dir([pv_data_path '*617100*.mat']);
V1_matfiles = {ses.name};

all_v1_traces = [];

% Loop through each matfile
for i = 1:length(V1_matfiles)
    matfile = V1_matfiles{i};
    
    data = load([pv_data_path matfile]);

    % Raw timestamps
    %figure;
    %% Loop through each trial
    %for j = 1:length(data.raw.trial)
    %    raw_trial = data.raw.trial{j};
    %    camera_start = raw_trial.raw_camera_start_time;

    %    % Plot Camera timestamps
    %    plot(raw_trial.raw_camera_frame_time - camera_start, repmat(j, 1, length(raw_trial.raw_camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    %% Aligned timestamps
    %figure;
    %for j = 1:length(data.align.trial)
    %    raw_trial = data.raw.trial{j};
    %    align_trial = data.align.trial{j};
    %    camera_start = align_trial.camera_start_time;
    %    
    %    plot(align_trial.camera_frame_time - camera_start, repmat(j, 1, length(align_trial.camera_frame_time)), 'k');
    %    hold on;
    %    plot(raw_trial.raw_stimulation_time - camera_start, repmat(j, 1, length(raw_trial.raw_stimulation_time)), 'g');
    %    hold on;
    %end
    %xline(1);
    %xlim([-1 4]);

    % Plot traces with stimulation 
    %figure;
    cur_fov_trace = [];
    for j = 1:length(data.align.trial)
        align_trial = data.align.trial{j};
        time = [5:2500]*1.2;
        min_trace = min(align_trial.detrend_traces);
        max_trace = max(align_trial.detrend_traces);
        norm_trace = (align_trial.detrend_traces - min_trace)./(max_trace - min_trace);
     %   plot(time, norm_trace + j);
     %   hold on;

        % collect for this neuron
        cur_fov_trace(:, end + 1) = align_trial.detrend_traces;
    end
    all_v1_traces(:, end + 1) = mean(cur_fov_trace, 2, 'omitnan');
    %xline(1000);

    figure;
    plot(time, mean(cur_fov_trace, 2, 'omitnan'));
    hold on;
    xline(1000);
    title(['M1 ' matfile]);
    % Plot the neurons average trace
end

% Plot the population V1 Vm
figure;
time = [5:2500]*1.2;
plot(time, mean(all_v1_traces, 2, 'omitnan'))
hold on;
xline(1000);
title('Pop average 140');


function [result] = norm(trace)
    min_trace = min(trace);
    max_trace = max(trace);
    result = (trace - min_trace)./(max_trace - min_trace);
end
