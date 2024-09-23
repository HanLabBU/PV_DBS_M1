% Using the local aligned data, perform an interpolation of the fluorescence traces at the theoretical sampling frequency
% And see if the population Vm average is different
clc;
clear all;
close all;

f = filesep;

% read in all timestamp matfiles

set(0,'DefaultFigureVisible','on');
addpath('..');

local_root_path = '~/Projects/';

ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

%ses = dir(['*617100*_140*']);
%ses = dir(['*617100*_40*']);
ses = dir(['*_V1_*_140*']);
%ses = dir(['*_V1_*_40*']);
all_matfiles = {ses.name};

stim_start_diff = [];
frame_start_diff = [];

num_frames_prestim = [];

% Theoretical sampling rate
Fs = 1/(1.2/1000);
period = 1.2;

front_frame_drop = 6;

probe_time = [[-1000:period:0], [0:period:2000] ];
pop_avg = [];

for i = 1:length(all_matfiles)
    matfile = all_matfiles{i};
    data = load(matfile);
    
    fname = erase(matfile, 'align_wspce_');
    ri = strsplit(fname, '_');
    
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

    % Note everything will be referenced by stimulation start,
    % This means that I will have 1 sec of trace points before and 2 sec of
    % points after stimulation onset for each trial (that way it will be consistented throughout)
    % ALL UNITS ARE IN MS

    cur_fov_trace = [];

    % Interpolate Raw Fluorescence traces
    figure;
    tiledlayout(length(data.align.trial), 1);
    for tr = trial_idxs
        align_trial = data.align.trial{tr}.camera_framerate;
        % Raw traces
        tr_trace = data.raw.trial{tr}.raw_traces';

        stim_start = data.raw.trial{tr}.raw_stimulation_time(1)*1000;
        frame_time = data.align.trial{tr}.camera_frame_time*1000 - stim_start;
        
        interp_trace = interp1(frame_time, tr_trace, probe_time);
        interp_trace(1:front_frame_drop) = [];

        [baseline, coeff] = Multi_func.exp_fit_Fx(interp_trace', round(Fs));
        detrend_subVm = interp_trace - baseline;
        cur_fov_trace(:, end + 1) = detrend_subVm';

        %nexttile;
        %% Original trace plot with original frame time
        %plot(frame_time, tr_trace);
        %hold on;
        %% Interpolated trace
        %plot(probe_time, interp_trace + 50);

    end
    pop_avg(:, end + 1) = mean(cur_fov_trace, 2);

end

close all;
figure;
probe_time(1:front_frame_drop) = [];
plot(probe_time, mean(pop_avg, 2));
hold on;
xline(0:7.14:1000, 'b');
%hold on;
%xline(0:25:1000, 'r');
