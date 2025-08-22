clear all;
close all;
clc;

% data path
%data_path = [server_root_path ];
data_path = '/home/pierfier/Projects/Pierre Fabris/PV DBS neocortex/CA1_Data/';

ses = dir([data_path '*mat']);
matfiles = {ses.name};

for cur_mat = matfiles %TODO just for testing here
    cur_mat = cur_mat{1}
    data = load([data_path cur_mat]);
    
    % If lfp does not exist, skip
    if ~isfield(data.result, 'lfp')
        continue;
    end
    
    % Replace the index for trials vectors
    trials = unique(data.result.trial_vec);
    if trials(1) ~= 1
        data.result.trial_vec = data.result.trial_vec - trials(1) + 1;
    end
    trials = unique(data.result.trial_vec);

    num_trace_pts = sum(data.result.trial_vec == data.result.trial_vec(1));
    
    ri = strsplit(cur_mat, '_');

    %TODO rsplit for the frequency points

    [frame_timestamps, stim_timestamps] = grab_camera_timestamps(data.result.lfp, num_trace_pts, str2num(ri{4}), 10000);

    if length(stim_timestamps) == 0
        continue;
    end
        
    trial_data = {};
    raw_trial_data = {};

    % Loop through each trials
    trials = unique(data.result.trial_vec);
    
    if length(frame_timestamps) < length(trials)
        continue
    end

    for tr = trials
        trial_data{tr}.camera_framerate =  data.result.lfp.camera_FS;   
        trial_data{tr}.camera_frame_time = frame_timestamps{tr};
        trial_data{tr}.spike_info = data.result.resultS{tr};
        
        raw_trial_data{tr}.raw_stimulation_time = stim_timestamps{tr};
        raw_trial_data{tr}.raw_traces = data.result.traces(data.result.trial_vec == tr);
        orient = size(raw_trial_data{tr}.raw_traces);
        if orient(2) ~= 1
            raw_trial_data{tr}.raw_traces = raw_trial_data{tr}.raw_traces';
        end

    end
    align.trial = trial_data;
    raw.trial = raw_trial_data;

    % Append value to matfile
    save([data_path cur_mat], 'align', 'raw' , '-append');
    disp(cur_mat);
    disp('Saved');
end

function [camera_frame_timestamps, stim_timestamps] = grab_camera_timestamps(ephs_data, num_trace_pts, stim_freq, ephys_rate)
    stim_timestamps = {};
    camera_frame_timestamps = {};
    full_ephys_time = [1:length(ephs_data.raw_camera_trig)]/ephys_rate;

    %figure;
    %plot(ephs_data.raw_camera_trig);
    %hold on;
    %plot(ephs_data.raw_stim_trig);

    % If there are no trigger data move on
    if length(ephs_data.stim_trigger_ind) < stim_freq || length(ephs_data.camera_trigger_ind) < num_trace_pts
        disp('Missing timestamps');
        return;
    end

    % Missing start time
    if ~isfield(ephs_data, 'start_only_cMOS')
        ephs_data.start_only_cMOS = find_trial_starts(full_ephys_time, ephs_data.camera_trigger_ind);

        %DEBUG
        figure;
        %TODO plot the indices instead!!!!
        plot(ephs_data.camera_trigger_ind, ones(size(ephs_data.camera_trigger_ind)), '.b');
        hold on;
        plot(ephs_data.stim_trigger_ind, ones(size(ephs_data.stim_trigger_ind)), '.g');
        hold on;
        plot(ephs_data.start_only_cMOS, ones(size(ephs_data.start_only_cMOS)),'|r');
    end

    % Calculate the approximate frame rate and see if it fits the past framerate
    first_trial_diff = diff(ephs_data.camera_trigger_ind(1:num_trace_pts- 50));
    emp_cal_rate = mean(first_trial_diff)/(1/ephs_data.camera_FS);

    % Check ephys rate difference for frame timestamps
    if abs(emp_cal_rate - ephys_rate) > 500
        disp('There is a difference in ephys sampling frequency!!!');
        disp(['Calculated ephys rate ' num2str(emp_cal_rate) ' passed rate ' num2str(ephys_rate)]);
        return;
    end

    % First stim difference of stimulation timestamps
    first_stim_diff = diff(full_ephys_time(ephs_data.stim_trigger_ind(1:stim_freq)));
    emp_stim_freq = 1/mean(first_stim_diff);
    
    if abs(emp_stim_freq - stim_freq) > 3
        disp(['Ephys stim frequency and specified frequency do not match']);
        disp(['Ephys stim ' num2str(emp_stim_freq) ' specified ' num2str(stim_freq)]);
        return;
    end

    for tr_start_idx = ephs_data.start_only_cMOS';
        tr_start_idx;
        trial_cam_idx = ephs_data.camera_trigger_ind(ephs_data.camera_trigger_ind > tr_start_idx);
        
        % Check if there are enough camera timestamps to look at trial
        if length(trial_cam_idx) < num_trace_pts
            continue;
        end

        trial_cam_idx = trial_cam_idx(1:num_trace_pts);
        trial_stim_idx = ephs_data.stim_trigger_ind(ephs_data.stim_trigger_ind > tr_start_idx);
        
        trial_stim_idx = trial_stim_idx(1:stim_freq);

        camera_frame_timestamps = [camera_frame_timestamps, {full_ephys_time(trial_cam_idx)'}];
        stim_timestamps = [stim_timestamps, {full_ephys_time(trial_stim_idx)'}];
    end
end

% Finds the trial start indices
function [trial_start_idx] = find_trial_starts(ephys_time, camera_frame_idx)
    % Account for frame rate error
    frame_diff = [0, diff(ephys_time(camera_frame_idx))];
    zscore_diff = zscore(frame_diff);
    trial_inflect_idx = find(zscore_diff > 2);

    trial_start_idx = [camera_frame_idx(1); camera_frame_idx(trial_inflect_idx) ] - 3;
end
