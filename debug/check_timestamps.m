clc;
clear all;
close all;

% read in all timestamp matfiles

set(0,'DefaultFigureVisible','off');

ses = dir(['*.mat']);
all_matfiles = {ses.name};

stim_start_diff = [];
frame_start_diff = [];

num_frames_prestim = [];

for i = 1:length(all_matfiles)
    matfile = all_matfiles{i};
    data = load(matfile);

    % Plot all of the start timestamps
    figure('Position', [0 0 1500 1200]);
    tiledlayout(3, 1);

    % Plotting start time
    ax1 = nexttile;
    plot(data.all_raw_cam_start_time, data.all_raw_cam_start_data);
    hold on;
    
    % Aligned points start
    for tr=1:length(data.align.trial)
        plot(data.align.trial{tr}.camera_start_time, 1, '|', 'LineWidth', 2);
        hold on;
    end
    title('Camera trial recording start');

    % Rise times
    %plot(data.all_cam_start_times, repmat(1, 1, length(data.all_cam_start_times)), '|', 'LineWidth', 2);
    
    % Plotting camera frame time
    ax2 = nexttile;
    
    plot(data.all_raw_cam_frame_time, data.all_raw_cam_frame_data);
    hold on;
    avg = mean([min(data.all_raw_cam_frame_data), max(data.all_raw_cam_frame_data)]);
    % Aligned points start
    for tr=1:length(data.align.trial)
        plot(data.align.trial{tr}.camera_frame_time, repmat(avg, 1, length(data.align.trial{tr}.camera_frame_time)), '|', 'LineWidth', 2);
        hold on;
    end
    
    title('Camera frame times');

    % Rise times
    %plot(data.all_cam_frame_times, repmat(1, 1, length(data.all_cam_frame_times)), '|', 'LineWidth', 2);
    
    % plotting stime time
    ax3 = nexttile;
    plot(data.all_raw_stim_time, data.all_raw_stim_data);
    hold on;

    % Aligned stim times
    for tr=1:length(data.raw.trial)
        plot(data.raw.trial{tr}.raw_stimulation_time, repmat(1, 1, length(data.raw.trial{tr}.raw_stimulation_time)), '|', 'LineWidth', 2);
        hold on;
    end

    title('Stimulation time');

    % Rise times
    %plot(data.all_stim_times, repmat(1, 1, length(data.all_stim_times)), '|', 'LineWidth', 2);
    linkaxes([ax1, ax2, ax3], 'x');
    sgtitle(matfile, 'Interpreter', 'none');


    % Checking initial start values

    % Aligned points start
    for tr=1:length(data.align.trial)
        
        % Skip trials taht do not have stimulation time
        if isempty(data.raw.trial{tr}.raw_stimulation_time)
            disp('No stim timestamps');
            disp(matfile);
            disp(['Trial #' num2str(tr)]);
            continue;
        end

        stim_start_diff(end + 1) = data.align.trial{tr}.camera_start_time - data.raw.trial{tr}.raw_stimulation_time(1);
        
        frame_start_diff(end + 1) = data.align.trial{tr}.camera_start_time - data.align.trial{tr}.camera_frame_time(1);

        unique(diff(data.align.trial{tr}.camera_frame_time)*1000);

        % Calculate number of frames before stim
        pre_idx = find(data.align.trial{tr}.camera_frame_time < data.raw.trial{tr}.raw_stimulation_time(1));
        num_frames_prestim(end + 1) = length(pre_idx);
    end
end

figure('visible', 'on');
plot(stim_start_diff*1000, '.', 'MarkerSize', 10);
title('Stimulation start from trial start');

figure('visible', 'on');
plot(frame_start_diff*1000, '.', 'MarkerSize', 10);
title('Frame start from trial start');

