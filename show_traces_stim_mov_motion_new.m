% Housekeeping clear all stuff
clc;
clear all;
close all;

f = filesep;

%USER modification
% Server root path
server_rootpath = 'Z:\';

% Maingear office computer
%server_rootpath = '/home/pierfier/handata_server/';

% Folder location of the saved and aligned data
data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

%END modification

% Get all FOV matfiles
matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);

    % Loop through each trial
    for tr=1:length(data.align.trial)
        
        % Setup figure to show alignment data for all trials
        figure('Position', [0, 0, 800, 1000]);
        tiledlayout(size(data.align.trial{tr}.detrend_traces, 2)*2 + 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Temporarily save this trial
        trial_data = data.align.trial{tr};

        % Loop through each ROI
        for roi=1:size(trial_data)
            
            % Plotting the raw trace
            raw_trace = data.raw.trial{tr}.raw_traces(:, roi);
            nexttile;
            plot(trial_data.camera_frame_time, raw_trace);
            title('Raw Trace');

            % Plotting detrended trace by exponential fit with stimulation pattern and spike detected points
            [x y] = exp_fit(raw_trace(:), trial_data.camera_framerate);
            detrend_trace = raw_trace - y';
            nexttile;
            plot(trial_data.camera_frame_time, detrend_trace, '-b');
            hold on
            % Plot the stimulation pattern
            plot(data.raw.trial{tr}.raw_stimulation_time, repmat(max(detrend_trace), ...
                length(data.raw.trial{tr}.raw_stimulation_time)), '|m');
            hold on;
            % Plot spikes detected
            spikes_idx = data.align.trial{tr}.spike_info.spike_idx{1};
            plot(trial_data.camera_frame_time(spikes_idx), detrend_trace(spikes_idx), '.r', 'MarkerSize', 4);
            title('Detrended Trace');
        end
        
        % Plotting motion correction vectors
        nexttile;
        plot(trial_data.camera_frame_time, trial_data.img_correct_vec);

        % Plot the movement
        nexttile;
        plot(trial_data.speed_timestamp, trial_data.speed);
        title('Raw Movement');

        % Set the title of the Figure
        sgtitle([matfile_names{i} ' trial: ' num2str(tr)], 'Interpreter', 'none');
    end
    

end

% Perform exponential fit
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
