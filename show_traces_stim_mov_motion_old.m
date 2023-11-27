% Housekeeping clear all stuff
clc;
clear all;
close all;

f = filesep;

%------------ USER modification
% Server root path
%server_rootpath = 'Z:\';

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';

% Server folder location of the saved and aligned data
%data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Local linux machine
%data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

% Data share on server
data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
% Data on local computer
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

% Determine whether to save the ignored traces (1) or not ignored traces (0)
show_ignored = 0;
% USER make sure this path changes based on the above line
save_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Kept_Traces' f];

% Filepath name for ignoring individual trial csv
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);


%------------- END modification


% Get all FOV matfiles
matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);

    %DEBUG
    disp(matfile_names{i});
        
    % Loop through each trial
    for tr=1:length(data.align.trial)
        if isempty(data.align.trial{tr})
            continue;
        end

        % Check if trial is in the ignore list
        ri = strsplit(matfile_names{i}, '_');
        try
            trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).('ROI1'); % Always assuming ROI 1
        catch
            trial_ignr_list = [];
        end

        % Check if current trial is in the ignore list
        if ismember(tr, trial_ignr_list) ~= show_ignored
            continue;
        end

        % Setup figure to show alignment data for all trials
        figure('visible', 'off', 'Position', [0, 0, 800, 1000], 'visible', 'off');
        tiledlayout((size(data.align.trial{tr}.detrend_traces, 2)*2) + 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Save current trial
        trial_data = data.align.trial{tr};
        
        %EBUG
        %disp((size(data.align.trial{tr}.detrend_traces, 2)*2) + 2);
        %disp(size(trial_data.detrend_traces, 2));

        % Loop through each ROI
        for roi=1:size(trial_data.detrend_traces, 2)
            
            % Check if trace has NaNs
            if sum(~isnan(data.raw.trial{tr}.raw_traces(:, roi))) == 0
                continue;
            end

            % Plotting the raw trace
            nexttile;
            raw_trace = data.raw.trial{tr}.raw_traces(:, roi);
            plot(trial_data.camera_frame_time, raw_trace);
            hold on;
            
            % Plot the stimulation pattern
            plot(data.raw.trial{tr}.raw_stimulation_time, repmat(max(raw_trace), ...
                length(data.raw.trial{tr}.raw_stimulation_time)), '|m');
            
            title('Raw Trace');

            % Plotting detrended trace by exponential fit with stimulation pattern and spike detected points
            nexttile;
            [y x] = Multi_func.exp_fit_Fx(raw_trace(:), round(trial_data.camera_framerate));
            detrend_trace = raw_trace - y';
            plot(trial_data.camera_frame_time, detrend_trace, '-b');
            hold on
            % Plot the stimulation pattern
            plot(data.raw.trial{tr}.raw_stimulation_time, repmat(max(detrend_trace), ...
                length(data.raw.trial{tr}.raw_stimulation_time)), '|m');
            hold on;
            % Plot spikes detected
            spikes_idx = data.align.trial{tr}.spike_info375.spike_idx{1};
            plot(trial_data.camera_frame_time(spikes_idx), detrend_trace(spikes_idx), 'or', 'MarkerSize', 4);
            
            % Grab neuron's SNRs
            snrs = trial_data.spike_info375.spike_snr{1};

            legend(['SNR ' num2str(nanmean(snrs)) ' Num spikes: ' num2str(length(snrs))]);    
            title('Detrended Trace');
        end
        
        % Plotting motion correction vectors
        nexttile;
        plot(trial_data.camera_frame_time, trial_data.img_correct_vec);
        title('Motion Correction Vector Error');

        % Plot the movement, if the variable exists
        if isfield(trial_data, 'speed')
            nexttile;
            plot(trial_data.speed_timestamp, trial_data.speed);
            title('Raw Movement');
        end
        
        % Set the title of the Figure
        sgtitle([matfile_names{i} ' trial: ' num2str(tr)], 'Interpreter', 'none');
        
        % Save figure as a jpeg
        saveas(gcf, [save_path matfile_names{i}(1:end-4) num2str(tr) '.png']);
    end
    

end

% Perform exponential fit
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
