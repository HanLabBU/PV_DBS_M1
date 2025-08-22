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
exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end

back_frame_drop = 2496;

% Get all FOV matfiles
matfile_names = dir([data_path '*V1*.mat']);
matfile_names = {matfile_names.name};

% Loop through each matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);

    %DEBUG
    disp(matfile_names{i});
        
    temp_idx = find(~cellfun(@isempty, data.align.trial));

    % Loop through each ROI
    for roi_idx =1:size(data.align.trial{temp_idx(1)}.detrend_traces, 2)

    % Make figure for each ROI
    figure('Position', [0, 0, 800, 1000]);
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    nexttile;

        % Loop through each trial
        for tr=1:length(data.align.trial)
            if isempty(data.align.trial{tr})
                continue;
            end

            trial_data = data.align.trial{tr};
            raw_tr_data = data.raw.trial{tr};

            % Check if trial is in the ignore list
            ri = strsplit(matfile_names{i}, '_');
            try
                trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
            catch
                trial_ignr_list = [];
            end

            % Check if current trial is in the ignore list
            if ismember(tr, trial_ignr_list) ~= show_ignored
                continue;
            end
            
            stim_start = raw_tr_data.raw_stimulation_time(1);

            % Check that there is interp_time
            if isfield(trial_data, 'interp_time')
                cur_trace_time = trial_data.interp_time(front_frame_drop:back_frame_drop)- stim_start;                
                cur_vm = trial_data.interp_detrend_traces(front_frame_drop:back_frame_drop, roi_idx);
            else
                cur_trace_time = trial_data.camera_frame_time - stim_start;
                cur_vm = data.align.trial{tr}.detrend_traces(:, roi_idx);
            end

            norm_vm = (cur_vm - min(cur_vm))./(max(cur_vm) - min(cur_vm));
            
            % Plot each roi 23072_V1_rec20220301_FOV1_140_280_.mattrace
            plot(cur_trace_time, norm_vm + (tr*1.2));
            hold on;

            % Plot the spikes
            spike_idx = data.align.trial{tr}.spike_info375.spike_idx{roi_idx};
            
            spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
            spike_idx = spike_idx - front_frame_drop + 1;

            plot(cur_trace_time(spike_idx), norm_vm(spike_idx) + (tr*1.2), 'or', 'MarkerSize', 4);

        end % End of looping through trial
        
        cur_stim_start = data.raw.trial{temp_idx(1)}.raw_stimulation_time(1);
        % Plot the stimulation pattern
        plot(data.raw.trial{temp_idx(1)}.raw_stimulation_time - cur_stim_start, repmat(1.2*tr + 1.5, ...
            length(data.raw.trial{temp_idx(1)}.raw_stimulation_time)), '|m');


        % Plot the motion correction error
        % It is not interpolated, so it wont be perfect
        nexttile;
        plot(trial_data.camera_frame_time, trial_data.img_correct_vec);
        sgtitle([matfile_names{i} ' roi: ' num2str(roi_idx)], 'Interpreter', 'none');

        % Save figure as a jpeg
        saveas(gcf, [save_path matfile_names{i}(1:end-4) '_roi_' num2str(roi_idx) '.png']);

    end % End of looping through ROI
end % End of looping through all matfiles
