clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
server_root_path = '~/Projects/';
% Windows server
%server_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine how much before and after a stimulation pulse to take the average
trace_sur = 10; % This is ~6ms before and after stim pulse

%%% END Modification

% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);

all_matfiles = {ses.name};

% Select matfiles by brain region
[region_matfiles] = Multi_func.find_region(all_matfiles);
region_data = struct();

%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%    %% Select matfiles by stim specific conditions for all regions
    [matfile_stim] = stim_cond(all_matfiles); %stim_cond(region_matfiles.(f_region).names);
    %% Loop through each field of the struct and concatenate everything together
    
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies
    all_Fs = [];

    % Loop through each stimulation condition
    for f_stim = fieldnames(matfile_stim)'
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    
    
        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();
        data_bystim.(f_stim).neuron_stim_avg = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_base_inter = [];

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    %Determine whether this roi is to be ignored for this particular trial
                    ri = strsplit(matfile{1}, '_');
                    try
                        trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                    catch
                        trial_ignr_list = [];
                    end
    
                    % Check if current trial is in the ignore list
                    if ismember(tr_idx, trial_ignr_list)
                        continue;
                    end
                    
                    % If the trial data is empty, that means it was skipped
                    if isempty(trial_data)
                        continue;
                    end

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;

                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace = trial_data.detrend_traces(front_frame_drop:back_frame_drop, roi_idx);
                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{5}));
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop);
                    
                    % Get the trace idx when there are stimulation pulses
                    % Cannot think of a clever way to vectorize this unfortunately
                    for pulse=cur_stim_time'
                        stim_center = find(min(abs(pulse - cur_trace_time)) == abs(pulse - cur_trace_time));
                        cur_fov_stim_avg = horzcat_pad(cur_fov_stim_avg, cur_trace(stim_center - trace_sur: stim_center + trace_sur));
                    end
                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_stim_avg(:))) || isempty(cur_fov_stim_avg)
                continue;
            end

            % Plot the average and all of the trace from the stim centers
            figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 500 1000]);
            plot(cur_fov_stim_avg, 'color', [95, 78, 78]./255);
            hold on;
            plot(mean(cur_fov_stim_avg, 2), 'r');
            hold on;
            xline(size(cur_fov_stim_avg, 1)/2, '--b');
            title([matfile{1}], 'Interpreter', 'none');
            saveas(gcf, [figure_path 'Stim_trig_avg' f matfile{1}(1:end-4) '_stim_avg.png']);
        end % End looping through FOVs of a condition
    end

%    % Save the VM to the specific region
%    region_data.(f_region) = data_bystim;

%end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%% Summary baseline vs stim vs offset properties
%stims = fieldnames(data_bystim);
%figure('Renderer', 'Painters', 'Position', [200 200 500 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%    title([f_stim{1}(3:end)], 'Interpreter', 'none');
%end
%sgtitle('Average baseline, stim, and offset Vm');

%% Specific functions for determining which FOVs to look at
% Return matfiles by stimulation condition
function [cond_struct] = stim_cond(matfile_names)
    cond_struct = struct();
    
    % Loop through each matfilename and group by stimulation conditions
    for i=1:length(matfile_names)
            file_parts = split(matfile_names{i}, '_');
            stim = file_parts{5};
            
            % Create stimulation field if it does not exist
            if ~isfield(cond_struct, ['f_' stim])
                cond_struct.(['f_' stim]).names = {};
            end

            cond_struct.(['f_' stim]).names{end+1} = matfile_names{i};
    end
end
