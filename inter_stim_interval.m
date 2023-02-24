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

% Smoothing parameter for spike rate
srate_win = 50;

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
        data_bystim.(f_stim).neuron_base_Vm = [];
        data_bystim.(f_stim).neuron_stim_Vm = [];
        data_bystim.(f_stim).neuron_offset_Vm = [];
        
        data_bystim.(f_stim).neuron_base_srate = [];
        data_bystim.(f_stim).neuron_stim_srate = [];
        data_bystim.(f_stim).neuron_offset_srate = [];
               
        %data_bystim.(f_stim).stim_timestamps = [];
        %data_bystim.(f_stim).trace_timestamps = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_base_srate = [];
            cur_fov_base_Vm = [];
            cur_fov_stim_srate = [];
            cur_fov_stim_Vm = [];
            cur_fov_offset_srate = [];
            cur_fov_offset_Vm = [];

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
                    cur_trace_ws = trial_data.spike_info.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    cur_stim_time = raw_trial_data.raw_stimulation_time;
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop);
                    
                    % Find the idxs of the baseline, stimulation, and offset periods
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    stim_idx = find(cur_trace_time >= cur_stim_time(1) & cur_trace_time <= cur_stim_time(end));
                    offset_idx = find(cur_trace_time > cur_stim_time(end));

                    cur_fov_base_Vm(end + 1) = nanmean(cur_trace_ws(baseline_idx));   
                    cur_fov_stim_Vm(end + 1) = nanmean(cur_trace_ws(stim_idx));   
                    cur_fov_offset_Vm(end + 1) = nanmean(cur_trace_ws(offset_idx));   
                    

                    % Calculate the spike rate during baseline and stimulation period
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                                       
                    cur_fov_base_srate(end + 1) = sum(cur_raster(baseline_idx))./range(cur_trace_time(baseline_idx));
                    cur_fov_stim_srate(end + 1) = sum(cur_raster(stim_idx))./range(cur_trace_time(stim_idx));
                    cur_fov_offset_srate(end + 1) = sum(cur_raster(offset_idx))./range(cur_trace_time(offset_idx));

                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_base_Vm(:))) || isempty(cur_fov_base_Vm)
                continue;
            end

            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_base_Vm;
            data_bystim.(f_stim).neuron_base_Vm = horzcat_pad(temp, nanmean(cur_fov_base_Vm));

        end % End looping through FOVs of a condition
    end

%    % Save the VM to the specific region
%    region_data.(f_region) = data_bystim;

%end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%% Summary baseline vs stim vs offset properties
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 500 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    violin([data_bystim.(f_stim{1}).neuron_base_Vm', data_bystim.(f_stim{1}).neuron_stim_Vm', data_bystim.(f_stim{1}).neuron_offset_Vm'], 'xlabel', {'Base', 'Stim','Offset'});
    hold on;
    plot(repmat(1, length(data_bystim.(f_stim{1}).neuron_base_Vm), 1), data_bystim.(f_stim{1}).neuron_base_Vm, 'ko');
    hold on;
    plot(repmat(2, length(data_bystim.(f_stim{1}).neuron_stim_Vm), 1), data_bystim.(f_stim{1}).neuron_stim_Vm, 'ko');
    hold on;
    plot(repmat(3, length(data_bystim.(f_stim{1}).neuron_offset_Vm), 1), data_bystim.(f_stim{1}).neuron_offset_Vm, 'ko');
  
    ylabel('Vm (A.U.)');
    legend('off');
    
    % Statisitics for subthreshold
    disp('Subthreshold Vm statistics');
    f_stim{1}
    [h,p,ci,stats] = ttest(data_bystim.(f_stim{1}).neuron_base_Vm', data_bystim.(f_stim{1}).neuron_stim_Vm')
    title([f_stim{1}(3:end) ' p-val: ' num2str(p)], 'Interpreter', 'none');
end
sgtitle('Average baseline, stim, and offset Vm');

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
