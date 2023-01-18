clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
server_root_path = '/home/pierfier/handata_server/';
% Windows server
%server_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = ['~/Projects' f 'PV DBS Project' f 'PV_Data' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate
srate_window = 26;

%%% END Modification

% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);

all_matfiles = {ses.name};

%% Select matfiles by specific conditions

[matfile_stim] = stim_cond(all_matfiles);

%% Loop through each field of the struct and concatenate everything together

% Store trace aspect data by each individual stimulation condition
data_bystim = struct();
% Store all of the calculated sampling frequencies
all_Fs = [];

% Loop through each stimulation condition
for field = fieldnames(matfile_stim)'
    field = field{1};
    matfiles = matfile_stim.(field).names;    

    % Initialize field subthreshold array
    data_bystim.(field) = struct();
    data_bystim.(field).sub_spectra_wt = [];
    data_bystim.(field).sub_spectra_f = [];
    

    % Loop through each matfile of the current stimulation condition
    for matfile = matfiles
        % Read in the mat file of the current condition
        data = load([pv_data_path matfile{1}]);
        

        cur_fov_wt = [];
        cur_fov_f = [];
        % Loop through each trial
        for tr_idx=1:length(data.align.trial)
            trial_data = data.align.trial{tr_idx};
            
            % If the trial data is empty, that means it was skipped
            if isempty(trial_data)
                continue;
            end
            
            % Store the camera framerate
            all_Fs(end+1) = trial_data.camera_framerate;
            
            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(matfile{1}, '_');
                try
                    roi_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                catch
                    roi_ignr_list = [];
                end
                if ismember(roi_idx, roi_ignr_list)
                    continue;
                end

                % Grab the subthreshold Vm
                % Chop the respective frames
                cur_trace_ws = trial_data.spike_info.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                [wt, f] = cwt(cur_trace_ws, all_Fs(end));
                
                cur_fov_wt = cat(3, cur_fov_wt, wt);
                cur_fov_f = [cur_fov_f, f];
                % Hopefully I do not need to pad the concatenation
                %cur_fov_subVm = horzcat_pad(cur_fov_subVm, cur_trace_ws');
            end
        end
        
        %EBUG
        %figure;
        %plot(nanmean(cur_fov_subVm, 2));

        % Average for each neuron and save the subthreshold Vm
        temp = data_bystim.(field).sub_spectra_wt;
        data_bystim.(field).sub_spectra_wt = cat(3, temp, cur_fov_wt);
        data_bystim.(field).sub_spectra_f = [data_bystim.(field).sub_spectra_f, cur_fov_f]; 

    end % End looping through FOVs of a condition
end 

% Plot all of the subthreshold Vm for each stimulation condition
%TODO find the stimulation onset with the timestamps from the stimulation and camera
stims = fieldnames(data_bystim);
avg_Fs = nanmean(all_Fs);
% Calculate the trial timeline, convert the idx's into timestamps
% Adding 4 here to account for the timestamps dropped by the alignment and motion correction
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;
figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for stim=stims'
    avg_wt = nanmean(data_bystim.(stim{1}).sub_spectra_wt, 3);
    avg_f = nanmean(data_bystim.(stim{1}).sub_spectra_f, 2);
    nexttile;
    imagesc([0, 3], [min(avg_f) max(avg_f)], abs(avg_wt));
    yticks(flip(avg_f(1:10:length(avg_f))));
    yticklabels(avg_f(1:10:length(avg_f)));
    title(stim, 'Interpreter', 'none');
end
sgtitle('Average of spectra from each individual trace');

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
