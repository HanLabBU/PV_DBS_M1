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

% Loop through each brain region

    % Loop through each stimulation condition
    for field = fieldnames(matfile_stim)'
        field = field{1};
        matfiles = matfile_stim.(field).names;    
    
        % Initialize field subthreshold array
        data_bystim.(field) = struct();
        data_bystim.(field).neuron_Vm = [];
        
    
        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
    
            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                cur_fov_subVm = [];
                
                % Store the camera framerate
                all_Fs(end+1) = trial_data.camera_framerate;
                
                % Loop through each trial                
                for tr_idx=1:length(data.align.trial)        
    
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
                    
                    trial_data = data.align.trial{tr_idx};
                    % If the trial data is empty, that means it was skipped
                    if isempty(trial_data)
                        continue;
                    end
    
                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, cur_trace_ws');
                    
                    % Keep track of the figure being used
                    if strcmp(ri{5}, '1000') == 1
                        figure('visible', 'off');
                        plot(cur_trace_ws)
                        title([matfile{1} ' ROI#' num2str(roi_idx) ' Trial ' num2str(tr_idx)], 'Interpreter', 'none');
                        saveas(gcf, ['debug/' matfile{1}(1:end-4) '.png']);
                        close gcf;
                    end
                end
            end
            
            %EBUG
            %figure;
            %plot(nanmean(cur_fov_subVm, 2));
    
            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(field).neuron_Vm;
            data_bystim.(field).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
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
    nexttile;
    plot(timeline, nanmean(data_bystim.(stim{1}).neuron_Vm, 2));
    title(stim, 'Interpreter', 'none');
end
sgtitle('Average Sub Vm by Stimulation condition', 'Interpreter', 'none');

% Quick FT check
figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
freqLimits = [0 90];
for stim=stims'
    if strcmp(stim{1}, 'f_1000') == 1
        freqLimits = [0 1000];
    end
    cur_subVm = nanmean(data_bystim.(stim{1}).neuron_Vm, 2);
    nexttile;
    fb = cwtfilterbank(SignalLength=length(cur_subVm),...
                       SamplingFrequency=avg_Fs,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(cur_subVm, FilterBank=fb);
    contourf(timeline, f, abs(wt), 'edgecolor', 'none');
    colorbar;
    if strcmp(stim{1}, 'f_40') == 1
        yline(40);
        hold on;
        yline(80);
    end
    title(stim, 'Interpreter', 'none');
    %nexttile;
    %imagesc(timeline, f, abs(wt));
    %
    %%TODO change scale for showing the frequencies
    %yticks(flip(f([1, 12, 24, 64])));
    %yticklabels(string(flip(f([1, 12, 24, 64] ))));
    %title(stim, 'Interpreter', 'none');
end
sgtitle('Spectra from averaged Sub Vm Trace');

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
