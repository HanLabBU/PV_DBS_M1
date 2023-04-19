clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate
srate_win = 100;

% Parameter to determine whether to combine all regions as one data
all_regions = 1;

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);

all_matfiles = {ses.name};

% Select matfiles by brain region
[region_matfiles] = Multi_func.find_region(all_matfiles);
region_data = struct();
all_Fs = [];
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    %% Select matfiles by stim specific conditions for all regions
    %[matfile_stim] = stim_cond(all_matfiles); 
    %% Select matfiles by stim condition for given region
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);

    %% Loop through each field of the struct and concatenate everything together
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies

    % Loop through each stimulation condition
    for f_stim = fieldnames(matfile_stim)'
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    
    
        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();
        data_bystim.(f_stim).neuron_Vm = [];
        data_bystim.(f_stim).neuron_srate = [];
        data_bystim.(f_stim).mouse_speed = [];
        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];
        data_bystim.(f_stim).speed_timestamps = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_subVm = [];
            cur_fov_srate = [];
            cur_fov_raster = [];
            cur_fov_speed = [];
            cur_fov_stim_time = [];
            cur_fov_trace_time = [];
            cur_fov_speed_time = [];

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
                    if isempty(trial_data) || sum(isnan(cur_fov_subVm(:))) > 0
                        continue;
                    end
                    % Check if trial contains movement data
                    if ~isfield(trial_data, 'speed')
                        continue;
                    end

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;

                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');
                    
                    % Calculate the spike rate
                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win, 1, 1);
                    cur_fov_srate = horzcat_pad(cur_fov_srate, cur_spikerate');            
                    
                    % Store the raster plot
                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');

                    % Grab movement data
                    cur_fov_speed = horzcat_pad(cur_fov_speed, trial_data.speed);

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
                    cur_fov_speed_time = horzcat_pad(cur_fov_speed_time, trial_data.speed_timestamp - stim_start);               
    
                    % DEBUG
                    if max(cur_fov_trace_time(:, end)) < 2
                        stim_trace_sz = size(cur_fov_stim_time)
                        trace_time_sz = size(cur_fov_trace_time)
                        pause;
                        matfile{1}
                    end
                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Store the speed data
            temp = data_bystim.(f_stim).mouse_speed;
            data_bystim.(f_stim).mouse_speed = horzcat_pad(temp, nanmean(cur_fov_speed, 2));

            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_Vm;
            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
            % Store average spike rate for each neuron
            temp = data_bystim.(f_stim).neuron_srate;
            data_bystim.(f_stim).neuron_srate = horzcat_pad(temp, nanmean(cur_fov_srate, 2));
            % Store the timestamp data
            temp = data_bystim.(f_stim).stim_timestamps;
            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
            temp = data_bystim.(f_stim).trace_timestamps;
            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
            temp = data_bystim.(f_stim).speed_timestamps;
            data_bystim.(f_stim).speed_timestamps = horzcat_pad(temp, nanmean(cur_fov_speed_time, 2));
 
        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

% Check if combining all of the regions or not
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%% Full collective spike rate over time
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);
    tiledlayout(2, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim{1}).neuron_srate, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile(find(strcmp(f_stim, stims) == 1));
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        set(gca, 'color', 'none');
        xlabel('Time from stim onset (S)');
        ylabel('Firing Rate (Hz)');
        legend(['Num FOVS ' num2str(num_neurons)])
        title(f_stim{1}(3:end), 'Interpreter', 'none');

        % Plot movement data
        nexttile(find(strcmp(f_stim, stims) == 1) + (length(stims)));
        speed_timeline = mean(data_bystim.(f_stim{1}).speed_timestamps, 2, 'omitnan');
        cur_speed = mean(data_bystim.(f_stim{1}).mouse_speed, 2, 'omitnan');
        std_speed = std(data_bystim.(f_stim{1}).mouse_speed, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).mouse_speed, 2);
        sem_speed = std_speed./sqrt(num_neurons);
        f = fill([speed_timeline; flip(speed_timeline)], [cur_speed + sem_speed; flipud(cur_speed - sem_speed)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(speed_timeline, cur_speed, 'k', 'linewidth', 1);
        set(gca, 'color', 'none');
        xlabel('time from stim onset (s)');
        ylabel('speed (cm/s)');
        title(['movement ' f_stim{1}(3:end)], 'interpreter', 'none');
    end
    sgtitle([f_region ' Average Spike rate with movement'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Movement/' f_region '_Movement_Continuous_FiringRate.png']);
    saveas(gcf, [figure_path 'Movement/' f_region '_Movement_Continuous_FiringRate.eps'], 'epsc');
end

%% Subthreshold Vm
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);
    tiledlayout(2, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_Vm = mean(data_bystim.(f_stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile(find(strcmp(f_stim, stims) == 1));
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        title(f_stim{1}(3:end), 'Interpreter', 'none');
        xlabel('Time from stim onset (S)');
        ylabel('Raw Vm (A.U.)');
        set(gca, 'color', 'none');

        % Plot movement data
        nexttile(find(strcmp(f_stim, stims) == 1) + (length(stims)));
        speed_timeline = mean(data_bystim.(f_stim{1}).speed_timestamps, 2, 'omitnan');
        cur_speed = mean(data_bystim.(f_stim{1}).mouse_speed, 2, 'omitnan');
        std_speed = std(data_bystim.(f_stim{1}).mouse_speed, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).mouse_speed, 2);
        sem_speed = std_speed./sqrt(num_neurons);
        f = fill([speed_timeline; flip(speed_timeline)], [cur_speed + sem_speed; flipud(cur_speed - sem_speed)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(speed_timeline, cur_speed, 'k', 'linewidth', 1);
        set(gca, 'color', 'none');
        xlabel('time from stim onset (s)');
        ylabel('speed (cm/s)');
        title(['movement ' f_stim{1}(3:end)], 'interpreter', 'none');
 
    end
    sgtitle([f_region ' Average subthreshold Vm with movement'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Movement/' f_region '_Movement_sub_thres.png']);
    saveas(gcf, [figure_path 'Movement/' f_region '_Movement_sub_thres.eps']);
end

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

