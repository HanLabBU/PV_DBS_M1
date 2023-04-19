clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine how much before and after a stimulation pulse to take the average
trace_sur = 10; % This is ~6ms before and after stim pulse

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
        data_bystim.(f_stim).all_first_pulse_spike_times = [];
        data_bystim.(f_stim).all_all_pulse_spike_times = [];
        data_bystim.(f_stim).neuron_avg_first_pulse_spike_times = [];
        data_bystim.(f_stim).neuron_avg_all_pulse_spike_times = [];


        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_first_pulse_spike_times = [];
            cur_fov_all_pulse_spike_times = [];

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

                    % For now, I will not have this criteria
                    % Only perform calculation for trials that have more than 5 spikes
                    %if length(trial_data.spike_info375.spike_idx{1}) < 5
                    %    continue;
                    %end

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;

                    % Grab the trace, raster, spike idxs, and timestamps
                    cur_trace = trial_data.detrend_traces(front_frame_drop:back_frame_drop, roi_idx);
                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spike_idx = trial_data.spike_info375.spike_idx{1};
                    cur_spike_idx(cur_spike_idx < front_frame_drop | cur_spike_idx > back_frame_drop) = [];
                    cur_spike_idx = cur_spike_idx - front_frame_drop;

                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{5}));
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop);

                    % Get the spike times from the first pulse
                    spike_times = cur_trace_time(cur_spike_idx);
                    spike_times(spike_times < cur_stim_time(1) | spike_times > cur_stim_time(2)) = [];
                    
                    first_to_pulse_time = spike_times - cur_stim_time(1);
                    first_to_pulse_time = first_to_pulse_time(:)';
                    
                    %DEBUG
                    if first_to_pulse_time >= diff([cur_stim_time(1), cur_stim_time(2)])
                        disp('Spike time to large!!');
                        pause;
                    end

                    cur_fov_first_pulse_spike_times = horzcat_pad(cur_fov_first_pulse_spike_times, first_to_pulse_time);

                    % Get the spike times from its closest preceding pulse
                    spike_times = cur_trace_time(cur_spike_idx);
                    spike_times(spike_times < cur_stim_time(1) | spike_times > cur_stim_time(end)) = [];
                    
                    % Skip if there are no spikes during stim period
                    if isempty(spike_times)
                        continue;
                    end

                    %loop through each spike time and get the minimum time from pulse
                    for spike_t = spike_times'
                        pulse_to_spike_time = spike_t - cur_stim_time;
                        pulse_to_spike_time(pulse_to_spike_time < 0) = [];
                        pulse_to_spike_time = min(pulse_to_spike_time);
                        cur_fov_all_pulse_spike_times = horzcat_pad(cur_fov_all_pulse_spike_times, pulse_to_spike_time);
                    end                  

                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_all_pulse_spike_times(:))) || isempty(cur_fov_all_pulse_spike_times)
                continue;
            end

            % Plot the average and all of the trace from the stim centers
            %figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 500 1000]);
            %histogram(cur_fov_base_inter);
            %hold on;
            %histogram(cur_fov_stim_inter);
            %legend('Base', 'Stim');
            %title([matfile{1}], 'Interpreter', 'none');
            %saveas(gcf, [figure_path 'Inter_Spike' f matfile{1}(1:end-4) '_inter_spike.png']);

            % Save each FOV pulse to spike times
            temp = data_bystim.(f_stim).all_first_pulse_spike_times;
            data_bystim.(f_stim).all_first_pulse_spike_times = horzcat_pad(temp, cur_fov_first_pulse_spike_times(:)');
            temp = data_bystim.(f_stim).all_all_pulse_spike_times;
            data_bystim.(f_stim).all_all_pulse_spike_times = horzcat_pad(temp, cur_fov_all_pulse_spike_times(:)');
            temp = data_bystim.(f_stim).neuron_avg_first_pulse_spike_times;
            data_bystim.(f_stim).neuron_avg_first_pulse_spike_times = horzcat_pad(temp, mean(cur_fov_first_pulse_spike_times(:), 'omitnan'));
            temp = data_bystim.(f_stim).neuron_avg_all_pulse_spike_times;
            data_bystim.(f_stim).neuron_avg_all_pulse_spike_times = horzcat_pad(temp, mean(cur_fov_all_pulse_spike_times(:), 'omitnan'));
            
            %DEBUG
            if 0 && length(cur_fov_first_pulse_spike_times) > 0
                figure;
                tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
                nexttile;
                plot(cur_stim_time, 2, 'b|');
                hold on;
                plot(spike_times, 1, 'r.');
                
                nexttile;
                plot(1, first_to_pulse_time, '.');
                sgtitle('First pulse to spike times');
                %pause;
            end

        end % End looping through FOVs of a condition
    end

%    % Save the VM to the specific region
%    region_data.(f_region) = data_bystim;

%end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%% Plot the histograms for all spike-intervals
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 500 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    histogram(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    histogram(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    avg_first_pts = 1000*mean(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:), 'omitnan');
    xline(avg_first_pts, 'o' );
    hold on;
    avg_all_pts = 1000*mean(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:), 'omitnan');
    xline(avg_all_pts, 'b' );
    xlabel('Time (ms)');
    legend('All Pulses', 'First Pulses', ['Avg first ' num2str(avg_first_pts)], ['Avg all ' num2str(avg_all_pts)]);
    title([f_stim{1}(3:end)], 'Interpreter', 'none');
end
sgtitle('Spike to Pulse Times');
%saveas(gcf, [figure_path 'Inter_Spike' f matfile{1}(1:end-4) '_all_interspike.png']);
%saveas(gcf, [figure_path 'Inter_Spike' f matfile{1}(1:end-4) '_all_interspike.eps']);

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
