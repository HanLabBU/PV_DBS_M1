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
% TODO could use a smaller window size
srate_win = 50;

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
        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];

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
            cur_fov_stim_time = [];
            cur_fov_trace_time = [];

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

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, ...
                                                    raw_trial_data.raw_stimulation_time(1:str2num(ri{5}) ) - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, ...
                                                    trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
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

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;


%% Spike rate aligned by the first stimulation pulse
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)*1000;
        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim{1}).neuron_srate, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
        sem_srate = cur_srate./sqrt(num_neurons);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        hold on;

        % Plot the DBS stimulation time pulses
        xline(nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        xlim([0 - 50, 0 + 50]);
        set(gca, 'color', 'none');
        ylabel('Firing Rate(Hz)');
        xlabel('Time from onset(ms)');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Population Spike rate with first pulse'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Trig_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Trig_FR.eps'], 'epsc');
end

%% Sub Vm averaged across all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        
        Vm_avg = [];
        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim{1}).neuron_Vm, 2)
            
            % TODO may be a good idea to just have a fixed width for all stimulations, independent of frequency

            % Calculate average DBS pulse time widths
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim{1}).stim_timestamps(:, nr) ), 'omitnan');

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim{1}).stim_timestamps(:, nr)'
                start_trace_idx = find(pulse_time <= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                start_trace_idx = start_trace_idx(1);
                end_trace_idx = find(pulse_time + nr_avg_pulse_width >= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                end_trace_idx = end_trace_idx(end);
                
                fr_pulse_width = data_bystim.(f_stim{1}).neuron_Vm(start_trace_idx:end_trace_idx, nr);
                Vm_avg = horzcat_pad(Vm_avg, fr_pulse_width);
            end
        end

        cur_subVm = mean(Vm_avg, 2, 'omitnan');
        std_subVm = std(Vm_avg, 0, 2, 'omitnan');
        num_pulses = size(Vm_avg, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_subVm, 1);
        sem_subVm = cur_subVm./sqrt(num_pulses);
        nexttile;

        %TODO the timeline and the stimulation pulses are not calculated by the timestamps
        f = fill([1:size(Vm_avg, 1), size(Vm_avg, 1):-1:1], [cur_subVm + sem_subVm; flipud(cur_subVm - sem_subVm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(1:size(Vm_avg, 1), cur_subVm, 'k', 'LineWidth', 1);
        hold on;

        % Plot the DBS stimulation time pulses
        xline([1, end_trace_idx - start_trace_idx], 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = 1;
        posy = 0;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.2, '50ms');

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        axis off;
        set(gca, 'color', 'none');
        ylabel('Vm');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Sub Vm all pulse average'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.eps'], 'epsc');
end

%TODO need to finish implementing this
%% Spike rate averaged across all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        
        FR_avg = [];
        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim{1}).neuron_srate, 2)
            
            % TODO may be a good idea to just have a fixed width for all stimulations, independent of frequency

            % Calculate average DBS pulse time widths
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim{1}).stim_timestamps(:, nr) ), 'omitnan');

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim{1}).stim_timestamps(:, nr)'
                start_trace_idx = find(pulse_time <= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                start_trace_idx = start_trace_idx(1);
                end_trace_idx = find(pulse_time + nr_avg_pulse_width >= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                end_trace_idx = end_trace_idx(end);
                
                fr_pulse_width = data_bystim.(f_stim{1}).neuron_srate(start_trace_idx:end_trace_idx, nr);
                FR_avg = horzcat_pad(FR_avg, fr_pulse_width);
            end
        end

        cur_srate = mean(FR_avg, 2, 'omitnan');
        std_srate = std(FR_avg, 0, 2, 'omitnan');
        num_pulses = size(FR_avg, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
        sem_srate = cur_srate./sqrt(num_pulses);
        nexttile;
        f = fill([1:size(FR_avg, 1), size(FR_avg, 1):-1:1], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(1:size(FR_avg, 1), cur_srate, 'k', 'LineWidth', 1);
        hold on;

        % Plot the DBS stimulation time pulses
        xline([1, end_trace_idx - start_trace_idx], 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = 1;
        posy = 0;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.2, '50ms');

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        axis off;
        set(gca, 'color', 'none');
        ylabel('Firing Rate(Hz)');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Firing Rate all pulse average'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.eps'], 'epsc');
end


%% Spike rate showing display all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim{1}).neuron_srate, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
        sem_srate = cur_srate./sqrt(num_neurons);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        hold on;

        % Plot the DBS stimulation time pulses
        stim_time = nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2);
        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = -.100;
        posy = 0;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.2, '50ms');

        % Increase timescale resolution
        xlim([0 - .10, max(stim_time) + .10]);
        a = gca;
        a.XAxis.Visible = 'off';
        set(gca, 'color', 'none');
        ylabel('Firing Rate(Hz)');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Population Spike rate with all pulse'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.eps'], 'epsc');
end


%% Subthreshold Vm first few pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)*1000; % Convert to ms
        cur_Vm = mean(data_bystim.(f_stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;

        % Standard Error
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        
        % Plot Vm
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        hold on;
        
        % Plot the DBS stimulation time pulses
        xline(nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        xlim([0 - 50, 0 + 50]);
        ylabel('Vm');
        xlabel('Time from onset(ms)');
        set(gca, 'color', 'none')
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Trig_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Trig_Vm.eps'], 'epsc');
end

%% Subthreshold Vm first all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_Vm = mean(data_bystim.(f_stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;

        % Standard Error
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        
        % Plot Vm
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        hold on;
        
        % Plot the DBS stimulation time pulses
        stim_time = nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2);
        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = -.100;
        posy = -3;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.5, '50ms');

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]);
        a = gca;
        a.XAxis.Visible = 'off';
        ylabel('Vm');
        set(gca, 'color', 'none')
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_Vm.eps'], 'epsc');
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
