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

%% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    for f_stim=stims'
        f_stim = f_stim{1};
        
        popul_data = data_bystim.(f_stim);

        % Grab only modulated neurons
        nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);

        % Setup figure plot
        figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
        tiledlayout(length(nr_idxs), 1, 'TileSpacing', 'none', 'Padding', 'none');
        %tiledlayout(length(nr_idxs), 1, 'TileSpacing', 'compact', 'Padding', 'compact',...
        %    'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
        
        % Loop through each neuron
        axes = [];
        for nr=nr_idxs'
            nr_stim_time = popul_data.stim_timestamps(:, nr);
            trace_time = popul_data.trace_timestamps(:, nr);
            nr_rasters =popul_data.all_trial_spike_rasters{nr}; 
            
            temp = nexttile;
            axes(end + 1) = temp;

            % Plot all of the DBS pulse times
            xline(nr_stim_time);
            hold on;

            plot(trace_time, nr_rasters.*[1:size(nr_rasters, 2)], '.');

            % Set the limits to better see the spike rasters
            xlim([0, 1]);
            ylim([.1, size(nr_rasters, 2)]);
            
        end
        linkaxes(axes);
        % Set the title
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');

    end
end

%% This section will loop through each pulse 
for f_region = fieldnames(region_data)' % {'r_V1'} %
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    for f_stim=stims'
        f_stim = f_stim{1};
        
        popul_data = data_bystim.(f_stim);

        % Grab only modulated neurons
        nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);

        % Setup figure plot
        %figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
        %tiledlayout(length(nr_idxs), 1, 'TileSpacing', 'none', 'Padding', 'none');
        %tiledlayout(length(nr_idxs), 1, 'TileSpacing', 'compact', 'Padding', 'compact',...
        %    'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
        
        all_lats = cell(1, size(popul_data.stim_timestamps, 1) - 1);

        % Loop through each neuron
        figure;
        for nr=nr_idxs'
            nr_stim_time = popul_data.stim_timestamps(:, nr);
            trace_time = popul_data.trace_timestamps(:, nr);
            nr_rasters = popul_data.all_trial_spike_rasters{nr};
            nr_rasters(nr_rasters == 0) = NaN;
            nr_rast_wtime = nr_rasters.*trace_time;

            %Convert spike raster into an array where each spike point is equal to the tracetimepoint
            % Maybe there is a way to sift through that to get the spike time for each spike in that raster
            % Keep the zeroes in this raster to keep window placement
            
            get_time_diff = @(diff_t_mat, pulse_diff) ...
                diff_t_mat(find(diff_t_mat >= 0 & diff_t_mat < pulse_diff));
            
            apply_diff = @(p_idx) get_time_diff(nr_rast_wtime - nr_stim_time(p_idx), ...
                nr_stim_time(p_idx + 1) - nr_stim_time(p_idx) );
            
            % Note this will only get number of stim points - 1 to store
            nr_spike_lats = arrayfun(apply_diff, 1:size(nr_stim_time, 1) - 1, 'UniformOutput', false);
            all_lats = cellfun(@(x, y) [x; y], all_lats, nr_spike_lats, 'UniformOutput', false);
            
            % DEBUG plot spikes for each neuron
            for p_num=1:size(nr_stim_time, 1) - 1
                if ~isempty(nr_spike_lats{p_num})
                    plot(p_num, nr_spike_lats{p_num}*1000, '.');
                    hold on;
                    %pause;
                end
            end
        end
        
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        ylabel('Latency (ms)');
        xlabel('Pulse Number');

       
        %TODO do a per neuron approach with each neuron's trial-averaged firing rate

        % Plot polar plots for both transient and sustained
        % figure;
        % tiledlayout(1, 2);
        % p_width = 1000*mean(diff(nr_stim_time, 1));
        % interv = floor(100/p_width);
        % for p_num=[1, 1+interv]  %1:interv:size(nr_stim_time, 1) - 1
        %     p_num
        %     if p_num+interv > size(nr_stim_time, 1) - 1
        %         interv =size(nr_stim_time, 1) - 1 - p_num;
        %     end
        %     
        %     conct_spikes = cat(1, all_lats{p_num:p_num+interv});
        %     if ~isempty(conct_spikes)
        %         nexttile;
        %         polarhistogram(2*pi*1000*conct_spikes/p_width, 50);
        %         title(num2str(p_num));
        %     end
        % end
        % sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');


        % Group by transient and sustained
        %figure;
        %data = [];
        %labels = [];
        %p_width = 1000*mean(diff(nr_stim_time, 1));
        %interv = floor(100/p_width);
        %for p_num=[1, 1+interv]  %1:interv:size(nr_stim_time, 1) - 1
        %    p_num
        %    if p_num+interv > size(nr_stim_time, 1) - 1
        %        interv =size(nr_stim_time, 1) - 1 - p_num;
        %    end
        %    
        %    conct_spikes = cat(1, all_lats{p_num:p_num+interv});
        %    if ~isempty(conct_spikes)
        %        data = cat(2, data, conct_spikes');
        %        labels = [labels, repmat({num2str(p_num+interv)}, 1, ...
        %            length(conct_spikes))];
        %    end
        %    interv = size(nr_stim_time, 1);
        %end
        %ViolinOpts = Multi_func.get_default_violin();
        %ViolinOpts.QuartileStyle = 'shadow';
        %group_order = unique(labels, 'stable');
        %violins = violinplot(data, labels, 'GroupOrder', group_order, ViolinOpts);



        % Group by time block
        figure;
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        ylabel('Latency (ms)');
        xlabel('Pulse Number');
        %interv = floor(1000*mean(diff(nr_stim_time, 1))/2);
        interv = floor(size(nr_stim_time, 1)/10);
        %interv = 40;
        data = [];
        labels = [];
        for p_num=1:interv:size(nr_stim_time, 1) - 1
            if p_num+interv > size(nr_stim_time, 1) - 1
                interv =size(nr_stim_time, 1) - 1 - p_num;
            end

            conct_spikes = cat(1, all_lats{p_num:p_num+interv});
            if ~isempty(conct_spikes)
                %errorbar(mean(p_num + interv), mean(conct_spikes*1000), ...
                %    std(conct_spikes*1000));
                %hold on;
                %plot(mean(p_num + interv), conct_spikes*1000, '.');
                %hold on;

                data = cat(2, data, conct_spikes');
                labels = [labels, repmat({[num2str(p_num) '-' num2str(p_num + interv)]}, 1, ...
                    length(conct_spikes))];
            end
        end
        ViolinOpts = Multi_func.get_default_violin();
        ViolinOpts.QuartileStyle = 'shadow';
        group_order = unique(labels, 'stable');
        violins = violinplot(data, labels, 'GroupOrder', group_order, ViolinOpts);


        % Plot each latency point without the box plots
        %figure;
        %for p_num=1:size(nr_stim_time, 1) - 1
        %    if ~isempty(all_lats{p_num})
        %        plot(p_num, all_lats{p_num}*1000, '.');
        %        hold on;
        %    end
        %end

        % Plot the timelines
        % figure;
        % data = [cat(1, all_lats{:})].*1000; % Converted to ms
        % % Constructing labels variable
        % labels = arrayfun(@(x) repmat({num2str(x)}, 1, length(all_lats{x})), ...
        %     1:size(nr_stim_time, 1) - 1, 'UniformOutput', false);

        % labels = [cat(2, labels{:})];
        % ViolinOpts = Multi_func.get_default_violin();
        % group_order = arrayfun(@(x) num2str(x), 1:size(nr_stim_time, 1) - 1, 'UniformOutput', false);
        % violins = violinplot(data, labels, 'GroupOrder', group_order, ViolinOpts);
        
        % ylabel('Latency (ms)');
        % xlabel('Pulse Number');
        % % Set the title
        % sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');

    end
end


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
                    cur_spike_idx = trial_data.spike_info375.spike_idx{1};
                    cur_spike_idx(cur_spike_idx < front_frame_drop | cur_spike_idx > back_frame_drop) = [];
                    cur_spike_idx = cur_spike_idx - front_frame_drop;

                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{5}));
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop);

                    % Get the spike times from the first pulse
                    spike_times = cur_trace_time(cur_spike_idx);
                    spike_times(spike_times <= cur_stim_time(1) | spike_times >= cur_stim_time(2)) = [];
                    
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
                    spike_times(spike_times <= cur_stim_time(1) | spike_times >= cur_stim_time(end)) = [];
                    
                    % Skip if there are no spikes during stim period
                    if isempty(spike_times)
                        continue;
                    end

                    %loop through each spike time and get the minimum time from pulse
                    for spike_t = spike_times'
                        pulse_to_spike_time = spike_t - cur_stim_time;
                        pulse_to_spike_time(pulse_to_spike_time <= 0) = [];
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

%% Plot the histograms for all pulse to spike times
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 500]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    histogram(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    histogram(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    avg_first_pts = 1000*data_bystim.(f_stim{1}).all_first_pulse_spike_times(:);
    xline(mean(avg_first_pts, 'omitnan'), 'o' );
    hold on;
    avg_all_pts = 1000*data_bystim.(f_stim{1}).all_all_pulse_spike_times(:);
    xline(mean(avg_all_pts, 'omitnan'), 'b' );
    xlabel('Time (ms)');
    legend('All Pulses', 'First Pulses', 'Avg first', 'Avg all');
    title([f_stim{1}(3:end)], 'Interpreter', 'none');
    
    % Print the descriptive statistics for the spike to pulse stuff
    disp(['Printing stats for ' f_stim{1}]);
    disp(['Mean±std: ' num2str(mean(avg_all_pts, 'omitnan')) '±' num2str(std(avg_all_pts, 'omitnan'))]);
    disp('');
end
sgtitle('Spike to Pulse Times');
saveas(gcf, [figure_path 'Inter_Spike' f 'Hist_Spike_to_pulse.png']);
saveas(gcf, [figure_path 'Inter_Spike' f 'Hist_Spike_to_pulse.eps'], 'epsc');

%% Plot the violing plots pulse to spike times
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 500]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    data = [data_bystim.(f_stim{1}).all_all_pulse_spike_times(:)*1000; data_bystim.(f_stim{1}).all_first_pulse_spike_times(:)*1000];
    num_all_pulses = length(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:));
    num_first_pulses = length(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:));
    labels = [repmat({'All pulses'}, num_all_pulses, 1); repmat({'First pulses'}, num_first_pulses, 1)];
    violinplot2(data, labels);
    ylabel('Time (ms)');
    title([f_stim{1}(3:end)], 'Interpreter', 'none');
    set(gca, 'color', 'none');

    % Print the descriptive statistics for the spike to pulse stuff
    disp(['Printing stats for ' f_stim{1}]);
    disp(['Mean±std: ' num2str(mean(avg_all_pts, 'omitnan')) '±' num2str(std(avg_all_pts, 'omitnan'))]);
    disp('');
end
sgtitle('Spike to Pulse Times');
saveas(gcf, [figure_path 'Inter_Spike' f 'Violin_Spike_to_pulse.png']);
saveas(gcf, [figure_path 'Inter_Spike' f 'Violin_Spike_to_pulse.eps'], 'epsc');

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
