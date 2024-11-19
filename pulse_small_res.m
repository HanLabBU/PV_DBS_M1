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

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate

% Do all region
all_region = 0;

%%% END Modification

% Seed number for random number generation
rng(100);

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Load all of the data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
%save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'ca1_data.mat'];

load(save_all_data_file);

% Combine all regions into 1
if all_region == 1
    region_data = Multi_func.combine_regions(region_data);
end

% Get the first region field
field1 = fieldnames(region_data);
field1 = 'r_M1';
avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

% Determine how much of extra Vm to use in the all pulse averaging
extra_trace = 3;

%% Spike rate showing first few pulses

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
nr_pop = 'all';
%nr_pop = 'etrain';
%nr_pop = 'non';

Fr_onset_time = struct();
stats_log = [figure_path 'Small_Res' f 'Fr_onset_sig_times_' nr_pop];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    for f_stim=stims'
        f_stim = f_stim{1};

        popul_data = data_bystim.(f_stim);

        % Check if there is an entrained field in the population data
        try
            switch nr_pop
                case 'etrain'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] > 0);
                case 'non'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] < 0);
                case 'all'
                    nr_idxs = 1:length(popul_data.plv_mod_stats);
            end
        catch ME
            disp(ME.message);
        end
        
        timeline = nanmean(popul_data.trace_timestamps, 2)*1000;
        cur_srate = mean(popul_data.neuron_srate_3(:, nr_idxs), 2, 'omitnan');
        std_srate = std(popul_data.neuron_srate_3(:, nr_idxs), 0, 2, 'omitnan');
        num_neurons = size(popul_data.neuron_srate_3(:, nr_idxs), 2);
        %num_points = size(popul_data.neuron_srate_3, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile;
        fill_h = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        hold on;

        % Plot the DBS stimulation time pulses
        xline(nanmean(popul_data.stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the semi-significant points on the subthreshold Vm
        base_srate = [];
        % Grab the baseline sub srate for all neurons
        for i = 1:size(popul_data.trace_timestamps, 2)
            baseline_idx = find(popul_data.trace_timestamps(:, i) < popul_data.stim_timestamps(1, i));
            base_srate = horzcat_pad(base_srate, popul_data.neuron_srate_3(baseline_idx, i) );
        end
        base_srate = mean(base_srate, 2, 'omitnan');
        std_baseline = std(base_srate, 0, 'omitnan');
        sig_idx = find(cur_srate > (2*std_baseline + mean(base_srate, 'omitnan')));

        % Calculate a good height for the points
        thres = 0.2;
        height = range(cur_srate);
        height = (1+thres)*height;
        %plot(timeline(sig_idx), repmat(height, 1, length(sig_idx)), '.b', 'MarkerSize', 10);
        
        % Plotting the 2 std line 
        yline(2*std_baseline, '--', 'Color', [0 0 0 0.5]);
        
        % Remove significant points that are during the baseline period
        sig_idx(find(timeline(sig_idx) < 0)) = [];
        % Plot and calculate the time to first significance
        if ~isempty(sig_idx)
            %DEBUG show where the first significant point is
            xline(timeline(sig_idx(1)), 'g');

            diary on;
            fprintf('\n\n');
            disp([f_region ' ' f_stim]);
            disp(['Sig FR first pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
            fprintf('\n\n');
            diary off;
        end

        % Increase timescale resolution
        xlim([0 - 50, 0 + 100]);
        x = gca; x = x.XAxis;
        Multi_func.set_spacing_axis(x, 50, 1);
        Multi_func.set_default_axis(gca);
        ylabel('Firing Rate(Hz)');
        xlabel('Time from onset(ms)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Population Spike rate with first pulse'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_First_Pulse_FR.png']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_First_Pulse_FR.pdf']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_First_Pulse_FR.eps'], 'epsc');
end

%% All regions overlay first few pulses Firing Rate
test_regions = fieldnames(region_data)' %{'r_CA1', 'r_M1'};
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    nexttile;

    %Loop through regions
    for f_region = test_regions
        f_region = f_region{1};
        stim_data = region_data.(f_region).(f_stim);

        timeline = nanmean(stim_data.trace_timestamps, 2)*1000;
        cur_srate = mean(stim_data.neuron_srate_3, 2, 'omitnan');
        std_srate = std(stim_data.neuron_srate_3, 0, 2, 'omitnan');
        num_neurons = size(stim_data.neuron_srate_3, 2);
        %num_points = size(stim_data.neuron_srate_3, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        fill_h = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5], 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        hold on;
        
        %TODO add V1 colors
        if strcmp(f_region, 'r_CA1') == 1
            cur_color = Multi_func.ca1_color;
            label = 'CA1';
        elseif strcmp(f_region, 'r_M1') == 1
            cur_color = Multi_func.m1_color;
            label = 'M1';
        elseif strcmp(f_region, 'r_V1') == 1
            cur_color = Multi_func.v1_color;
            label = 'V1';
        end

        plot(timeline, cur_srate, 'Color', cur_color, 'LineWidth', 1, 'DisplayName', label);
        hold on;

        % Plot the DBS stimulation time pulses
        xline(nanmean(stim_data.stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5, 'HandleVisibility', 'off');
        hold on;
    end

    % Increase timescale resolution
    xlim([0 - 50, 0 + 100]);
    x = gca; x = x.XAxis;
    Multi_func.set_spacing_axis(x, 50, 1);
    Multi_func.set_default_axis(gca);
    legend('Location', 'northeast');

    ylabel('Firing Rate (Hz)');
    xlabel('Time from onset(ms)');
    title([f_stim(3:end) ' Hz'], 'Interpreter', 'none');
end
sgtitle('Onset Fr Plot by Region Overlay');
saveas(gcf, [figure_path 'Small_Res' f 'Onset_Overlay_Fr.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Onset_Overlay_Fr.pdf']);

%% All Vm pulse-triggered averaged across all pulses
% This calculates a shuffled distribution 

extra_trace = 3;
num_iter = 1;

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non_entr';

% Flag for which statistical method to use
%stat_met = 'shuff';
stat_met = 'sign';

% Flag for removing non-modulated neurons from analysis
remove_nonmod_nrs = 1;

vm_trig_avg_time = struct();
stats_log = [figure_path 'Small_Res' f 'Vm_pulse_triggered_time_to_final_average_' nr_pop];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Set up all of the figures
    all_pulse_fig = figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 30.59 27.94]);
    %all_tt = tiledlayout(all_pulse_fig, length(stims), 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3*3.62, 5.16]); % Use for fourth of a figure
    %all_tt = tiledlayout(all_pulse_fig, length(stims), 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3*2.27, 5.16]); % Use for a sixth of a figure
    all_tt = tiledlayout(all_pulse_fig, length(stims), 3, 'TileSpacing', 'compact', 'Padding', 'compact'); % Use to help with zooming in
    
    tilenum = 1;

    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Check if there is an entrained field in the population data
        try
            switch nr_pop
                case 'etrain'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] > 0);
                case 'non_entr'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] < 0);
                case 'all'
                    nr_idxs = 1:length(popul_data.plv_mod_stats);
            end
        catch ME
            disp(ME.message);
        end
        
        % Filter out the neurons that were non-modulated at all
        if remove_nonmod_nrs == 1
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_idxs = nr_idxs(~ismember(nr_idxs, non_mod_nr));
        end

        num_neurons = length(nr_idxs);

        all_pulse_vm_by_neuron = cell(num_neurons, 1);
        all_pulse_vm = [];
        norm_vms = popul_data.neuron_RawVm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        all_pulse_trace_width = [];
        all_pulse_width_time = [];

        % Save the pulse data for transient and sustained period  
        trans_pulse_vm = [];
        trans_pulse_vm_by_neuron = cell(num_neurons, 1);
        sus_pulse_vm = [];
        sus_pulse_vm_by_neuron = cell(num_neurons, 1); 

        % Looping through each neuron
        for nr = nr_idxs
            
            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width_time = mean(diff(popul_data.stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(popul_data.trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width_time./nr_avg_trace_time);
        
            all_pulse_trace_width = horzcat_pad(all_pulse_trace_width, trace_width(:));
            all_pulse_width_time = horzcat_pad(all_pulse_width_time, nr_avg_pulse_width_time(:));

            % Use fixed length of points for all frequencies
            %trace_width = 21;

            % Loop through each stimulation time pulses
            for pulse_time = popul_data.stim_timestamps(:, nr)'
                follow_trace_idx = find(pulse_time <= popul_data.trace_timestamps(:, nr));
                
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                
                % Calculate end trace idx using the average inter-pulse time
                %end_trace_idx = find(pulse_time + trace_width >= data_bystim.(f_stim).trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                % Using fixed number of pulse times
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;

                % "Baseline" subtract from the first few trace points
                %Vm_pulse_width = norm_vms(start_trace_idx:end_trace_idx, nr) - nanmean(norm_vms(start_trace_idx:start_trace_idx + extra_trace, nr));

                % store the pulse width
                Vm_pulse_width = norm_vms(start_trace_idx:end_trace_idx, find(nr == nr_idxs));

                % Store all of the Vm pulses per neuron preserving relative stimulation period shapes
                all_pulse_vm_by_neuron{find(nr == nr_idxs)} = horzcat_pad(all_pulse_vm_by_neuron{find(nr == nr_idxs)}, Vm_pulse_width');

                % Store windows depending on period of pulse
                %Transient pulse
                if pulse_time >= Multi_func.trans_ped(1)/1000 & pulse_time < Multi_func.trans_ped(2)/1000
                    trans_pulse_vm_by_neuron{find(nr == nr_idxs)} = horzcat_pad(trans_pulse_vm_by_neuron{find(nr == nr_idxs)},...
                        Vm_pulse_width');
                end

                %Sustained pulse
                if pulse_time >= Multi_func.sus_ped(1)/1000 & pulse_time < Multi_func.sus_ped(2)/1000
                    sus_pulse_vm_by_neuron{find(nr == nr_idxs)} = horzcat_pad(sus_pulse_vm_by_neuron{find(nr == nr_idxs)},...
                        Vm_pulse_width');
                end

                % Subtract from the pulse onset
                Vm_pulse_width_on_sub = norm_vms(start_trace_idx:end_trace_idx, find(nr == nr_idxs)) ...
                    - norm_vms(start_trace_idx + extra_trace, find(nr == nr_idxs));

                % Store all of the onset subtracted Vm pulses in one array
                all_pulse_vm = horzcat_pad(all_pulse_vm, Vm_pulse_width_on_sub);

                % Store onset SUBTRACTED windows depending on period of pulse
                %Transient pulse
                if pulse_time >= Multi_func.trans_ped(1)/1000 & pulse_time < Multi_func.trans_ped(2)/1000
                    trans_pulse_vm = horzcat_pad(trans_pulse_vm,...
                        Vm_pulse_width_on_sub);
                end
                
                %Sustained pulse
                if pulse_time >= Multi_func.sus_ped(1)/1000 & pulse_time < Multi_func.sus_ped(2)/1000
                    sus_pulse_vm = horzcat_pad(sus_pulse_vm,...
                        Vm_pulse_width_on_sub);
                end
            end
        end
        
        % Save all of the pulse triggered Vm 

        % -- Plotting the whole period pulse average
        cur_subVm = mean(all_pulse_vm, 2, 'omitnan');
        std_subVm = std(all_pulse_vm, 0, 2, 'omitnan');
        
        disp([f_region ' ' f_stim]);
        num_pulses = size(all_pulse_vm, 2)
        %num_points = size(data_bystim.(f_stim).neuron_subVm, 1);
        sem_subVm = std_subVm./sqrt(num_pulses);
        
        nexttile(all_tt, tilenum);
        
        % Use the average trace framerate to calculate the timelinefigure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
        timeline = [ [0:size(all_pulse_vm, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        fill_h = fill([timeline, flip(timeline)], [cur_subVm + sem_subVm; flipud(cur_subVm - sem_subVm)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, cur_subVm, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the whole period pulse average
        
        % -- Plotting the transient pulse average
        trans_cur_subVm = mean(trans_pulse_vm, 2, 'omitnan');
        trans_sem_subVm = std(trans_pulse_vm, 0, 2, 'omitnan')./sqrt(size(trans_pulse_vm, 2));
        
        nexttile(all_tt, tilenum + 1);

        timeline = [ [0:size(trans_pulse_vm, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        fill_h = fill([timeline, flip(timeline)], [trans_cur_subVm + trans_sem_subVm; flipud(trans_cur_subVm - trans_sem_subVm)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, trans_cur_subVm, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the transient period pulse average

        % -- Plotting the sustained pulse average
        sus_cur_subVm = mean(sus_pulse_vm, 2, 'omitnan');
        sus_sem_subVm = std(sus_pulse_vm, 0, 2, 'omitnan')./sqrt(size(sus_pulse_vm, 2));
        
        nexttile(all_tt, tilenum + 2);

        timeline = [ [0:size(sus_pulse_vm, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        fill_h = fill([timeline, flip(timeline)], [sus_cur_subVm + sus_sem_subVm; flipud(sus_cur_subVm - sus_sem_subVm)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, sus_cur_subVm, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the transient period pulse average

        switch stat_met
            case 'shuff'
                % This shuffle takes random window samples within the stim period, should be across the pulse times to indicate if the pulse timing is significant
                wind_size = size(all_pulse_vm, 1);
                shuf_val_dist = [];
                shuf_trans_dist = [];
                shuf_sus_dist = [];
                tic;

                for i=1:num_iter
                    iter_wind = [];
                    iter_trans_wind = [];
                    iter_sus_wind = [];

                    for j=1:size(all_pulse_vm, 2) %length(data_bystim.(f_stim).stim_timestamps(:, 1))
                        % Randomly grab a neuron
                        rand_nr = randi([1, length(all_pulse_vm_by_neuron)]);
                        % Grab a window sample that fully fits within the stimulation period
                        rand_win_start = randi([1, length(all_pulse_vm_by_neuron{rand_nr}) - wind_size]);
                        iter_wind = horzcat_pad(iter_wind, ...
                            all_pulse_vm_by_neuron{rand_nr}(rand_win_start:rand_win_start+wind_size - 1)' ...
                            - all_pulse_vm_by_neuron{rand_nr}(rand_win_start + extra_trace));
                            
                        % Same as above, but for both stimulation periods, can use the same neuron by the way
                        trans_rand_win_start = randi([1,...
                            length(trans_pulse_vm_by_neuron{rand_nr}) - wind_size]);
                        iter_trans_wind = horzcat_pad(iter_trans_wind, ...
                            trans_pulse_vm_by_neuron{rand_nr}(trans_rand_win_start:trans_rand_win_start+wind_size - 1)' ...
                            - trans_pulse_vm_by_neuron{rand_nr}(trans_rand_win_start + extra_trace));

                        sus_rand_win_start = randi([1,...
                            length(sus_pulse_vm_by_neuron{rand_nr}) - wind_size]);
                        iter_sus_wind = horzcat_pad(iter_sus_wind, ...
                            sus_pulse_vm_by_neuron{rand_nr}(sus_rand_win_start:sus_rand_win_start+wind_size - 1)' ...
                            - sus_pulse_vm_by_neuron{rand_nr}(sus_rand_win_start + extra_trace));
                    end

                    end_iter = toc;
                    shuf_val_dist = horzcat_pad(shuf_val_dist, mean(iter_wind, 2));
                    shuf_trans_dist = horzcat_pad(shuf_trans_dist, mean(iter_trans_wind, 2));
                    shuf_sus_dist = horzcat_pad(shuf_sus_dist, mean(iter_sus_wind, 2));
                end
                end_shuf = toc
                
                % Save the shuffle distribution data
                region_data.(f_region).(f_stim).vm_shuf_val_dist = shuf_val_dist;
                
                % -- Plotting the all pulse shuffled data
                % Calculate percentile values
                low_perc = prctile(shuf_val_dist, 2.5, 2);
                high_perc = prctile(shuf_val_dist, 97.5, 2);
                shuf_mean = mean(shuf_val_dist, 2);

                nexttile(all_tt, tilenum);
                % Plot the shuffled values as dashed horizontal lines
                plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color); 
                hold on;
                plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
                hold on;

                % Plotting the significant time and peak
                sig_idx = find(cur_subVm > high_perc);
                sig_idx(timeline(sig_idx) <= 0 | timeline(sig_idx) >= 1000*nr_avg_pulse_width_time) = [];

                peak_idx = find(cur_subVm == max(cur_subVm(timeline' >= 0 & timeline' <= nr_avg_pulse_width_time*1000))); % The messsy condition is for inside the pulse window and not the extra traces
                if ~isempty(sig_idx)
                    % DEBUG show the significant point
                    xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
                    
                    diary on;
                    fprintf('\n\n');
                    disp([f_region ' ' f_stim]);
                    disp(['Sig Vm all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
                    disp(['Peak Vm all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
                    fprintf('\n\n');
                    diary off;

                    % Save the peak time for the region and stimulation
                    region_data.(f_region).(f_stim).all_pulse_vm_pk_time = timeline(peak_idx(1));
                end
        
            case 'sign'
                nexttile(all_tt, tilenum);
                
                % Alpha is 0.05, even for two tailed
                adj_alpha = 0.05/size(all_pulse_vm, 1);
                p_vals = arrayfun(@(row) signtest(all_pulse_vm(row, :) ) , 1:size(all_pulse_vm, 1));
                
                sig_idx = find(p_vals < adj_alpha);
                hold on;

                plot(timeline(sig_idx), max(cur_subVm + sem_subVm)*ones(size(sig_idx)), ...
                    'r.', 'MarkerSize', 6)
                hold on;
        end


        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Normalized Vm');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);
        
        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title(f_stim(3:end), 'Interpreter', 'none');
        % -- End plotting the all pulse shuffled data
    
        switch stat_met
            case 'shuff'
                % -- Plotting the trans pulse shuffled data %TODO  double check that I am using the correct trans or sus variable
                % Calculate percentile values
                low_perc = prctile(shuf_trans_dist, 2.5, 2);
                high_perc = prctile(shuf_trans_dist, 97.5, 2);
                shuf_mean = mean(shuf_trans_dist, 2);
                
                nexttile(all_tt, tilenum + 1);
                % Plot the shuffled values as dashed horizontal lines
                plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color); 
                hold on;
                plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
                hold on;

                % Plotting the significant time and peak
                sig_idx = find(sus_cur_subVm > high_perc);
                sig_idx(timeline(sig_idx) <= 0 | timeline(sig_idx) >= 1000*nr_avg_pulse_width_time) = [];

                peak_idx = find(trans_cur_subVm == max(trans_cur_subVm(timeline' >= 0 & timeline' <= nr_avg_pulse_width_time*1000))); % The messsy condition is for inside the pulse window and not the extra traces
                if ~isempty(sig_idx)
                    % DEBUG show the significant point
                    xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
                    
                    diary on;
                    fprintf('\n\n');
                    disp([f_region ' ' f_stim ' for the transient period']);
                    disp(['Sig Vm all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
                    disp(['Peak Vm all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
                    fprintf('\n\n');
                    diary off;

                    % Save the peak time for the region and stimulation
                    region_data.(f_region).(f_stim).trans_pulse_vm_pk_time = timeline(peak_idx(1));
                end
            case 'sign'
                nexttile(all_tt, tilenum + 1);
                
                % Alpha is 0.05, even for two tailed
                adj_alpha = 0.05/size(trans_pulse_vm, 1);
                p_vals = arrayfun(@(row) signtest(trans_pulse_vm(row, :) ) , 1:size(trans_pulse_vm, 1));
                
                sig_idx = find(p_vals < adj_alpha);
                hold on;
                plot(timeline(sig_idx), max(trans_cur_subVm + trans_sem_subVm)*ones(size(sig_idx)), ...
                    'r.', 'MarkerSize', 6)
                hold on;
        end

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Normalized Vm');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);
        
        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title([f_stim(3:end) ' transient'], 'Interpreter', 'none');
        % -- End plotting the trans pulse shuffled data
        
        switch stat_met
            case 'shuff'

                % -- Plotting the sus pulse shuffled data
                % Calculate percentile values
                low_perc = prctile(shuf_sus_dist, 2.5, 2);
                high_perc = prctile(shuf_sus_dist, 97.5, 2);
                shuf_mean = mean(shuf_sus_dist, 2);
                
                nexttile(all_tt, tilenum + 2);
                % Plot the shuffled values as dashed horizontal lines
                plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color); 
                hold on;
                plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
                hold on;

                % Plotting the significant time and peak
                sig_idx = find(sus_cur_subVm > high_perc);
                sig_idx(timeline(sig_idx) <= 0 | timeline(sig_idx) >= 1000*nr_avg_pulse_width_time) = [];

                peak_idx = find(sus_cur_subVm == max(sus_cur_subVm(timeline' >= 0 & timeline' <= nr_avg_pulse_width_time*1000))); % The messsy condition is for inside the pulse window and not the extra traces
                if ~isempty(sig_idx)
                    % DEBUG show the significant point
                    xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
                    
                    diary on;
                    fprintf('\n\n');
                    disp([f_region ' ' f_stim ' for the sustained period']);
                    disp(['Sig Vm all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
                    disp(['Peak Vm all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
                    fprintf('\n\n');
                    diary off;

                    % Save the peak time for the region and stimulation
                    region_data.(f_region).(f_stim).sus_pulse_vm_pk_time = timeline(peak_idx(1));
                end
            case 'sign'
                nexttile(all_tt, tilenum + 2);
                
                % Alpha is 0.05, even for two tailed
                adj_alpha = 0.05/size(sus_pulse_vm, 1);
                p_vals = arrayfun(@(row) signtest(sus_pulse_vm(row, :) ) , 1:size(sus_pulse_vm, 1));
                
                sig_idx = find(p_vals < adj_alpha);
                hold on;
                plot(timeline(sig_idx), max(sus_cur_subVm + sus_sem_subVm)*ones(size(sig_idx)), ...
                    'r.', 'MarkerSize', 6)
                hold on;
        end

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Normalized Vm');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);
        
        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title([f_stim(3:end) ' sustained'], 'Interpreter', 'none');
        % -- End plotting the sus pulse shuffled data

        % Increase tile number indicator
        tilenum = tilenum + 3;

        % Save the pulse and shuffled distribution data
        region_data.(f_region).(f_stim).(nr_pop).all_pulse_trig_Vm = all_pulse_vm;
        region_data.(f_region).(f_stim).(nr_pop).all_pulse_vm_by_neuron = all_pulse_vm_by_neuron;
        
        region_data.(f_region).(f_stim).(nr_pop).trans_pulse_trig_Vm = trans_pulse_vm;
        region_data.(f_region).(f_stim).(nr_pop).trans_pulse_vm_by_neuron = trans_pulse_vm_by_neuron;
        region_data.(f_region).(f_stim).(nr_pop).sus_pulse_trig_Vm = sus_pulse_vm;
        region_data.(f_region).(f_stim).(nr_pop).sus_pulse_vm_by_neuron = sus_pulse_vm_by_neuron;
        region_data.(f_region).(f_stim).(nr_pop).num_neurons = num_neurons; 

        region_data.(f_region).(f_stim).(nr_pop).shuf_all_trig_Vm = shuf_val_dist;
        region_data.(f_region).(f_stim).(nr_pop).shuf_trans_trig_Vm = shuf_trans_dist;
        region_data.(f_region).(f_stim).(nr_pop).shuf_sus_trig_Vm = shuf_sus_dist;

        % DO I need to save this???
        %region_data.(f_region).(f_stim).all_pulse_trig_Vm.(nr_pop) = all_pulse_vm;
        %region_data.(f_region).(f_stim).all_pulse_trace_width.(nr_pop) = all_pulse_trace_width;
        %region_data.(f_region).(f_stim).all_pulse_width_time.(nr_pop) = all_pulse_width_time;
    end

    sgtitle([f_region(3:end) ' Sub Vm Pulse-triggered average ' nr_pop], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_Pulse_Avg_Vm.png']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_Pulse_Avg_Vm.pdf']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_Pulse_Avg_Vm.eps'], 'epsc');
end

%% Spike rate pulse-trggered averaged across all pulses
% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non_entr';

fr_trig_avg_time = struct();
stats_log = [figure_path 'Small_Res' f '_Fr_pulse_triggered_times_final_average_' nr_pop ];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    all_tt = tiledlayout(length(stims), 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    all_tt = tiledlayout(length(stims), 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3*2.27, 5.16]); % Use for a sixth of a figure
    tilenum = 1;

    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        %Maybe use spike-plv_mod instead for this part of the plots
        % Check if there is an entrained field in the population data
        try
            switch nr_pop
                case 'etrain'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] > 0);
                case 'non_entr'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] < 0);
                case 'all'
                    nr_idxs = 1:length(popul_data.plv_mod_stats);
            end
        catch ME
            disp(ME.message);
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        
        num_neurons = length(nr_idxs);
        FR_avg = [];
        extra_trace = 3;
        all_pulse_fr = [];
        all_pulse_fr_by_neuron = cell(num_neurons, 1);
            
        % Save transient and sustained period pulse data
        trans_pulse_fr = [];
        trans_pulse_fr_by_neuron = cell(num_neurons, 1);
        sus_pulse_fr = [];
        sus_pulse_fr_by_neuron = cell(num_neurons, 1);

        % Looping through each neuron
        for nr = nr_idxs

            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width_time = mean(diff(popul_data.stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(popul_data.trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width_time./nr_avg_trace_time);
                
            %Neuron framerate
            nr_rate = popul_data.framerate(nr);

            % Loop through each stimulation time pulses
            for pulse_time = popul_data.stim_timestamps(:, nr)'
                follow_trace_idx = find(pulse_time <= popul_data.trace_timestamps(:, nr));
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;
                %end_trace_idx = find(pulse_time + nr_avg_pulse_width_time >= popul_data.trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                % Store the pulse width for the given neuron preserving relative stimulation period shape
                fr_pulse_width = popul_data.neuron_spikecounts_raster(start_trace_idx:end_trace_idx, nr)*nr_rate;
                all_pulse_fr_by_neuron{find(nr == nr_idxs)} = horzcat_pad(all_pulse_fr_by_neuron{find(nr == nr_idxs)}, fr_pulse_width');
                
                % Store windows depending on period of pulse, still preserving stimulation shape
                %Transient pulse
                if pulse_time >= Multi_func.trans_ped(1)/1000 & pulse_time < Multi_func.trans_ped(2)/1000
                    trans_pulse_fr_by_neuron{find(nr == nr_idxs)} = horzcat_pad(trans_pulse_fr_by_neuron{find(nr == nr_idxs)},...
                        fr_pulse_width');
                end

                %Sustained pulse
                if pulse_time >= Multi_func.sus_ped(1)/1000 & pulse_time < Multi_func.sus_ped(2)/1000
                    sus_pulse_fr_by_neuron{find(nr == nr_idxs)} = horzcat_pad(sus_pulse_fr_by_neuron{find(nr == nr_idxs)},...
                        fr_pulse_width');
                end

                % Pulse onset subtraction
                fr_pulse_width_on_sub = fr_pulse_width - popul_data.neuron_spikecounts_raster(follow_trace_idx(1), nr)*nr_rate;

                all_pulse_fr = horzcat_pad(all_pulse_fr, fr_pulse_width_on_sub);

                % Store onset SUBTRACTED windows depending on period of pulse
                %Transient pulse
                if pulse_time >= Multi_func.trans_ped(1)/1000 & pulse_time < Multi_func.trans_ped(2)/1000
                    trans_pulse_fr = horzcat_pad(trans_pulse_fr,...
                        fr_pulse_width_on_sub);
                end
                
                %Sustained pulse
                if pulse_time >= Multi_func.sus_ped(1)/1000 & pulse_time < Multi_func.sus_ped(2)/1000
                    sus_pulse_fr = horzcat_pad(sus_pulse_fr,...
                        fr_pulse_width_on_sub);
                end
            end
        end
        
        %TODO save individual neurons transient and sustained stuff
        %TODO save the two entrained and non-entrained pulse triggered to show in the single cell heatemaps

        % Save the all pulse stuff to region data
        region_data.(f_region).(f_stim).all_pulse_trig_fr = all_pulse_fr;
        region_data.(f_region).(f_stim).trans_pulse_trig_fr = trans_pulse_fr;
        region_data.(f_region).(f_stim).sus_pulse_trig_fr = sus_pulse_fr;
        

        % -- Plotting the whole period pulse average
        cur_srate = mean(all_pulse_fr, 2, 'omitnan');
        std_srate = std(all_pulse_fr, 0, 2, 'omitnan');
        
        disp([f_region ' ' f_stim]);
        num_pulses = size(all_pulse_fr, 2)
        %num_points = size(popul_data.neuron_srate_20, 1);
        sem_srate = std_srate./sqrt(num_pulses);
        
        timeline = [ [0:size(all_pulse_fr, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        nexttile(all_tt, tilenum);
        fill_h = fill([timeline, flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the whole period average    

        % -- Plotting the transient pulse average
        trans_cur_srate = mean(trans_pulse_fr, 2, 'omitnan');
        trans_std_srate = std(trans_pulse_fr, 0, 2, 'omitnan');
        
        disp([f_region ' ' f_stim]);
        num_pulses = size(trans_pulse_fr, 2)
        %num_points = size(popul_data.neuron_srate_20, 1);
        trans_sem_srate = trans_std_srate./sqrt(num_pulses);
        
        timeline = [ [0:size(trans_pulse_fr, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        nexttile(all_tt, tilenum + 1);
        fill_h = fill([timeline, flip(timeline)], [trans_cur_srate + trans_sem_srate; flipud(trans_cur_srate - trans_sem_srate)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, trans_cur_srate, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the transient pulse average    

        % -- Plotting the sustained pulse average
        sus_cur_srate = mean(sus_pulse_fr, 2, 'omitnan');
        sus_std_srate = std(sus_pulse_fr, 0, 2, 'omitnan');
        
        disp([f_region ' ' f_stim]);
        num_pulses = size(sus_pulse_fr, 2)
        %num_points = size(popul_data.neuron_srate_20, 1);
        sus_sem_srate = sus_std_srate./sqrt(num_pulses);
        
        timeline = [ [0:size(sus_pulse_fr, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        nexttile(all_tt, tilenum + 2);
        fill_h = fill([timeline, flip(timeline)], [sus_cur_srate + sus_sem_srate; flipud(sus_cur_srate - sus_sem_srate)], [0.5 0.5 0.5]);
        if ~isempty(fill_h)
            Multi_func.set_fill_properties(fill_h);
        end
        hold on;
        plot(timeline, sus_cur_srate, 'k', 'LineWidth', 1);
        hold on;
        % -- End plotting the sustained pulse average    

        
        %Shuffling data for significance testing
        wind_size = size(all_pulse_fr, 1);
        shuf_val_dist = [];
        shuf_trans_dist = [];
        shuf_sus_dist = [];
        tic;
        for i=1:1000
            iter_wind = [];
            iter_trans_wind = [];
            iter_sus_wind = [];

            for j=1:size(all_pulse_fr, 2)
                
                % Grab a neuron
                rand_nr = randi([1, length(all_pulse_fr_by_neuron)]);
                
                % Grab a window sample
                rand_win_start = randi([1, length(all_pulse_fr_by_neuron{rand_nr}) - wind_size]);
                iter_wind = horzcat_pad(iter_wind, ...
                    all_pulse_fr_by_neuron{rand_nr}(rand_win_start:rand_win_start + wind_size - 1)' ... 
                    - all_pulse_fr_by_neuron{rand_nr}(rand_win_start + extra_trace));
                
                % Same as above, but for both stimulation periods, can use the same neuron by the way
                % transient
                trans_rand_win_start = randi([1,...
                    length(trans_pulse_fr_by_neuron{rand_nr}) - wind_size]);
                iter_trans_wind = horzcat_pad(iter_trans_wind, ...
                    trans_pulse_fr_by_neuron{rand_nr}(trans_rand_win_start:trans_rand_win_start+wind_size - 1)' ...
                    - trans_pulse_fr_by_neuron{rand_nr}(trans_rand_win_start + extra_trace));
                
                % sustained
                sus_rand_win_start = randi([1,...
                    length(sus_pulse_fr_by_neuron{rand_nr}) - wind_size]);
                iter_sus_wind = horzcat_pad(iter_sus_wind, ...
                    sus_pulse_fr_by_neuron{rand_nr}(sus_rand_win_start:sus_rand_win_start+wind_size - 1)' ...
                    - sus_pulse_fr_by_neuron{rand_nr}(sus_rand_win_start + extra_trace));

            end
            end_iter = toc;
            shuf_val_dist = horzcat_pad(shuf_val_dist, mean(iter_wind, 2));
            shuf_trans_dist = horzcat_pad(shuf_trans_dist, mean(iter_trans_wind, 2));
            shuf_sus_dist = horzcat_pad(shuf_sus_dist, mean(iter_sus_wind, 2));
        end
        end_shuf = toc;

        % Save the firing rate shuffled distribution values
        region_data.(f_region).(f_stim).fr_shuf_val_dist = shuf_val_dist;

        % -- Plotting the all pulse shuffled data
        % Calculate percentile values
        low_perc = prctile(shuf_val_dist, 2.5, 2);
        high_perc = prctile(shuf_val_dist, 97.5, 2);
        shuf_mean = mean(shuf_val_dist, 2);
    
        nexttile(all_tt, tilenum);
        % Plot the shuffled values
        % Plotting using dashed lines
        plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color);
        hold on;
        plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
        hold on;

        %shade_yvals = [repmat(high_perc, 1, length(timeline)), repmat(low_perc, 1, length(timeline))];
        %f = fill([timeline, flip(timeline)], shade_yvals, [0.62 0.71 1]);
        %Multi_func.set_fill_properties(f);
        %hold on;
        %yline(shuf_mean, '--');

        % Plotting the significant time and peak
        sig_idx = find(cur_srate > high_perc);
        sig_idx(sig_idx <= extra_trace) = [];
        if ~isempty(sig_idx)
            peak_idx = find(cur_srate == max(cur_srate(sig_idx)));
            % DEBUG show the significant point
            xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
            
            diary on;
            fprintf('\n\n');
            disp([f_region ' ' f_stim]);
            disp(['Sig FR all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
            disp(['Peak FR all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
            fprintf('\n\n');
            diary off;

            region_data.(f_region).(f_stim).all_pulse_trig_fr_pk_time = timeline(peak_idx(1));
        end

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the timescale bar
        %posx = 1;
        %posy = 0;
        %plot([posx, posx + 4], [posy posy], 'k', 'LineWidth', 2);
        %text(posx, posy - 0.1, [num2str(4) 'ms']);

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        x = gca; x = x.XAxis;
        %Multi_func.set_spacing_axis(x, 6, 1);
        Multi_func.set_default_axis(gca);
        ylabel('Firing Rate(Hz)');
        xlabel('Time from pulse(ms)');
        
        %ylim([-1 10]);

        title([f_stim(3:end) ' all'], 'Interpreter', 'none');
        % -- End plotting the all pulse shuffled data

        % -- Plotting the trans pulse shuffled data
        % Calculate percentile values
        low_perc = prctile(shuf_trans_dist, 2.5, 2);
        high_perc = prctile(shuf_trans_dist, 97.5, 2);
        shuf_mean = mean(shuf_trans_dist, 2);
        
        nexttile(all_tt, tilenum + 1);
        % Plot the shuffled values as dashed horizontal lines
        plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color); 
        hold on;
        plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
        hold on;

        % Plotting the significant time and peak
        sig_idx = find(trans_cur_srate > high_perc);
        sig_idx(timeline(sig_idx) <= 0 | timeline(sig_idx) >= 1000*nr_avg_pulse_width_time) = [];

        peak_idx = find(trans_cur_srate == max(trans_cur_srate(timeline' >= 0 & timeline' <= nr_avg_pulse_width_time*1000))); % The messsy condition is for inside the pulse window and not the extra traces
        if ~isempty(sig_idx)
            % DEBUG show the significant point
            xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
            
            diary on;
            fprintf('\n\n');
            disp([f_region ' ' f_stim ' for the transient period']);
            disp(['Sig Fr all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
            disp(['Peak Fr all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
            fprintf('\n\n');
            diary off;

            % Save the peak time for the region and stimulation
            region_data.(f_region).(f_stim).trans_pulse_fr_pk_time = timeline(peak_idx(1));
        end

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Normalized Fr');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);
        
        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title([f_stim(3:end) ' transient'], 'Interpreter', 'none');
        % -- End plotting the trans pulse shuffled data

        
        % -- Plotting the sus pulse shuffled data
        % Calculate percentile values
        low_perc = prctile(shuf_sus_dist, 2.5, 2);
        high_perc = prctile(shuf_sus_dist, 97.5, 2);
        shuf_mean = mean(shuf_sus_dist, 2);
        
        nexttile(all_tt, tilenum + 2);
        % Plot the shuffled values as dashed horizontal lines
        plot(timeline, [low_perc, high_perc], '--', 'Color', Multi_func.shuf_color); 
        hold on;
        plot(timeline, shuf_mean, 'Color', Multi_func.shuf_color);
        hold on;

        % Plotting the significant time and peak
        sig_idx = find(sus_cur_srate > high_perc);
        sig_idx(timeline(sig_idx) <= 0 | timeline(sig_idx) >= 1000*nr_avg_pulse_width_time) = [];

        peak_idx = find(sus_cur_srate == max(sus_cur_srate(timeline' >= 0 & timeline' <= nr_avg_pulse_width_time*1000))); % The messsy condition is for inside the pulse window and not the extra traces
        peak_idx = peak_idx(find(timeline(peak_idx) >= 0 & timeline(peak_idx) <= nr_avg_pulse_width_time*1000)); % Keep only the idxs that are within in the time window
        if ~isempty(sig_idx)
            % DEBUG show the significant point
            xline([timeline(peak_idx(1)), timeline(sig_idx(1))], 'g');
            
            diary on;
            fprintf('\n\n');
            disp([f_region ' ' f_stim ' for the sustained period']);
            disp(['Sig Fr all pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
            disp(['Peak Fr all pulse at ' num2str(timeline(peak_idx(1))) 'ms']);
            fprintf('\n\n');
            diary off;

            % Save the peak time for the region and stimulation
            region_data.(f_region).(f_stim).sus_pulse_fr_pk_time = timeline(peak_idx(1));
        end

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
        hold on;

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Normalized Fr');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);
        
        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title([f_stim(3:end) ' sustained'], 'Interpreter', 'none');
        % -- End plotting the sus pulse shuffled data

        % Increase the tile number indicator
        tilenum = tilenum + 3;
    end
    sgtitle([f_region(3:end) ' Firing Rate pulse-triggered average ' nr_pop], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_Pulse_Avg_FR.png']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_Pulse_Avg_FR.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.eps'], 'epsc');
end

%% All region All Pulse Vm overlay
test_regions = fieldnames(region_data)';
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    nexttile;

    %Loop through regions
    for f_region = test_regions
        f_region = f_region{1};
        stim_data = region_data.(f_region).(f_stim);

        timeline = [ [0:size(stim_data.all_pulse_trig_Vm, 1) - 1] - extra_trace]*1000./avg_Fs;
        cur_trig_vm = nanmean(stim_data.all_pulse_trig_Vm, 2);
        std_trig_vm = std(stim_data.all_pulse_trig_Vm, 0, 2, 'omitnan');
        num_pulses = size(stim_data.all_pulse_trig_Vm, 2);
        %num_points = size(data_bystim.(f_stim).neuron_trig_vm, 1);
        sem_trig_vm = std_trig_vm./sqrt(num_pulses);

        fill_h = fill([timeline, flip(timeline)], [cur_trig_vm + sem_trig_vm; flipud(cur_trig_vm - sem_trig_vm)], [0.5 0.5 0.5], 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        hold on;
    
        if strcmp(f_region, 'r_CA1') == 1
            cur_color = Multi_func.ca1_color;
            label = 'CA1';
        elseif strcmp(f_region, 'r_M1') == 1
            cur_color = Multi_func.m1_color;
            label = 'M1';
        elseif strcmp(f_region, 'r_V1') == 1
            cur_color = Multi_func.v1_color;
            label = 'V1';
        end

        plot(timeline, cur_trig_vm, 'Color', cur_color, 'LineWidth', 1, 'DisplayName', label);
        hold on;
        
        % Plot the shuffled distribution
        shuf_val_dist = stim_data.vm_shuf_val_dist; %TODO this shuffled depends on which modulated population of neurons is being used
        low_perc = prctile(shuf_val_dist, 2.5, 2);
        high_perc = prctile(shuf_val_dist, 97.5, 2);
        shuf_mean = mean(shuf_val_dist, 2);

        % Plot the shuffled values as dashed horizontal lines
        plot(timeline, [low_perc, high_perc], '--', 'Color', cur_color, 'HandleVisibility', 'off');
        hold on;
        plot(timeline, shuf_mean, 'Color', cur_color, 'HandleVisibility', 'off');
        hold on;
    end

    % Plot the pulses
    pulse_times = [0, mean(stim_data.all_pulse_width_time, 2)*1000]
    xline(pulse_times, 'Color', Multi_func.pulse_color, 'LineWidth', 2, 'HandleVisibility', 'off');
    hold on;

    legend('Location', 'northwest');
    Multi_func.set_default_axis(gca);
    title(f_stim, 'Interpreter', 'none');
    xlabel('time (ms)');
    ylabel('Normalized Vm Change');
end
sgtitle('All Pulse Vm Plot by Region Overlay');
saveas(gcf, [figure_path 'Small_Res' f 'Pulse_Triggered_Region_Overlay_Vm.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Pulse_Triggered_Region_Overlay_Vm.pdf']);

%% Compare the heights of each Vm fluctuation based on the time of peak M1
stats_log = [figure_path 'Small_Res' f 'Pulse_triggered_Vm_heights_from_peak_time'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 20, 20]);
stims = fieldnames(region_data.r_M1)';
for f_stim = stims
    f_stim = f_stim{1};
    m1_pulses_vm = region_data.r_M1.(f_stim).all_pulse_trig_Vm;
    ca1_pulses_vm = region_data.r_CA1.(f_stim).all_pulse_trig_Vm;
    timeline = [[0:size(m1_pulses_vm, 1) - 1] - extra_trace]*1000./avg_Fs;

    pk_t_comp = region_data.r_M1.(f_stim).all_pulse_vm_pk_time;
    pk_idx = find(timeline == pk_t_comp);
    
    % Plot the violin of all the M1 height at time of peak
    m1_vm_heights = m1_pulses_vm(pk_idx, :);
    ca1_vm_heights = ca1_pulses_vm(pk_idx, :);
    nexttile;

    % Setup variables for violinplots
    data = [m1_vm_heights, ca1_vm_heights];
    labels = [repmat({'M1'}, 1, length(m1_vm_heights)), repmat({'CA1'}, 1, length(ca1_vm_heights))];
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data, labels, 'GroupOrder', {'M1', 'CA1'}, ViolinOpts);

    %violins(1).ViolinColor = {'k'};
    %violins(2).ViolinColor = {'g'};
    violins(1).ViolinColor = {Multi_func.m1_color};
    violins(2).ViolinColor = {Multi_func.ca1_color};

    % Perform the statisical for significant difference
    diary on;
    disp(['Pulse Vm heights at peak for ' num2str(f_stim)]);
    [p, h, stats] = ranksum(m1_vm_heights, ca1_vm_heights);
    disp(['CA1 heights ' num2str(nanmean(ca1_vm_heights))]);
    disp(['M1 heights ' num2str(nanmean(m1_vm_heights))]);
    fprintf('\n\n');
    diary off;
    legend(['p= ' num2str(p)]);    

    title([f_stim(3:end) ' M1 and CA1 Vm hieghts comparison'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Vm_Heights_ca1_vs_m1.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Vm_Heights_ca1_vs_m1.pdf']);

%% Compare the heights of each firing rate fluctuation based on the time of peak M1
stats_log = [figure_path 'Small_Res' f 'Pulse_triggered_Fr_heights_from_peak_time'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 20, 20]);
stims = fieldnames(region_data.r_M1)';
for f_stim = stims
    f_stim = f_stim{1};
    m1_pulses_fr = region_data.r_M1.(f_stim).all_pulse_trig_fr;
    ca1_pulses_fr = region_data.r_CA1.(f_stim).all_pulse_trig_fr;
    timeline = [[0:size(m1_pulses_fr, 1) - 1] - extra_trace]*1000./avg_Fs;

    pk_t_comp = region_data.r_M1.(f_stim).all_pulse_trig_fr_pk_time;
    pk_idx = find(timeline == pk_t_comp);
    
    % Plot the violin of all the M1 height at time of peak
    m1_fr_heights = m1_pulses_fr(pk_idx, :);
    ca1_fr_heights = ca1_pulses_fr(pk_idx, :);
    nexttile;

    % Setup variables for violinplots
    data = [m1_fr_heights, ca1_fr_heights];
    labels = [repmat({'M1'}, 1, length(m1_fr_heights)), repmat({'CA1'}, 1, length(ca1_fr_heights))];
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data, labels, 'GroupOrder', {'M1', 'CA1'}, ViolinOpts);

    %violins(1).ViolinColor = {'k'};
    %violins(2).ViolinColor = {'g'};
    violins(1).ViolinColor = {Multi_func.m1_color};
    violins(2).ViolinColor = {Multi_func.ca1_color};

    % Perform the statisical for significant difference
    diary on;
    disp(['Pulse Fr heights at peak for ' num2str(f_stim)]);
    [p, h, stats] = ranksum(m1_fr_heights, ca1_fr_heights);
    disp(['CA1 heights ' num2str(nanmean(ca1_fr_heights))]);
    disp(['M1 heights ' num2str(nanmean(m1_fr_heights))]);
    fprintf('\n\n');
    diary off;
    legend(['p= ' num2str(p)]);    

    title([f_stim(3:end) ' M1 and CA1 Fr hieghts comparison'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Fr_Heights_ca1_vs_m1.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Fr_Heights_ca1_vs_m1.pdf']);

%% Time to peak Vm during pulse-triggered average 
stats_log = [figure_path 'Small_Res' f 'Vm_neuronwise_pulse_triggered_region_times'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;

stats_data = struct();

figure('Position', [140 150 , 1000, 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 3.5, 5]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
for f_stim = stims'
    f_stim = f_stim{1};
    stats_data.(f_stim) = struct();
    stats_data.(f_stim).reg_order = {};
    stats_data.(f_stim).data = [];

    diary on;
    disp([f_stim(3:end) ' Hz Vm']);
    diary off;

    nexttile;
    % Loop through each brain region
    for f_region = fieldnames(region_data)'
        f_region = f_region{1};
        % Store the time to peak values for both regions
        reg_pk_time = [];
        
        reg_pulses = region_data.(f_region).(f_stim).all_pulse_trig_Vm(extra_trace + 1:end - extra_trace - 1, :);
        timeline = [0:size(reg_pulses) - 1]*1000./avg_Fs;
        for i=1:size(reg_pulses, 2)
            peak_idx = find(reg_pulses(:, i) == max(reg_pulses(:, i)));
            if timeline(peak_idx) == 0
                continue;
            end
            reg_pk_time(end + 1) = timeline(peak_idx);
        end

        diary on;
        num_pk_times = length(reg_pk_time);
        disp([f_region(3:end) ' ' num2str(nanmean(reg_pk_time)) '+-' ... 
                 num2str(nanstd(reg_pk_time)/sqrt(num_pk_times)) ...
                 ' n=' num2str(num_pk_times)]);
        diary off;
        histogram(reg_pk_time, 100, 'DisplayName', f_region(3:end)); 
        hold on;

        % Save the peak times
        stats_data.(f_stim).reg_order = cat(2, stats_data.(f_stim).reg_order, repmat({f_region}, 1, length(reg_pk_time)));
        stats_data.(f_stim).data = [stats_data.(f_stim).data, reg_pk_time];
    end
    diary on;
    fprintf('\n\n');
    diary off;
    legend();
    title([num2str(f_stim) ' Vm time to peak histogram'], 'Interpreter', 'none');
    diary on;
    [p, h, stats] = kruskalwallis(stats_data.(f_stim).data, stats_data.(f_stim).reg_order)
    disp('Group Columns');
    disp(stats.gnames');
    disp('The first two columns indicate the pairs to compare between regions as indicated above.');
    disp('Last column is the p-value');
    c = multcompare(stats, 'CriticalValueType', 'dunn-sidak')
    diary off;
end
saveas(gcf, [figure_path 'Small_Res' f 'Histogram_Pulse_Triggered_Peak_Time_Vm_Regions.pdf']);
saveas(gcf, [figure_path 'Small_Res' f 'Histogram_Pulse_Triggered_Peak_Time_Vm_Regions.png']);

%% All region all pulse firing rate overlay
test_regions = fieldnames(region_data)';
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    nexttile;

    %Loop through regions
    for f_region = test_regions
        f_region = f_region{1};
        stim_data = region_data.(f_region).(f_stim);

        timeline = [ [0:size(stim_data.all_pulse_trig_fr, 1) - 1] - extra_trace]*1000./avg_Fs;
        cur_trig_fr = nanmean(stim_data.all_pulse_trig_fr, 2);
        std_trig_fr = std(stim_data.all_pulse_trig_fr, 0, 2, 'omitnan');
        num_pulses = size(stim_data.all_pulse_trig_fr, 2);
        %num_points = size(data_bystim.(f_stim).neuron_trig_vm, 1);
        sem_trig_fr = std_trig_fr./sqrt(num_pulses);

        fill_h = fill([timeline, flip(timeline)], [cur_trig_fr + sem_trig_fr; flipud(cur_trig_fr - sem_trig_fr)], [0.5 0.5 0.5], 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        hold on;

        if strcmp(f_region, 'r_CA1') == 1
            cur_color = Multi_func.ca1_color;
            label = 'CA1';
        elseif strcmp(f_region, 'r_M1') == 1
            cur_color = Multi_func.m1_color;
            label = 'M1';
        elseif strcmp(f_region, 'r_V1') == 1
            cur_color = Multi_func.v1_color;
            label = 'V1';
        end
        plot(timeline, cur_trig_fr, 'Color', cur_color, 'LineWidth', 1, 'DisplayName', label);
        hold on;
            
        % Plot the shuffled distribution
        shuf_val_dist = stim_data.fr_shuf_val_dist;
        low_perc = prctile(shuf_val_dist, 2.5, 2);
        high_perc = prctile(shuf_val_dist, 97.5, 2);
        shuf_mean = mean(shuf_val_dist, 2);

        % Plot the shuffled values as dashed horizontal lines
        plot(timeline, [low_perc, high_perc], '--', 'Color', cur_color, 'HandleVisibility', 'off');
        hold on;
        plot(timeline, shuf_mean, 'Color', cur_color, 'HandleVisibility', 'off');
        hold on;
    end
 
    % Plot the time pulses
    pulse_times = [0, mean(stim_data.all_pulse_width_time, 2)*1000]
    % Display the pulses
    xline(pulse_times, 'Color', Multi_func.pulse_color, 'LineWidth', 2, 'HandleVisibility', 'off');
    hold on;
    
    legend('Location', 'northeast');
    Multi_func.set_default_axis(gca);
    title(f_stim, 'Interpreter', 'none');
    xlabel('time (ms)');
    ylabel('Firing Rate Change (Hz)');
end
sgtitle('All Pulse Fr Plot by Region Overlay');
saveas(gcf, [figure_path 'Small_Res' f 'Pulse_Triggered_Region_Overlay_Fr.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Pulse_Triggered_Region_Overlay_Fr.pdf']);

%% -- Time to peak pulse triggered Firing Rate between regions and compare statistical test
stats_log = [figure_path 'Small_Res' f 'FR_neuronwise_pulse_triggered_region_times'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

% Store all of the statistical data
stats_data = struct();

figure('Position', [140 150 , 1000, 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
for f_stim = stims'
    f_stim = f_stim{1};
    stats_data.(f_stim) = struct();
    stats_data.(f_stim).reg_order = {};
    stats_data.(f_stim).data = [];

    diary on;
    disp([f_stim(3:end) ' Hz FR']);
    diary off;
    nexttile;

    for f_region = fieldnames(region_data)'
        f_region = f_region{1};
        
        % Store the time to peak values for both regions
        reg_pk_time = [];
    
        reg_pulses = region_data.(f_region).(f_stim).all_pulse_trig_fr(extra_trace + 1:end - extra_trace - 1, :);
        timeline = [0:size(reg_pulses) - 1]*1000./avg_Fs;
        for i=1:size(reg_pulses, 2)
            peak_idx = find(reg_pulses(:, i) == max(reg_pulses(:, i)));
            if length(peak_idx) > 1
                peak_idx = peak_idx(1);
            end
            if timeline(peak_idx) == 0
                continue;
            end
            reg_pk_time(end + 1) = timeline(peak_idx);
        end

        % Print data to file and commandline
        diary on;
        num_pk_times = length(reg_pk_time);
        disp([f_region(3:end) ' ' num2str(nanmean(reg_pk_time)) '+-' ... 
                 num2str(nanstd(reg_pk_time)/sqrt(num_pk_times)) ...
                 ' n=' num2str(num_pk_times)]);
        diary off;
        histogram(reg_pk_time, 100, 'DisplayName', f_region(3:end)); 
        hold on;
    
        % Save the peak times to a structure
        stats_data.(f_stim).reg_order = cat(2, stats_data.(f_stim).reg_order, repmat({f_region}, 1, length(reg_pk_time)));
        stats_data.(f_stim).data = [stats_data.(f_stim).data, reg_pk_time];
    end

    diary on;
    fprintf('\n\n');
    diary off;
    legend();
    title([num2str(f_stim) ' FR time to peak histogram'], 'Interpreter', 'none');
    diary on;
    [p, h, stats] = kruskalwallis(stats_data.(f_stim).data, stats_data.(f_stim).reg_order)
    disp('Group Columns');
    disp(stats.gnames');
    c = multcompare(stats, 'CriticalValueType', 'dunn-sidak')
    diary off;
end
saveas(gcf, [figure_path 'Small_Res' f 'Histogram_Pulse_Triggered_Peak_Time_Firing_Rate_Region.pdf']);
saveas(gcf, [figure_path 'Small_Res' f 'Histogram_Pulse_Triggered_Peak_Time_Firing_Rate_Region.png']);

%% Firing rate regression CA1 vs. M1
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(2, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
tiledlayout(length(stims), 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 20, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    num_stim = str2num(f_stim(3:end));
    nexttile;
    
    % Grab all pulse firing rate stuff 
    ca1_pulses_fr = region_data.r_CA1.(f_stim).all_pulse_trig_fr((extra_trace+1):(end - extra_trace), :);
    m1_pulses_fr = region_data.r_M1.(f_stim).all_pulse_trig_fr((extra_trace+1):(end - extra_trace), :);

    % First test averages
    ca1_fr_avg = nanmean(ca1_pulses_fr, 2);
    m1_fr_avg = nanmean(m1_pulses_fr, 2);

    % Plot the scatter of Ca1 FR to M1 FR
    plot(ca1_fr_avg, m1_fr_avg, '.');
    hold on;
    
    fitresults = polyfit(ca1_fr_avg, m1_fr_avg, 1);
    x_vals = ca1_fr_avg;
    y_orig = m1_fr_avg;
    fit_y = polyval(fitresults, x_vals);
    plot(x_vals, fit_y, '-');
    pearson_coeff = corrcoef(y_orig, fit_y);
    legend(['r= ' num2str(pearson_coeff(1, 2))]);

    xlabel('CA1 Firing rate');
    ylabel('M1 Firing rate');
    title([f_stim ' using average pulse'], 'Interpreter', 'none');

    nexttile;
    
    % Use the pulse trace width
    p_trace_width = size(ca1_pulses_fr, 1);
    
    ca1_pulses_fr = reshape(ca1_pulses_fr, p_trace_width*num_stim, []);
    m1_pulses_fr = reshape(m1_pulses_fr, p_trace_width*num_stim, []);
    ca1_pulse_nr  = reshape(nanmean(ca1_pulses_fr, 2), p_trace_width, []);
    m1_pulse_nr  = reshape(nanmean(m1_pulses_fr, 2), p_trace_width, []);

    ca1_size = size(ca1_pulse_nr)
    m1_size = size(m1_pulse_nr)

    % Plot the scatter of Ca1 FR to M1 FR
    plot(ca1_pulse_nr, m1_pulse_nr, '.');
    hold on;

    linearize_ca1 = ca1_pulse_nr(:);
    linearize_m1 = m1_pulse_nr(:);

    fitresults = polyfit(linearize_ca1, linearize_m1, 1);
    x_vals = linearize_ca1;
    y_orig = linearize_m1;
    fit_y = polyval(fitresults, x_vals);
    plot(x_vals, fit_y, '-');
    pearson_coeff = corrcoef(y_orig, fit_y)
    legend(['r= ' num2str(pearson_coeff(1, 2))])

    xlabel('CA1 Firing rate');
    ylabel('M1 Firing rate');
    title([f_stim ' regressed over individual pulses and neuronwise average'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Regress_Region_Pulse_Triggered_Fr.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Regress_Region_Pulse_Triggered_Fr.pdf']);

%% Vm regression CA1 vs. M1
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(2, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
tiledlayout(length(stims), 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 20, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    num_stim = str2num(f_stim(3:end));
    nexttile;
    
    % Grab all pulse firing rate stuff 
    ca1_pulses_vm = region_data.r_CA1.(f_stim).all_pulse_trig_Vm((extra_trace+1):(end - extra_trace), :);
    m1_pulses_vm = region_data.r_M1.(f_stim).all_pulse_trig_Vm((extra_trace+1):(end - extra_trace), :);

    % First test averages
    ca1_vm_avg = nanmean(ca1_pulses_vm, 2);
    m1_vm_avg = nanmean(m1_pulses_vm, 2);

    % Plot the scatter of Ca1 FR to M1 FR
    plot(ca1_vm_avg, m1_vm_avg, '.');
    hold on;
    
    fitresults = polyfit(ca1_vm_avg, m1_vm_avg, 1);
    x_vals = ca1_vm_avg;
    y_orig = m1_vm_avg;
    fit_y = polyval(fitresults, x_vals);
    plot(x_vals, fit_y, '-');
    pearson_coeff = corrcoef(y_orig, fit_y);
    legend(['r= ' num2str(pearson_coeff(1, 2))]);

    xlabel('CA1 Vm');
    ylabel('M1 Vm');
    title([f_stim ' using average pulse'], 'Interpreter', 'none');

    nexttile;
    
    % Use the pulse trace width
    p_trace_width = size(ca1_pulses_vm, 1);
    
    ca1_pulses_vm = reshape(ca1_pulses_vm, p_trace_width*num_stim, []);
    m1_pulses_vm = reshape(m1_pulses_vm, p_trace_width*num_stim, []);
    ca1_pulse_nr  = reshape(nanmean(ca1_pulses_vm, 2), p_trace_width, []);
    m1_pulse_nr  = reshape(nanmean(m1_pulses_vm, 2), p_trace_width, []);

    ca1_size = size(ca1_pulse_nr)
    m1_size = size(m1_pulse_nr)

    % Plot the scatter of Ca1 Vm to M1 Vm
    plot(ca1_pulse_nr, m1_pulse_nr, '.');
    hold on;

    linearize_ca1 = ca1_pulse_nr(:);
    linearize_m1 = m1_pulse_nr(:);

    fitresults = polyfit(linearize_ca1, linearize_m1, 1);
    x_vals = linearize_ca1;
    y_orig = linearize_m1;
    fit_y = polyval(fitresults, x_vals);
    plot(x_vals, fit_y, '-');
    pearson_coeff = corrcoef(y_orig, fit_y)
    legend(['r= ' num2str(pearson_coeff(1, 2))])

    xlabel('CA1 Vm');
    ylabel('M1 Vm');
    title([f_stim ' regressed over individual pulses and neuronwise average'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Regress_Region_Pulse_Triggered_Vm.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Regress_Region_Pulse_Triggered_Vm.pdf']);

%%---- Cross-correlation of the final all pulse Vm plot!! ----
figure('Position', [0 0 , 1000, 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 3.5, 5]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    nexttile;
    % Chopping off the extra points flanking the window
    ca1_all_pulse_vm = region_data.r_CA1.(f_stim).all_pulse_trig_Vm(extra_trace + 1:end - extra_trace, :);
    ca1_pulse_vm = nanmean(ca1_all_pulse_vm, 2);
        
    m1_all_pulse_vm = region_data.r_M1.(f_stim).all_pulse_trig_Vm(extra_trace + 1:end - extra_trace, :);
    m1_pulse_vm = nanmean(m1_all_pulse_vm, 2);

    [c, lags] = xcorr(ca1_pulse_vm, m1_pulse_vm);
    lags = (lags/avg_Fs)*1000; % Convert lags to ms
    plot(lags, c);
    hold on;
    [~, peak_idx] = findpeaks(c);
    labels = cellstr(strsplit(num2str(lags(peak_idx)), ' '));
    xline(lags(peak_idx), '--',  labels,'LabelOrientation', 'horizontal');
    xlabel('time lag (ms)');
    ylabel('Correlation');

    title(f_stim, 'Interpreter', 'none');
end
sgtitle('Time laggs cross-correlation');

%-- Perform lag cross-correlation across all pairs
%figure('Position', [300 300, 1000, 1000]);
%%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 3.5, 5]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
%% First loop through stimulations
%for f_stim=stims
%    f_stim = f_stim{1};
%    nexttile;
%    % Chopping off the extra points flanking the window
%    ca1_all_pulse_vm = region_data.r_CA1.(f_stim).all_pulse_trig_Vm(extra_trace + 1:end - extra_trace, :);
%    ca1_pulse_vm = nanmean(ca1_all_pulse_vm, 2);
%        
%    m1_all_pulse_vm = region_data.r_M1.(f_stim).all_pulse_trig_Vm(extra_trace + 1:end - extra_trace, :);
%    m1_pulse_vm = nanmean(m1_all_pulse_vm, 2);
%
%    all_cc = [];
%    all_lags = [];
%
%    %Suggestion, I should create the repeated arrays so that I can just do arrayfun(I) with xcorr at the end
%    % This maxes out array size
%    tic
%    %Start CA1 pulses loop
%    for ca1_i=1:size(ca1_all_pulse_vm, 2)
%        %M1 pulses loop
%        for m1_i=1:size(m1_all_pulse_vm, 2)
%            [c, lags] = xcorr(ca1_all_pulse_vm(:, ca1_i), m1_all_pulse_vm(:, m1_i));
%            all_cc = horzcat_pad(all_cc, c(:));
%            all_lags = horzcat_pad(all_lags, lags(:));
%        end
%    end
%    toc
%    
%    [c, lags] = xcorr(ca1_pulse_vm, m1_pulse_vm);
%    lags = (lags/avg_Fs)*1000; % Convert lags to ms
%    plot(lags, c);
%    hold on;
%    [~, peak_idx] = findpeaks(c);
%    labels = cellstr(strsplit(num2str(lags(peak_idx)), ' '));
%    xline(lags(peak_idx), '--',  labels,'LabelOrientation', 'horizontal');
%    xlabel('time lag (ms)');
%    ylabel('Correlation');
%
%    title(f_stim, 'Interpreter', 'none');
%end
%sgtitle('Time laggs cross-correlation all pairs');
% -- End of cross-correlation for all pairs

%% Clustering method for comparing significance between the regions
%Suggestion, maybe do my own kind of clustering and then compare those results to permuted cluster

% Or just do summed t between observed and permutated distribution

% Using the t-statistic did not seem to work
% Clustered permutation test
%num_shuf = 50; %Suggestion, change this to like a thousand
%data1 = region_data.r_CA1.f_40.all_pulse_trig_Vm(extra_trace:end - extra_trace, :);
%data2 = region_data.r_M1.f_40.all_pulse_trig_Vm(extra_trace:end - extra_trace, :);
%alpha = 0.5;
%
%t_vals = [];
%% Iterate over time points
%for time_idx = 1:size(data1, 1)
%
%    % Perform two-sample t-test
%    [~, ~, ~, stats] = ttest2(data1(time_idx, :), data1(time_idx, :))
%    
%    % Store the t-values for each time point
%    t_vals(time_idx) = stats.tstat;
%end
%clusters = bwlabel(t_vals > 0.5);
%
%%DEBUG
%figure();
%plot(clusters)
%
%cluster_stats = accumarray(clusters, t_vals, [], @sum);


%% Spike rate showing display all pulses
%for f_region = fieldnames(region_data)'
%    f_region = f_region{1};
%    data_bystim = region_data.(f_region);
%    stims = fieldnames(data_bystim);
%    
%    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%    for f_stim=stims'
%        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
%        cur_srate = mean(data_bystim.(f_stim).neuron_srate_20, 2, 'omitnan');
%        std_srate = std(data_bystim.(f_stim).neuron_srate_20, 0, 2, 'omitnan');
%        num_neurons = size(data_bystim.(f_stim).neuron_srate_20, 2);
%        %num_points = size(data_bystim.(f_stim).neuron_srate, 1);
%        sem_srate = std_srate./sqrt(num_neurons);
%        nexttile;
%        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
%        Multi_func.set_fill_properties(f);
%        hold on;
%        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
%        hold on;
%
%        % Plot the DBS stimulation time pulses
%        stim_time = nanmean(data_bystim.(f_stim).stim_timestamps, 2);
%        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
%        hold on;
%
%        % Plot the timescale bar
%        posx = -.100;
%        posy = 0;
%        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
%        text(posx, posy - 0.2, '50ms');
%
%        % Increase timescale resolution
%        xlim([0 - .10, max(stim_time) + .10]);
%        a = gca;
%        a.XAxis.Visible = 'off';
%        set(gca, 'color', 'none');
%        ylabel('Firing Rate(Hz)');
%        title(f_stim(3:end), 'Interpreter', 'none');
%    end
%    sgtitle([f_region '_' num2str(srate_win)  ' Population Spike rate with all pulse'], 'Interpreter', 'none');
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.png']);
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.pdf']);
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.eps'], 'epsc');
%end

%% Subthreshold Vm first few pulses and middle sustained period
% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

vm_onset_time = struct();
stats_log = [figure_path 'Small_Res' f 'Vm_onset_sig_times_' nr_pop];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    %tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.40, 4.96]); % Use for a 4th of a figure length
    %tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 2.27, 4.96]); % Use for a 6th of the figure
    tiledlayout(length(stims), 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 9, 4.96]); % Use for matching the same length as all DBS pulses
    
    for f_stim=stims'
        f_stim = f_stim{1};
        
        popul_data = data_bystim.(f_stim);
        
        % Check if there is an entrained field in the population data
        try
            switch nr_pop
                case 'etrain'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] > 0);
                case 'non'
                    nr_idxs = find([popul_data.plv_mod_stats.mod] < 0);
                case 'all'
                    nr_idxs = 1:length(popul_data.plv_mod_stats);
            end
        catch ME
            disp(ME.message);
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        norm_vms = popul_data.neuron_RawVm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(popul_data.neuron_RawVm, 1);
        
        % Grab the min and max value for later plotting
        min_yval = min(cur_Vm - sem_Vm);
        max_yval = max(cur_Vm + sem_Vm);

        nexttile;

        % Standard Error
        fill_h = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        
        % Plot Vm
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        hold on;
        
        % Plot the DBS stimulation time pulses
        xline(nanmean(popul_data.stim_timestamps, 2), 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the semi-significant points on the subthreshold Vm

        base_Vm = [];
        % Grab the baseline sub Vm for all neurons
        for i = nr_idxs
            baseline_idx = find(popul_data.trace_timestamps(:, find(i == nr_idxs)) < popul_data.stim_timestamps(1, find(i == nr_idxs)));
            base_Vm = horzcat_pad(base_Vm, norm_vms(baseline_idx, find(i == nr_idxs) ) );
        end
        
        % Calculate each Vm point's significance based on the population average of all neurons
        base_Vm = mean(base_Vm, 2, 'omitnan');
        std_baseline = std(base_Vm, 0, 'omitnan'); 
        sig_idx = find(cur_Vm > (2*std_baseline + mean(base_Vm, 'omitnan')));
            
        sig_idx((timeline(sig_idx) < 0) | (timeline(sig_idx) > 100)) = [];

        % Calculate a good height for marking significance
        thres = 0.2;
        height = range(cur_Vm);
        height = (1+thres)*height;
        %plot(timeline(sig_idx), repmat(height, 1, length(sig_idx)), '.b', 'MarkerSize', 8);
        
        % Plotting the 2 std line 
        yline(2*std_baseline, '--', 'Color', [0 0 0 0.5]);

        % Plot and calculate the time to first significance
        if ~isempty(sig_idx)
            %DEBUG show where the first significant point is
            xline(timeline(sig_idx(1)), 'g');

            diary on;
            fprintf('\n\n');
            disp([f_region ' ' f_stim]);
            disp(['Sig Vm first pulse at ' num2str(timeline(sig_idx(1))) 'ms']);
            fprintf('\n\n');
            diary off;
        end

        % Increase timescale resolution
        %TODO figure out a better way to get the Y axis limits
        xlim([Multi_func.onset_ped(1), Multi_func.onset_ped(2)]./1000);
        ylim([min_yval, max_yval]);
        ylabel('Normalized Vm');
        xlabel('Time from onset(ms)');

        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' onset wind'], 'Interpreter', 'none');
        
        % Copy the current plot and switch to another window snippet
        ax = gca;
        nexttile;
        copyobj(allchild(ax), gca);
        ylim(get(ax, 'YLim'));
        Multi_func.set_default_axis(gca);

        ylabel('Normalized Vm');
        xlabel('Time from onset(ms)');
        
        xlim([Multi_func.mid_stim_ped(1), Multi_func.mid_stim_ped(2)]./1000);
        title([f_stim(3:end) ' sus wind'], 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Average subthreshold Vm ' nr_pop], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_onset_and_sus_wind_pulse_Vm.png']);
    saveas(gcf, [figure_path 'Small_Res' f f_region '_' nr_pop '_onset_and_sus_wind_pulse_Vm.pdf']);
end

%% All regions overlay first few pulses Vm 
test_regions = {'r_CA1', 'r_M1', 'r_V1'};
stims = fieldnames(region_data.r_M1)';
figure('Position', [0 0 , 1000, 1000]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 8.3, 3.5]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);
% First loop through stimulations
for f_stim=stims
    f_stim = f_stim{1};
    nexttile;

    %Loop through regions
    for f_region = test_regions
        f_region = f_region{1};
        stim_data = region_data.(f_region).(f_stim);

        timeline = nanmean(stim_data.trace_timestamps, 2)*1000; % Convert to ms
        norm_vms = stim_data.neuron_Vm./stim_data.neuron_spike_amp;
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(stim_data.neuron_Vm, 1);
 
        % Standard Error
        fill_h = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5], 'HandleVisibility','off');
        Multi_func.set_fill_properties(fill_h);
        hold on;
        
        if strcmp(f_region, 'r_CA1') == 1
            cur_color = Multi_func.ca1_color;
            label = 'CA1';
        elseif strcmp(f_region, 'r_M1') == 1
            cur_color = Multi_func.m1_color;
            label = 'M1';
        elseif strcmp(f_region, 'r_V1') == 1
            cur_color = Multi_func.v1_color;
            label = 'V1';
        end

        % Plot Vm
        plot(timeline, cur_Vm, 'Color', cur_color, 'LineWidth', 1, 'DisplayName', label);
        hold on;
    end
        
    % Plot the DBS stimulation time pulses
    xline(nanmean(stim_data.stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5, 'HandleVisibility','off');
    hold on;

    % Increase timescale resolution
    xlim([0 - 50, 0 + 100]);
    ylabel('Normalized Vm');
    xlabel('Time from onset(ms)');
    
    legend('Location', 'southeast');
    Multi_func.set_default_axis(gca);
    title([f_stim(3:end) ' Hz'], 'Interpreter', 'none');
end

sgtitle(['Onset Vm Subthreshold Region Compare'], 'Interpreter', 'none');
saveas(gcf, [figure_path 'Small_Res' f 'Onset_Overlay_Vm.png']);
saveas(gcf, [figure_path 'Small_Res' f 'Onset_Overlay_Vm.pdf']);

%% Time to peak Vm onset between regions and statistical test
stats_log = [figure_path 'Small_Res' f 'Vm_pop_onset_ca1_vs_m1_times'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;
test_regions = {'r_CA1', 'r_M1'};
stims = fieldnames(region_data.r_M1)';
figure('Position', [140 150 , 1000, 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
for f_stim = stims
    f_stim = f_stim{1};
    
    % Store the time to peak values for both regions
    ca1_pk_time = [];
    m1_pk_time = [];
    ca1_vm = [];
    m1_vm = [];
    
    % Grab just the first 100 ms onset trace times
    timeline = nanmean(stim_data.trace_timestamps, 2)*1000;
    onset_idx = find(timeline >=0 & timeline <= 100); % Only grab the first 100 ms
    timeline = timeline(onset_idx);
    ca1_vm = region_data.r_CA1.(f_stim).neuron_Vm(onset_idx, :)./region_data.r_CA1.(f_stim).neuron_spike_amp;
    m1_vm = region_data.r_M1.(f_stim).neuron_Vm(onset_idx, :)./region_data.r_M1.(f_stim).neuron_spike_amp;
    
    for i=1:size(ca1_vm, 2)
        peak_idx = find(ca1_vm(:, i) == max(ca1_vm(:, i)));
        if length(peak_idx) > 1
            peak_idx = peak_idx(1);
        end
        if timeline(peak_idx) == 0
            continue;
        end
        ca1_pk_time(end + 1) = timeline(peak_idx);
    end

    for i=1:size(m1_vm, 2)
        peak_idx = find(m1_vm(:, i) == max(m1_vm(:, i)));
        if length(peak_idx) > 1
            peak_idx = peak_idx(1);
        end
        if timeline(peak_idx) == 0
            continue;
        end
        m1_pk_time(end + 1) = timeline(peak_idx);
    end

    diary on;
    disp(['Time to peak ' f_stim]);
    [p, h, stats] = ranksum(m1_pk_time, ca1_pk_time)
    num_ca1_pk_times = length(ca1_pk_time);
    num_m1_pk_times = length(m1_pk_time);
    disp(['CA1 num traces ' num2str(num_ca1_pk_times)]);
    disp(['CA1 mean:' num2str(nanmean(ca1_pk_time))]);
    disp(['CA1 sem: ' num2str(nanstd(ca1_pk_time)/sqrt(num_ca1_pk_times))]);
    disp(['CA1 std:' num2str(nanstd(ca1_pk_time))]);
    disp(['M1 num traces ' num2str(num_m1_pk_times)]);
    disp(['M1 mean:' num2str(nanmean(m1_pk_time))]);
    disp(['M1 sem: ' num2str(nanstd(m1_pk_time)/sqrt(num_m1_pk_times))]);
    disp(['M1 std:' num2str(nanstd(m1_pk_time))]);
    fprintf('\n\n');
    diary off;
        
    nexttile;
    data = [m1_pk_time, ca1_pk_time];
    labels = [repmat({'M1'}, 1, length(m1_pk_time)), repmat({'CA1'}, 1, length(ca1_pk_time))];
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data, labels, 'GroupOrder', {'M1', 'CA1'}, ViolinOpts);

    %violins(1).ViolinColor = {'k'};
    %violins(2).ViolinColor = {'g'};
    violins(1).ViolinColor = {Multi_func.m1_color};
    violins(2).ViolinColor = {Multi_func.ca1_color};

    % Originally were histograms
    %histogram(ca1_pk_time, 100, 'DisplayName', 'CA1'); 
    %hold on;
    %histogram(m1_pk_time, 100, 'DisplayName', 'M1');
    title([num2str(f_stim) ' Vm onset time to peak violins'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Onset_Peak_Time_Vm_CA1_vs_M1.pdf']);
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Onset_Peak_Time_Vm_CA1_vs_M1.png']);

%% Time to peak firing rate between regions
stats_log = [figure_path 'Small_Res' f 'Fr_pop_onset_ca1_vs_m1_times'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;
test_regions = {'r_CA1', 'r_M1'};
stims = fieldnames(region_data.r_M1)';
figure('Position', [140 150 , 1000, 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 15, 20]);
for f_stim = stims
    f_stim = f_stim{1};
    
    % Store the time to peak values for both regions
    ca1_pk_time = [];
    m1_pk_time = [];
    ca1_fr = [];
    m1_fr = [];
    
    % Grab just the first 100 ms onset trace times
    timeline = nanmean(stim_data.trace_timestamps, 2)*1000;
    onset_idx = find(timeline >=0 & timeline <= 100); % Only grab the first 100 ms
    timeline = timeline(onset_idx);
    ca1_fr = region_data.r_CA1.(f_stim).neuron_srate_3(onset_idx, :)./region_data.r_CA1.(f_stim).neuron_spike_amp;
    m1_fr = region_data.r_M1.(f_stim).neuron_srate_3(onset_idx, :)./region_data.r_M1.(f_stim).neuron_spike_amp;
    
    for i=1:size(ca1_fr, 2)
        peak_idx = find(ca1_fr(:, i) == max(ca1_fr(:, i)));
        if length(peak_idx) > 1
            peak_idx = peak_idx(1);
        end
        if timeline(peak_idx) == 0
            continue;
        end
        ca1_pk_time(end + 1) = timeline(peak_idx);
    end

    for i=1:size(m1_fr, 2)
        peak_idx = find(m1_fr(:, i) == max(m1_fr(:, i)));
        if length(peak_idx) > 1
            peak_idx = peak_idx(1);
        end
        if timeline(peak_idx) == 0
            continue;
        end
        m1_pk_time(end + 1) = timeline(peak_idx);
    end

    diary on;
    disp(['Time to peak ' f_stim]);
    [p, h, stats] = ranksum(m1_pk_time, ca1_pk_time)
    num_ca1_pk_times = length(ca1_pk_time);
    num_m1_pk_times = length(m1_pk_time);
    disp(['CA1 num traces ' num2str(num_ca1_pk_times)]);
    disp(['CA1 mean:' num2str(nanmean(ca1_pk_time))]);
    disp(['CA1 sem: ' num2str(nanstd(ca1_pk_time)/sqrt(num_ca1_pk_times))]);
    disp(['CA1 std:' num2str(nanstd(ca1_pk_time))]);
    disp(['M1 num traces ' num2str(num_m1_pk_times)]);
    disp(['M1 mean:' num2str(nanmean(m1_pk_time))]);
    disp(['M1 sem: ' num2str(nanstd(m1_pk_time)/sqrt(num_m1_pk_times))]);
    disp(['M1 std:' num2str(nanstd(m1_pk_time))]);
    fprintf('\n\n');
    diary off;
        
    nexttile;
    data = [m1_pk_time, ca1_pk_time];
    labels = [repmat({'M1'}, 1, length(m1_pk_time)), repmat({'CA1'}, 1, length(ca1_pk_time))];
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data, labels, 'GroupOrder', {'M1', 'CA1'}, ViolinOpts);

    %violins(1).ViolinColor = {'k'};
    %violins(2).ViolinColor = {'g'};
    violins(1).ViolinColor = {Multi_func.m1_color};
    violins(2).ViolinColor = {Multi_func.ca1_color};

    % Originally were histograms
    %histogram(ca1_pk_time, 100, 'DisplayName', 'CA1'); 
    %hold on;
    %histogram(m1_pk_time, 100, 'DisplayName', 'M1');
    title([num2str(f_stim) ' Vm onset time to peak violins'], 'Interpreter', 'none');
end
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Onset_Peak_Time_FR_CA1_vs_M1.pdf']);
saveas(gcf, [figure_path 'Small_Res' f 'Violin_Onset_Peak_Time_FR_CA1_vs_M1.png']);
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
