clear all;
close all;
clc;
f = filesep;

%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end

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

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

% This is mostly for debugging purposes
display_names = 0;

%% END Modification

%% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
    %save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'ca1_data.mat'];
end
%Load the data
load(save_all_data_file);

%% Calculate the Vm modulation stats from onset
wind_dist = 100/1000; % numerator in ms, window from onset to include for baseline and stim comparison
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        Vm_trans_mod_stats = struct;

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            
            % Calculate the average amplitude Vm for each trial
            Vm_mean_pre = mean(popul_data.all_trial_SubVm{nr}(base_idx, :), 1);
            Vm_mean_stim = mean(popul_data.all_trial_SubVm{nr}(stim_idx, :), 1);

            [p_sign, ~, stats] = signrank(Vm_mean_pre, Vm_mean_stim);
            
            % Check the modulation level
            if p_sign < 0.05 && mean(Vm_mean_stim - Vm_mean_pre) > 0
                mod = 1;
            elseif p_sign < 0.05 && mean(Vm_mean_stim - Vm_mean_pre) < 0
                mod = -1;
            else
                mod = 0;
            end

            % Append all of the modulation data
            Vm_trans_mod_stats(nr).p_sign = p_sign;
            Vm_trans_mod_stats(nr).sign_stats = stats;
            Vm_trans_mod_stats(nr).mod = mod;

        end
        region_data.(f_region).(f_stim).Vm_trans_mod_stats = Vm_trans_mod_stats;
    end
end


%% Plot the Vm modulation heatmap and modulated averages based on the transient phase
[heatmap_color] = (cbrewer('div', 'RdBu',500));
heatmap_color(heatmap_color > 1) = 1;
heatmap_color(heatmap_color < 0) = 0;
heatmap_color = flipud(heatmap_color);

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Get the neuron idx of neurons that were modulated
        nr_act = find([popul_data.Vm_trans_mod_stats.mod] > 0);
        nr_sup = find([popul_data.Vm_trans_mod_stats.mod] < 0);
        nr_non = find([popul_data.Vm_trans_mod_stats.mod] == 0);

        % Grab idxs for different parts of the signal
        base_idx = find(mean(popul_data.trace_timestamps, 2) < 0);

        stim_idx = find(mean(popul_data.trace_timestamps, 2) > 0 &...
                    mean(popul_data.trace_timestamps, 2) < 0 + wind_dist);

        % Need to sort across each groups
        [~, act_i] = sort(mean(popul_data.neuron_Vm(stim_idx, nr_act), 1), 'descend');
        [~, sup_i] = sort(mean(popul_data.neuron_Vm(stim_idx, nr_sup), 1), 'descend');
        [~, non_i] = sort(mean(popul_data.neuron_Vm(stim_idx, nr_non), 1), 'descend');
        act_i = nr_act(act_i);
        sup_i = nr_sup(sup_i);
        non_i = nr_non(non_i);

        % Save the trial-averaged neuron Vm
        pop_act_vm = popul_data.neuron_Vm(:, act_i);
        pop_non_vm = popul_data.neuron_Vm(:, non_i);
        pop_sup_vm = popul_data.neuron_Vm(:, sup_i);

        % Zscore all of the trial-averaged neuron Vm
        pop_act_vm = zscore(pop_act_vm, [], 1);
        pop_non_vm = zscore(pop_non_vm, [], 1);
        pop_sup_vm = zscore(pop_sup_vm, [], 1);

        % Zscore the trial-averaged neuron Vm using baseline mean and standard deviation
        %pop_act_vm = (pop_act_vm - mean(pop_act_vm(base_idx, :), 1))./std(pop_act_vm(base_idx, :), 0, 1);
        %
        %pop_non_vm = (pop_non_vm - mean(pop_non_vm(base_idx, :), 1))./std(pop_non_vm(base_idx, :), 0, 1);
        %
        %pop_sup_vm = (pop_sup_vm - mean(pop_sup_vm(base_idx, :), 1))./std(pop_sup_vm(base_idx, :), 0, 1);

        % Baseline subtract the vm data
        pop_act_vm = pop_act_vm - mean(pop_act_vm(base_idx, :), 1);
        pop_non_vm = pop_non_vm - mean(pop_non_vm(base_idx, :), 1);
        pop_sup_vm = pop_sup_vm - mean(pop_sup_vm(base_idx, :), 1);

        % Construct Vm heatmap
        vm_heatmap = pop_act_vm;
        vm_heatmap = horzcat_pad(vm_heatmap, pop_non_vm);
        vm_heatmap = horzcat_pad(vm_heatmap, pop_sup_vm)';

        % Calculate the population average of each groups
        act_Vm_avg = mean(pop_act_vm, 2);
        act_Vm_sem = std(pop_act_vm, 0, 2)./sqrt(size(pop_act_vm, 2));
        non_Vm_avg = mean(pop_non_vm, 2);
        non_Vm_sem = std(pop_non_vm, 0, 2)./sqrt(size(pop_non_vm, 2));
        sup_Vm_avg = mean(pop_sup_vm, 2);
        sup_Vm_sem = std(pop_sup_vm, 0, 2)./sqrt(size(pop_sup_vm, 2));

        % Perform the plotting
        figure('Position', [0, 0, 800, 1000]);
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);

        % Plot the heatmap
        nexttile([2, 1]);
        timeline = mean(popul_data.trace_timestamps, 2);
        imagesc('XData', timeline, 'YData', 1:size(vm_heatmap, 1), 'CData', vm_heatmap);
        hold on;
        yline(0.5 + [length(act_i), length(act_i) + length(non_i), length(act_i) + length(non_i) + length(sup_i)]);
        hold on;
        xline([0 1]);
        hold on;
        Multi_func.set_default_axis(gca);

        %DEBUG
        if display_names == 1
            i=1;
            for nr=[act_i, non_i, sup_i]
                text(2.2, i, [num2str(popul_data.Vm_mod_stats(nr).mod) ' '...
                            popul_data.neuron_name{nr}], ...
                            'Interpreter', 'none');
                i = i+1;
            end
            xlim([-1, 4.5]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end
        ylim([0.5 size(vm_heatmap, 1) + 0.5]);
        xlabel('Time from onset (S)');

        colormap(heatmap_color);
        c = colorbar;
        c.Label.String = 'Vm (z-scored)';

        % Plot the population average for each group activated
        nexttile([1 1]);

        fill_h = fill([timeline; flip(timeline)], ...
                [non_Vm_avg + non_Vm_sem; flipud(non_Vm_avg - non_Vm_sem)], ...
                Multi_func.non_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, non_Vm_avg, 'DisplayName', 'Non-modulated', 'Color', Multi_func.non_color);
        hold on;

        fill_h = fill([timeline; flip(timeline)], ...
                [act_Vm_avg + act_Vm_sem; flipud(act_Vm_avg - act_Vm_sem)], ...
                Multi_func.act_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, act_Vm_avg, 'DisplayName', 'Activated', 'Color', Multi_func.act_color);
        hold on;
        
        fill_h = fill([timeline; flip(timeline)], ...
                [sup_Vm_avg + sup_Vm_sem; flipud(sup_Vm_avg - sup_Vm_sem)], ...
                Multi_func.sup_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, sup_Vm_avg, 'DisplayName', 'Suppressed', 'Color', Multi_func.sup_color);
        hold on;

        Multi_func.set_default_axis(gca);
        hold on;
        xline([0 1], 'HandleVisibility', 'off');
        legend('Location', 'northeastoutside');
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        
        % Set x-axis limits
        if display_names == 0
            xlim([-.7 2.05]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end

        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_Vm_mod_plots.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_Vm_mod_plots.pdf']);
    end
end

%% Calculate the Vm modulation stats for sustained period
per_start = Multi_func.sus_ped(1)/1000; % convert to sec
per_stop = Multi_func.sus_ped(2)/1000;

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        Vm_sus_mod_stats = struct;

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - range([per_start, per_stop]) & nr_timeline < 0);
            stim_idx = find(nr_timeline > per_start & nr_timeline < per_stop);
            
            % Calculate the average amplitude Vm for each trial
            Vm_mean_pre = mean(popul_data.all_trial_SubVm{nr}(base_idx, :), 1);
            Vm_mean_stim = mean(popul_data.all_trial_SubVm{nr}(stim_idx, :), 1);

            [p_sign, ~, stats] = signrank(Vm_mean_pre, Vm_mean_stim);
            
            % Check the modulation level
            if p_sign < 0.05 && mean(Vm_mean_stim - Vm_mean_pre) > 0
                mod = 1;
            elseif p_sign < 0.05 && mean(Vm_mean_stim - Vm_mean_pre) < 0
                mod = -1;
            else
                mod = 0;
            end

            % Append all of the modulation data
            Vm_sus_mod_stats(nr).p_sign = p_sign;
            Vm_sus_mod_stats(nr).sign_stats = stats;
            Vm_sus_mod_stats(nr).mod = mod;

        end
        region_data.(f_region).(f_stim).Vm_sus_mod_stats = Vm_sus_mod_stats;
    end
end


%% Calculate modulation of spike rate in the transient period
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        fr_trans_mod_stats = struct;

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            
            % Calculate the average amplitude Vm for each trial
            fr_mean_pre = sum(popul_data.all_trial_spike_rasters{nr}(base_idx, :), 1);
            fr_mean_stim = sum(popul_data.all_trial_spike_rasters{nr}(stim_idx, :), 1);

            [p_sign, ~, stats] = signrank(fr_mean_pre, fr_mean_stim);
            
            % Check the modulation level
            if p_sign < 0.05 && mean(fr_mean_stim - fr_mean_pre) > 0
                mod = 1;
            elseif p_sign < 0.05 && mean(fr_mean_stim - fr_mean_pre) < 0
                mod = -1;
            else
                mod = 0;
            end

            %DEBUG
            %figure;
            %imagesc(popul_data.all_trial_spike_rasters{nr}');
            %hold on;
            %plot(base_idx, repmat(1, length(base_idx), 1), 'r');
            %hold on;
            %plot(stim_idx, repmat(1, length(stim_idx), 1), 'g');
            %legend([num2str(p_sign) ' ' num2str(mean(fr_mean_stim - fr_mean_pre))]);
            %title(popul_data.neuron_name{nr}, 'Interpreter', 'none');

            % Append all of the modulation data
            fr_trans_mod_stats(nr).p_sign = p_sign;
            fr_trans_mod_stats(nr).sign_stats = stats;
            fr_trans_mod_stats(nr).mod = mod;

        end
        region_data.(f_region).(f_stim).fr_trans_mod_stats = fr_trans_mod_stats;
    end
end

%% Plot the firing rate modulation heatmap and modulated averages
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Get the neuron idx of neurons that were modulated
        nr_act = find([popul_data.fr_trans_mod_stats.mod] > 0);
        nr_sup = find([popul_data.fr_trans_mod_stats.mod] < 0);
        nr_non = find([popul_data.fr_trans_mod_stats.mod] == 0);

        % Grab idxs for different parts of the signal
        base_idx = find(mean(popul_data.trace_timestamps, 2) < 0);

        stim_idx = find(mean(popul_data.trace_timestamps, 2) > 0 &...
                    mean(popul_data.trace_timestamps, 2) < 0 + wind_dist);

        % Need to sort across each groups
        [~, act_i] = sort(mean(popul_data.neuron_srate_50(stim_idx, nr_act), 1), 'descend');
        [~, sup_i] = sort(mean(popul_data.neuron_srate_50(stim_idx, nr_sup), 1), 'descend');
        [~, non_i] = sort(mean(popul_data.neuron_srate_50(stim_idx, nr_non), 1), 'descend');
        act_i = nr_act(act_i);
        sup_i = nr_sup(sup_i);
        non_i = nr_non(non_i);

        % Save the trial-averaged neuron firing rate
        
        pop_act_fr = popul_data.neuron_srate_50(:, act_i);
        pop_non_fr = popul_data.neuron_srate_50(:, non_i);
        pop_sup_fr = popul_data.neuron_srate_50(:, sup_i);

        % Zscore all of the trial-averaged neuron firing rate
        pop_act_fr = zscore(pop_act_fr, [], 1);
        pop_non_fr = zscore(pop_non_fr, [], 1);
        pop_sup_fr = zscore(pop_sup_fr, [], 1);

        % Zscore the trial-averaged neuron fr using baseline mean and standard deviation
        % Check for zero std
        %if std(pop_act_fr(base_idx, :), 0, 1) == 0
        %    pop_act_fr = (pop_act_fr - mean(pop_act_fr(base_idx, :), 1));
        %else
        %    pop_act_fr = (pop_act_fr - mean(pop_act_fr(base_idx, :), 1))./std(pop_act_fr(base_idx, :), 0, 1);
        %end

        %if std(pop_non_fr(base_idx, :), 0, 1) == 0
        %    pop_non_fr = (pop_non_fr - mean(pop_non_fr(base_idx, :), 1));
        %else
        %    pop_non_fr = (pop_non_fr - mean(pop_non_fr(base_idx, :), 1))./std(pop_non_fr(base_idx, :), 0, 1);
        %end
        %
        %if std(pop_sup_fr(base_idx, :), 0, 1) == 0
        %    pop_sup_fr = (pop_sup_fr - mean(pop_sup_fr(base_idx, :), 1));
        %else
        %    pop_sup_fr = (pop_sup_fr - mean(pop_sup_fr(base_idx, :), 1))./std(pop_sup_fr(base_idx, :), 0, 1);
        %end
        

        % Baseline subtract the fr data
        pop_act_fr = pop_act_fr - mean(pop_act_fr(base_idx, :), 1);
        pop_non_fr = pop_non_fr - mean(pop_non_fr(base_idx, :), 1);
        pop_sup_fr = pop_sup_fr - mean(pop_sup_fr(base_idx, :), 1);

        % Construct fr heatmap
        fr_heatmap = pop_act_fr;
        fr_heatmap = horzcat_pad(fr_heatmap, pop_non_fr);
        fr_heatmap = horzcat_pad(fr_heatmap, pop_sup_fr)';

        % Calculate the population average of each groups
        act_fr_avg = mean(pop_act_fr, 2);
        act_fr_sem = std(pop_act_fr, 0, 2)./sqrt(size(pop_act_fr, 2));
        non_fr_avg = mean(pop_non_fr, 2);
        non_fr_sem = std(pop_non_fr, 0, 2)./sqrt(size(pop_non_fr, 2));
        sup_fr_avg = mean(pop_sup_fr, 2);
        sup_fr_sem = std(pop_sup_fr, 0, 2)./sqrt(size(pop_sup_fr, 2));

        % Perform the plotting
        figure('Position', [0, 0, 800, 1000]);
        fontsize = 11;
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);

        % Plot the heatmap
        nexttile([2 1]);
        timeline = mean(popul_data.trace_timestamps, 2);
        imagesc('XData', timeline, 'YData', 1:size(fr_heatmap, 1), 'CData', fr_heatmap);
        hold on;
        yline(0.5 + [length(act_i), length(act_i) + length(non_i), length(act_i) + length(non_i) + length(sup_i)]);
        hold on;
        xline([0 1]);
        hold on;
        Multi_func.set_default_axis(gca);

        %DEBUG
        if display_names == 1
            i=1;
            for nr=[act_i, non_i, sup_i]
                text(2.2, i, [num2str(popul_data.fr_trans_mod_stats(nr).mod) ' '...
                            popul_data.neuron_name{nr}], ...
                            'Interpreter', 'none');
                i = i+1;
            end
            xlim([-1, 4.5]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end
        ylim([0.5 size(fr_heatmap, 1) + 0.5]);
        xlabel('Time from onset (S)');


        colormap(heatmap_color);
        c = colorbar;
        caxis([-5 5]);
        c.Label.String = 'spike rate (z-scored)';

        % Plot the population average for each group activated
        nexttile;

        fill_h = fill([timeline; flip(timeline)], ...
                [non_fr_avg + non_fr_sem; flipud(non_fr_avg - non_fr_sem)], ...
                Multi_func.non_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, non_fr_avg, 'DisplayName', 'Non-modulated', 'Color', Multi_func.non_color);
        hold on;

        fill_h = fill([timeline; flip(timeline)], ...
                [act_fr_avg + act_fr_sem; flipud(act_fr_avg - act_fr_sem)], ...
                Multi_func.act_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, act_fr_avg, 'DisplayName', 'Activated', 'Color', Multi_func.act_color);
        hold on;
        
        fill_h = fill([timeline; flip(timeline)], ...
                [sup_fr_avg + sup_fr_sem; flipud(sup_fr_avg - sup_fr_sem)], ...
                Multi_func.sup_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, sup_fr_avg, 'DisplayName', 'Suppressed', 'Color', Multi_func.sup_color);
        hold on;

        Multi_func.set_default_axis(gca);
        hold on;
        xline([0 1], 'HandleVisibility', 'off');
        legend('Location', 'northeastoutside');
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        
        % Set x-axis limits
        if display_names == 0
            xlim([-.7 2.05]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end

        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_fr_mod_plots.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_fr_mod_plots.pdf']);
    end
end


%% Calculate modulation of spike rate in the sustained period
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        fr_sus_mod_stats = struct;

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            
            % Calculate the average amplitude Vm for each trial
            fr_mean_pre = sum(popul_data.all_trial_spike_rasters{nr}(base_idx, :), 1);
            fr_mean_stim = sum(popul_data.all_trial_spike_rasters{nr}(stim_idx, :), 1);

            [p_sign, ~, stats] = signrank(fr_mean_pre, fr_mean_stim);
            
            % Check the modulation level
            if p_sign < 0.05 && mean(fr_mean_stim - fr_mean_pre) > 0
                mod = 1;
            elseif p_sign < 0.05 && mean(fr_mean_stim - fr_mean_pre) < 0
                mod = -1;
            else
                mod = 0;
            end

            %DEBUG
            figure;
            imagesc(popul_data.all_trial_spike_rasters{nr}');
            hold on;
            plot(base_idx, repmat(1, length(base_idx), 1), 'r');
            hold on;
            plot(stim_idx, repmat(1, length(stim_idx), 1), 'g');
            legend([num2str(p_sign) ' ' num2str(mean(fr_mean_stim - fr_mean_pre))]);
            title(popul_data.neuron_name{nr}, 'Interpreter', 'none');

            % Append all of the modulation data
            fr_sus_mod_stats(nr).p_sign = p_sign;
            fr_sus_mod_stats(nr).sign_stats = stats;
            fr_sus_mod_stats(nr).mod = mod;

        end
        region_data.(f_region).(f_stim).fr_sus_mod_stats = fr_sus_mod_stats;
    end
end

%% Calculate single cell PLV as well as shuffled distribution
num_iter = 500; 
wind_dist = 1000/1000; %ms
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Stim freq for PLV
        plv_freq = str2num(erase(f_stim, 'f_'))
        avg_Fs = mean(popul_data.framerate, 'omitnan');
        plv_mod_stats = struct;

        figure('Position', [0, 0, 800, 1000]);
        tiledlayout(length(popul_data.neuron_name), 2, 'TileSpacing', 'compact', 'Padding', 'compact');

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);

            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, plv_freq, avg_Fs)));
            
            % Function that waits for a function and matrix and then returns a function waiting for a column value
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            % Aply the filtering function with the whole trial matrix, so all is left is waiting for a column number
            partial_apply = applyFunToColi(filt_trial, popul_data.all_trial_SubVm{nr});
            % Iterate through each column trial and concatenate all o fthe results
            vm_phases = arrayfun(partial_apply, [1:size(popul_data.all_trial_SubVm{nr}, 2)]' , 'UniformOutput', false); 
            % 
            vm_phases = cat(3, vm_phases{:});

            % Calculate stimulation raster with dimensions of the whole trace
            nr_stim_time = reshape(popul_data.stim_timestamps(:, nr), 1, []);
            nr_frame_time = popul_data.trace_timestamps(:, nr);
            diffs = abs(nr_stim_time - nr_frame_time);
            [~, stim_idx_i] = min(diffs, [], 1);
            
            dbs_raster = zeros(size(nr_frame_time));
            dbs_raster(stim_idx_i) = 1;
            dbs_rasters = repmat(dbs_raster, 1, size(vm_phases, 3));

            stim_rasters = zeros(size(dbs_rasters));
            stim_rasters(stim_idx) = dbs_rasters(stim_idx);
            
            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, stim_rasters, 0, 10);             
            obs_PLV = PLV;
            obs_PLV2 = PLV2;
            obs_norm_vecs= norm_vecs;

            % Find the upper limit for randomization
            filt_time = nr_frame_time(nr_stim_time(1) > nr_frame_time);
            % Compute absolute differences
            diffs = abs(filt_time - nr_stim_time(1));
            % Find index of minimum difference in filtered_array
            [minDiff, minIndex] = min(diffs);
            % Find index of corresponding element in original array
            last_index = find(nr_frame_time == filt_time(minIndex), 1);

            % Shuffled PLV data
            shuf_plv = [];
            shuf_plv_adj = [];

            % Shuffle the start of the dbs timepoints and recalculate the PLV
            % Need to randomize where the DBS points are along the rasters
            for i=1:num_iter
                % Randomly select stim onset time
                onset_rand = randi(last_index, 1, size(vm_phases, 3));
                
                % Reset the starting stim_idx for each trial and keep the spacing between indices the same
                rand_stim_idx = repmat(stim_idx_i' - stim_idx_i(1), 1, size(vm_phases, 3)) + onset_rand;
                col_ind_mat = repmat(1:size(rand_stim_idx, 2), size(rand_stim_idx, 1), 1);
                stim_rasters = zeros(size(dbs_rasters));   
                stim_rasters(sub2ind(size(stim_rasters),...
                                     rand_stim_idx(:),...
                                     col_ind_mat(:))) = 1;
            
                [sh_PLV, sh_PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, stim_rasters, 0, 10);            
                shuf_plv(end + 1) = sh_PLV;
                shuf_plv_adj(end + 1) = sh_PLV2;
            end
        
            % Calculate the 95th percentile
            high_prc = prctile(shuf_plv_adj, 95);

            % Plot the trial-averaged Vm
            nexttile;
            plot(nr_timeline, popul_data.neuron_Vm(:, nr));

            % plot the shuffled distribution
            nexttile;
            histogram(shuf_plv_adj, 1000);
            hold on;
            xline(high_prc, '-b');
            hold on;
            Multi_func.set_default_axis(gca);
            
            % Check if significant
            if obs_PLV2 > high_prc
                xline(obs_PLV2, '-g');
                plv_mod_stats(nr).mod = 1;
            else
                xline(obs_PLV2, '-r');
                plv_mod_stats(nr).mod = -1;
            end

            sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');

            % Store the PLV data into the neuron structure
            plv_mod_stats(nr).obs_PLV = obs_PLV;
            plv_mod_stats(nr).obs_PLV2 = obs_PLV2;
            plv_mod_stats(nr).shuf_PLV = sh_PLV;
            plv_mod_stats(nr).shuf_PLV2 = sh_PLV2;

        end %neuron loop
        % Add plv mod stats to structure
        region_data.(f_region).(f_stim).plv_mod_stats = plv_mod_stats;
    end % Stim freq loop
end % Region loop

%% Plot stuff based on stim-Vm PLV
%plot_mode = 'pow'; % Use stimulation frequency power spectra for each neuron
%plot_mode = 'Vm'; % Use the trial averaged Vm
plot_mode = 'pulse'; % use Pulse-triggered average
freqs = Multi_func.entr_freqs;
extra_trace = 3;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        avg_Fs = mean(popul_data.framerate, 'omitnan');

        cur_freq = str2num(erase(f_stim, 'f_'));

        % Get the neuron idx of neurons that were modulated
        nr_etrain = find([popul_data.plv_mod_stats.mod] > 0);
        nr_non = find([popul_data.plv_mod_stats.mod] < 0);

        % Grab idxs for different parts of the signal
        base_idx = find(mean(popul_data.trace_timestamps, 2) < 0);

        stim_idx = find(mean(popul_data.trace_timestamps, 2) > 0 &...
                    mean(popul_data.trace_timestamps, 2) < 0 + wind_dist);

        % Need to sort across each groups
        [~, etrain_i] = sort([popul_data.plv_mod_stats(nr_etrain).obs_PLV2], 'descend');
        [~, non_i] = sort([popul_data.plv_mod_stats(nr_non).obs_PLV2], 'descend');
        etrain_i = nr_etrain(etrain_i);
        non_i = nr_non(non_i);

        switch plot_mode
            case 'pow'

                % Calulcate the trial-averaged power spectra at the given stimulation frequency for all neurons
                pop_entrain_pow = cellfun(@(x) mean(abs(x(cur_freq, :, :)), 3)',...
                                            popul_data.neuron_hilbfilt(etrain_i),...
                                            'UniformOutput', false); 

                pop_non_pow = cellfun(@(x) mean(abs(x(cur_freq, :, :)), 3)',...
                                            popul_data.neuron_hilbfilt(non_i),...
                                            'UniformOutput', false); 
                % zscore everything
                pop_entrain_pow = zscore([pop_entrain_pow{:}], [], 1);
                pop_non_pow = zscore([pop_non_pow{:}], [], 1);

                % Baseline subtract power data
                pop_entrain_pow = pop_entrain_pow - mean(pop_entrain_pow(base_idx, :), 1);
                pop_non_pow = pop_non_pow - mean(pop_non_pow(base_idx, :), 1);

                % Construct heatmap
                pow_heatmap = pop_entrain_pow;
                pow_heatmap = horzcat_pad(pow_heatmap, pop_non_pow)';
            
            case 'Vm'
                % Save the trial-averaged neuron Vm
                pop_etrain_vm = popul_data.neuron_Vm(:, etrain_i);
                pop_non_vm = popul_data.neuron_Vm(:, non_i);

                % Zscore all of the trial-averaged neuron Vm
                pop_etrain_vm = zscore(pop_etrain_vm, [], 1);
                pop_non_vm = zscore(pop_non_vm, [], 1);
                
                % Zscore the trial-averaged neuron Vm using baseline mean and standard deviation
                %pop_etrain_vm = (pop_etrain_vm - mean(pop_etrain_vm(base_idx, :), 1))./std(pop_etrain_vm(base_idx, :), 0, 1);
                %
                %pop_non_vm = (pop_non_vm - mean(pop_non_vm(base_idx, :), 1))./std(pop_non_vm(base_idx, :), 0, 1);
                %
                %pop_sup_vm = (pop_sup_vm - mean(pop_sup_vm(base_idx, :), 1))./std(pop_sup_vm(base_idx, :), 0, 1);

                % Baseline subtract the vm data
                pop_etrain_vm = pop_etrain_vm - mean(pop_etrain_vm(base_idx, :), 1);
                pop_non_vm = pop_non_vm - mean(pop_non_vm(base_idx, :), 1);

                % Construct Vm heatmap
                vm_heatmap = pop_etrain_vm;
                vm_heatmap = horzcat_pad(vm_heatmap, pop_non_vm)';

                % Calculate the population average of each groups
                etrain_Vm_avg = mean(pop_etrain_vm, 2);
                etrain_Vm_sem = std(pop_etrain_vm, 0, 2)./sqrt(size(pop_etrain_vm, 2));
                non_Vm_avg = mean(pop_non_vm, 2);
                non_Vm_sem = std(pop_non_vm, 0, 2)./sqrt(size(pop_non_vm, 2));
            
            case 'pulse'
                pop_pulse_vm = [];

                % Need to loop through each neuron and calculate its pulse triggered average
                for nr=1:length(popul_data.neuron_name)
                    % Calculate the number of trace idxs between pulses
                    nr_avg_pulse_width_time = mean(diff(data_bystim.(f_stim).stim_timestamps(:, nr) ), 'omitnan');
                    nr_avg_trace_time = mean(diff(data_bystim.(f_stim).trace_timestamps(:, nr) ), 'omitnan');
                    trace_width = ceil(nr_avg_pulse_width_time./nr_avg_trace_time);
                        
                    % Grab a single trial's pulse window idxs
                    get_follow_trace_idx = @(pulse_time) ...
                        find(pulse_time < popul_data.trace_timestamps(:, nr), 1);
                    get_width_idx = @(pulse_time) [get_follow_trace_idx(pulse_time) - extra_trace:...
                                    get_follow_trace_idx(pulse_time) + trace_width + extra_trace]';
                    pulse_trial_wind_idxs = arrayfun(get_width_idx,...
                        popul_data.stim_timestamps(:, nr)', 'UniformOutput', false);

                    % Grab the window from the trial, could probably do baseline subtraction from here
                    apply_trial = @(trial) @(indices) trial(indices) - trial(indices(1 + extra_trace));
                    index_indTrial_pulses = @(trial) cell2mat(cellfun(apply_trial(trial), pulse_trial_wind_idxs,...
                        'UniformOutput', false));

                    % Convert matrix into cell array by columns
                    all_trial = popul_data.all_trial_SubVm{nr};
                    all_trial_SubVm_bycol = mat2cell(all_trial, size(all_trial, 1),...
                        ones(1, size(all_trial, 2) ) );
                    
                    nr_norm_vm = cellfun(index_indTrial_pulses, all_trial_SubVm_bycol, 'UniformOutput', false);
                    pop_pulse_vm(:, end + 1) = mean(cell2mat(nr_norm_vm), 2);
                end
            
                % Separate all pulse stuff between entrained and not entrained neurons
                pop_etrain_pulse = pop_pulse_vm(:, etrain_i);
                pop_non_pulse = pop_pulse_vm(:, non_i);

                % Zscore the pulse-triggered averages
                %pop_etrain_pulse = zscore(pop_etrain_pulse, [], 1);
                %pop_non_pulse = zscore(pop_non_pulse, [], 1);

                pulse_heatmap = pop_etrain_pulse;
                pulse_heatmap = horzcat_pad(pulse_heatmap, pop_non_pulse)';

        end

        % Compute the PLV population average between entrained and non entrained neurons
        etrain_plv_avg = mean(popul_data.stim_dbsvm_plvs_adj(:, etrain_i), 2);
        etrain_plv_sem = std(popul_data.stim_dbsvm_plvs_adj(:, etrain_i), 0, 2)./sqrt(length(etrain_i));
        
        non_plv_avg = mean(popul_data.stim_dbsvm_plvs_adj(:, non_i), 2);
        non_plv_sem = std(popul_data.stim_dbsvm_plvs_adj(:, non_i), 0, 2)./sqrt(length(non_i));

        % Perform the plotting
        figure('Position', [0, 0, 800, 1000]);
        tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);

        switch plot_mode
            case 'pow'
                nexttile;
                timeline = mean(popul_data.trace_timestamps, 2);
                
                imagesc('XData', timeline, 'YData', 1:size(pow_heatmap, 1), 'CData', pow_heatmap);
                hold on;
                yline(0.5 + [length(etrain_i), length(etrain_i) + length(non_i)]);
                hold on;
                xline([0 1]);
                hold on;
                Multi_func.set_default_axis(gca);

                 %DEBUG
                if display_names == 1
                    i=1;
                    for nr=[etrain_i, non_i]
                        text(2.2, i, [num2str(popul_data.plv_mod_stats(nr).mod) ' '...
                                    popul_data.neuron_name{nr}], ...
                                    'Interpreter', 'none');
                        i = i+1;
                    end
                    xlim([-1, 4.5]);
                else
                    xlim([min(timeline) - .01, max(timeline)]);
                end
                ylim([0.5 size(pow_heatmap, 1) + 0.5]);
                xlabel('Time from onset (S)');

                colormap(heatmap_color);
                c = colorbar;
                c.Label.String = 'Power Spectra (z-scored)';

            case 'Vm'
                %---------- Plot the population Vm heatmap
                nexttile;
                timeline = mean(popul_data.trace_timestamps, 2);

                imagesc('XData', timeline, 'YData', 1:size(vm_heatmap, 1), 'CData', vm_heatmap);
                hold on;
                yline(0.5 + [length(etrain_i), length(etrain_i) + length(non_i)]);
                hold on;
                xline([0 1]);
                hold on;
                Multi_func.set_default_axis(gca);

                 %DEBUG
                if display_names == 1
                    i=1;
                    for nr=[etrain_i, non_i, sup_i]
                        text(2.2, i, [num2str(popul_data.plv_mod_stats(nr).mod) ' '...
                                    popul_data.neuron_name{nr}], ...
                                    'Interpreter', 'none');
                        i = i+1;
                    end
                    xlim([-1, 4.5]);
                else
                    xlim([min(timeline) - .01, max(timeline)]);
                end
                ylim([0.5 size(vm_heatmap, 1) + 0.5]);
                xlabel('Time from onset (S)');

                colormap(heatmap_color);
                c = colorbar;
                c.Label.String = 'Vm (z-scored)';
                %---------- End Vm heatmap plotting
            case 'pulse'
                nexttile;
                
                timeline = ([0:size(pulse_heatmap, 2)]- extra_trace)*1000./avg_Fs;

                imagesc('XData', timeline, 'YData', 1:size(pulse_heatmap, 1), 'CData', pulse_heatmap);
                hold on;
                yline(0.5 + [length(etrain_i), length(etrain_i) + length(non_i)]);
                hold on;
                xline([0 round(1000/cur_freq)]);
                hold on;
                Multi_func.set_default_axis(gca);

                 %DEBUG
                if display_names == 1
                    i=1;
                    for nr=[etrain_i, non_i, sup_i]
                        text(2.2, i, [num2str(popul_data.plv_mod_stats(nr).mod) ' '...
                                    popul_data.neuron_name{nr}], ...
                                    'Interpreter', 'none');
                        i = i+1;
                    end
                    xlim([-1, 4.5]);
                else
                    xlim([min(timeline) - .01, max(timeline)]);
                end
                ylim([0.5 size(pulse_heatmap, 1) + 0.5]);
                xlabel('Time from onset (S)');

                colormap(heatmap_color);
                c = colorbar;
                clim = caxis;
                caxis([-1 1]*range(clim)/2);
                c.Label.String = 'Vm (z-scored)';

        end

        % Plot the population average for each group entrainment
        nexttile;
        
        fill_h = fill([freqs'; flip(freqs)'], ...
                [non_plv_avg + non_plv_sem; flipud(non_plv_avg - non_plv_sem)], ...
                Multi_func.non_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        fill_h.EdgeAlpha = 0;
        hold on;
        plot(freqs', non_plv_avg, 'DisplayName', 'Non-entrained', 'Color', Multi_func.non_color);

        hold on;
        fill_h = fill([freqs'; flip(freqs)'], ...
                [etrain_plv_avg + etrain_plv_sem; flipud(etrain_plv_avg - etrain_plv_sem)], ...
                Multi_func.act_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(freqs', etrain_plv_avg, 'DisplayName', 'Entrained', 'Color', Multi_func.act_color);
        hold on;

        set(gca, 'XScale', 'log');
        Multi_func.set_default_axis(gca);
        legend();
        ylabel('DBS-Vm PLV^2');
        xlabel('Frequency (Hz)');

        % Plot each neuron's PLV for all frequencies
        nexttile;
            
        % Zscore values
        %pop_etrain_plv = zscore(popul_data.stim_dbsvm_plvs_adj(:, etrain_i), [], 1);
        %pop_non_plv = zscore(popul_data.stim_dbsvm_plvs_adj(:, non_i), [], 1);

        % regular values
        pop_etrain_plv = popul_data.stim_dbsvm_plvs_adj(:, etrain_i);
        pop_non_plv = popul_data.stim_dbsvm_plvs_adj(:, non_i);

        plv_heatmap = pop_etrain_plv;
        plv_heatmap = horzcat_pad(plv_heatmap, pop_non_plv)';
        
        % Maybe linear interpolate data along the logged axis?
        [X, Y] = meshgrid(freqs, 1:size(plv_heatmap, 1));
        surface(X, Y, plv_heatmap, ...
            'CDataMapping', 'scale', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
        hold on;
        yline([length(etrain_i), length(etrain_i) + length(non_i)]);

        colormap(heatmap_color);
        c = colorbar;
        caxis([-1 1]);
        xlim([0 200]);
        ylim([1 size(plv_heatmap, 1) ]);
        xlabel('Frequency (Hz)');
        set(gca, 'XScale', 'log');
        c.Label.String = 'DBS-Vm PLV (Z-scored)';

        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        
        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_' num2str(wind_dist) '_stimVm_plv_mod_plots.png']);
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_' num2str(wind_dist) '_stimVm_plv_mod_plots_mode' plot_mode '.pdf']);
    end
end



%% Display the number of modulated stuff
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        disp([f_region ' ' f_stim ' ' num2str(length([popul_data.Vm_mod_stats.mod]))]);

        % Display the modulation of Vm
        disp('Vm');
        disp(['Activated: ' num2str(length(find([popul_data.Vm_mod_stats.mod] > 0 )))]);
        disp(['Non-modulated: ' num2str(length(find([popul_data.Vm_mod_stats.mod] == 0 )))]);
        disp(['Suppression: ' num2str(length(find([popul_data.Vm_mod_stats.mod] < 0 )))]);
        fprintf('\n');
        
        % Display the modulation of Firing Rate
        disp('Firing Rate');
        disp(['Activated: ' num2str(length(find([popul_data.fr_mod_stats.mod] > 0 )))]);
        disp(['Non-modulated: ' num2str(length(find([popul_data.fr_mod_stats.mod] == 0 )))]);
        disp(['Suppression: ' num2str(length(find([popul_data.fr_mod_stats.mod] < 0 )))]);
        fprintf('\n');
    
        % Display the modulation of Firing Rate
        disp('PLV');
        disp(['Entrained: ' num2str(length(find([popul_data.plv_mod_stats.mod] > 0 )))]);
        disp(['Non-entrained: ' num2str(length(find([popul_data.plv_mod_stats.mod] < 0 )))]);
        fprintf('\n');

        fprintf('\n');
        fprintf('\n');
    end
end

%% Display the number of neurons per modulation
stats_t = table();
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
       
        % Grab the neuron number of each Vm trans positive modulation
        vm_act_nr = find([popul_data.Vm_trans_mod_stats.mod] > 0);
        vm_non_nr = find([popul_data.Vm_trans_mod_stats.mod] == 0);
        vm_sup_nr = find([popul_data.Vm_trans_mod_stats.mod] < 0);

        % Grab the firing rate modulated trans data
        fr_act_nr = find([popul_data.fr_trans_mod_stats.mod] > 0);
        fr_non_nr = find([popul_data.fr_trans_mod_stats.mod] == 0);
        fr_sup_nr = find([popul_data.fr_trans_mod_stats.mod] < 0);

        stats_t([f_region f_stim], 'Vm Trans Act') = {length(vm_act_nr)};
        stats_t([f_region f_stim], 'Vm Trans Non') = {length(vm_non_nr)};
        stats_t([f_region f_stim], 'Vm Trans Sup') = {length(vm_sup_nr)};
    
        stats_t([f_region f_stim], 'FR Trans Act') = {length(fr_act_nr)};
        stats_t([f_region f_stim], 'FR Trans Non') = {length(fr_non_nr)};
        stats_t([f_region f_stim], 'FR Trans Sup') = {length(fr_sup_nr)};

        %----- Modulated data for sustained period
        % Grab the neuron number of each Vm sus positive modulation
        vm_act_nr = find([popul_data.Vm_sus_mod_stats.mod] > 0);
        vm_non_nr = find([popul_data.Vm_sus_mod_stats.mod] == 0);
        vm_sup_nr = find([popul_data.Vm_sus_mod_stats.mod] < 0);

        % Grab the firing rate modulated sus data
        fr_act_nr = find([popul_data.fr_sus_mod_stats.mod] > 0);
        fr_non_nr = find([popul_data.fr_sus_mod_stats.mod] == 0);
        fr_sup_nr = find([popul_data.fr_sus_mod_stats.mod] < 0);

        stats_t([f_region f_stim], 'Vm sus Act') = {length(vm_act_nr)};
        stats_t([f_region f_stim], 'Vm sus Non') = {length(vm_non_nr)};
        stats_t([f_region f_stim], 'Vm sus Sup') = {length(vm_sup_nr)};
    
        stats_t([f_region f_stim], 'FR sus Act') = {length(fr_act_nr)};
        stats_t([f_region f_stim], 'FR sus Non') = {length(fr_non_nr)};
        stats_t([f_region f_stim], 'FR sus Sup') = {length(fr_sup_nr)};
        

        %--- Determine the entrainment PLV mod stuff
        etrain_nr = find([popul_data.plv_mod_stats.mod] > 0);
        non_nr = find([popul_data.plv_mod_stats.mod] < 0);
        stats_t([f_region f_stim], 'PLV Etrain') = {length(etrain_nr)};
        stats_t([f_region f_stim], 'PLV non-etrain') = {length(non_nr)};
        
        %--- Add total number of neurons
        stats_t([f_region f_stim], 'Total Neurons') = {length([popul_data.plv_mod_stats.mod])};

        % Loop and create percentages
        for row = stats_t.Properties.RowNames'
            stats_t{row{1}, :} = stats_t{row{1}, :}*100 / stats_t{row{1}, 'Total Neurons'};
        end
    end
end
disp(stats_t);
writetable(stats_t, [figure_path 'Neuronwise' f 'modulation_stats.csv'], 'WriteRowNames', true);

%% Display the number of neurons that demonstrated similar modulation between Vm and entrainment
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Grab the neuron number of each Vm positive modulation
        vm_act_nr = find([popul_data.Vm_mod_stats.mod] > 0);
        vm_non_nr = find([popul_data.Vm_mod_stats.mod] == 0);
        vm_sup_nr = find([popul_data.Vm_mod_stats.mod] < 0);

        % Grab the firing rate modulated data
        fr_act_nr = find([popul_data.fr_mod_stats.mod] > 0);
        fr_non_nr = find([popul_data.fr_mod_stats.mod] == 0);
        fr_sup_nr = find([popul_data.fr_mod_stats.mod] < 0);

        % Grab the PLV modulation
        etrain_nr = find([popul_data.plv_mod_stats.mod] > 0);
        non_etrain_nr = find([popul_data.plv_mod_stats.mod] < 0);

        sets = {vm_act_nr, vm_non_nr, vm_sup_nr,...
            fr_act_nr, fr_non_nr, fr_sup_nr, etrain_nr, non_etrain_nr};

        setLabels = {'Vm activated', 'Vm non-modul', 'Vm suppr', 'FR activated',...
            'FR non-modul', 'FR suppr', 'Entrained', 'Non-entrained'};
        DrawVennDiag(length(sets), sets, setLabels);
    end
end
