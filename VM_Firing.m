clear all;
close all;
f = filesep;
warning('off', 'all');
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

% -- Depracated --
% Check if combining all of the regions or not
%if all_regions == 1
%    region_data = Multi_func.combine_regions_old(region_data);
%end

if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

field1 = fieldnames(region_data);
field1 = field1(1);
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');

%%Compact full collective spike rate over time

% Filter neurons based on this variable
nr_pop = 'all_mod';
%nr_pop = 'non_mod';

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 17.56, 3]);
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data =data_bystim.(f_stim);
        
        % Determine specific neuron population to plot
        try
            switch nr_pop
                case 'all_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);
                case 'non_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) == 0);
            end
        catch ME
            disp(ME.message);
        end
        
        % 50 Has been used for the paper
        cur_win_srate = 3;
        timeline = nanmean(popul_data.trace_timestamps, 2);
        cur_srate = mean(popul_data.neuron_srate_3(:, nr_idxs), 2, 'omitnan');
        std_srate = std(popul_data.neuron_srate_3(:, nr_idxs), 0, 2, 'omitnan');
        num_neurons = size(popul_data.neuron_srate_3(:, nr_idxs), 2);
        %num_points = size(data_bystim.(f_stim).neurons_srate_3, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile;
    
        % Plot the transient and sustained rectangle squares
        trans_fill = fill([Multi_func.trans_ped(1), Multi_func.trans_ped(2), ...
                        Multi_func.trans_ped(2), Multi_func.trans_ped(1)  ]./1000, ...
                        [min(cur_srate - sem_srate), min(cur_srate - sem_srate), ...
                        max(cur_srate + sem_srate), max(cur_srate + sem_srate) ], ...
                        Multi_func.trans_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(trans_fill);
        trans_fill.EdgeAlpha = 0;

        hold on;
        sus_fill = fill([Multi_func.sus_ped(1), Multi_func.sus_ped(2), ...
                        Multi_func.sus_ped(2), Multi_func.sus_ped(1)  ]./1000, ...
                        [min(cur_srate - sem_srate), min(cur_srate - sem_srate), ...
                        max(cur_srate + sem_srate), max(cur_srate + sem_srate) ], ...
                        Multi_func.sus_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(sus_fill);
        sus_fill.EdgeAlpha = 0;
        hold on;

        fill_h = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;
 
        %DEBUG find the bottom plot of the points
        if strcmp(f_region, 'r_V1') == 1 &&  strcmp(f_stim, 'f_40') == 1
        %    min_val = min(cur_srate - sem_srate);
        %    yline(min_val, 'r|', 'LineWidth', 4);
        %    y.Limits
        end

        % Some axis stuff
        a = gca; y = a.YAxis;
        Multi_func.set_default_axis(gca);

 
        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1 && cur_win_srate == 50
            Multi_func.plot_dbs_bar([0, 1], 7, [f_stim(3:end) ' Hz DBS']);
            y.Limits = [-2 8];
            Multi_func.set_spacing_axis(y, 2, 1);
        elseif strcmp(f_region, 'r_M1') == 1 && cur_win_srate == 50
            Multi_func.plot_dbs_bar([0, 1], 4, [f_stim(3:end) ' Hz DBS']);
            y.Limits = [-2 5];
            Multi_func.set_spacing_axis(y, 2, 1);
        elseif strcmp(f_region, 'r_V1') == 1 && cur_win_srate == 50
            Multi_func.plot_dbs_bar([0, 1], 18, [f_stim(3:end) ' Hz DBS']);
            y.Limits = [-5 20];
            Multi_func.set_spacing_axis(y, 5, 1);
        
        elseif strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_40') == 1
            min_val = min(cur_srate - sem_srate);
            y.Limits(1) = min_val;
        end

        if cur_win_srate == 3
            Multi_func.plot_dbs_bar([0, 1], y.Limits(2), [f_stim(3:end) ' Hz DBS']);
        end
        


        xlabel('Time from stim onset (s)');
        
        % Plot timeline depending on which chopping was done
        if exclude_200ms
            xlim([-0.85 2.05]);
        else
            xlim([-1 2.05]);
        end
        ylabel('Firing Rate Change (Hz)');
        %title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Average Spike rate'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate' num2str(cur_win_srate) '_' nr_pop '.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate' num2str(cur_win_srate) '_' nr_pop '.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.eps']);
end

%% Compact population Vm

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
%nr_pop = 'etrain';
%nr_pop = 'non';
nr_pop = 'all_mod';
%nr_pop = 'non_mod';

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 17.56, 2.80]);
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
                case 'all_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);
                case 'non_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) == 0);
            end
        catch ME
            disp(ME.message);
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        norm_vms = popul_data.neuron_RawVm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        % Baseline subtract each neuron's baseline from itself
        base_idx = find(mean(popul_data.trace_timestamps, 2) < 0);
        pop_base = mean(norm_vms(base_idx, :), 1);
        norm_vms = norm_vms - pop_base;
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_RawVm, 1);

        % Print out the number of neurons for each condition
        disp([f_region ' ' f_stim]);
        num_neurons

        nexttile;
        
        % Plot the transient and sustained rectangle squares
        trans_fill = fill([Multi_func.trans_ped(1), Multi_func.trans_ped(2), ...
                        Multi_func.trans_ped(2), Multi_func.trans_ped(1)  ]./1000, ...
                        [min(cur_Vm - sem_Vm), min(cur_Vm - sem_Vm), ...
                        max(cur_Vm + sem_Vm), max(cur_Vm + sem_Vm) ], ...
                        Multi_func.trans_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(trans_fill);
        trans_fill.EdgeAlpha = 0;

        hold on;
        sus_fill = fill([Multi_func.sus_ped(1), Multi_func.sus_ped(2), ...
                        Multi_func.sus_ped(2), Multi_func.sus_ped(1)  ]./1000, ...
                        [min(cur_Vm - sem_Vm), min(cur_Vm - sem_Vm), ...
                        max(cur_Vm + sem_Vm), max(cur_Vm + sem_Vm) ], ...
                        Multi_func.sus_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(sus_fill);
        sus_fill.EdgeAlpha = 0;
        hold on;

        fill_h = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5], 'linewidth', 0.2);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_Vm, 'k', 'LineWidth', 0.3);
        hold on;

        % Grab axis handles
        a = gca; y = a.YAxis; x = a.XAxis;
        Multi_func.set_default_axis(gca);

        if exclude_200ms
            xlim([-0.85 2.05]);
        else
            xlim([-1 2.05]);
        end
        
        % Plot a DBS bar generically
        Multi_func.plot_dbs_bar([0, 1], max(cur_Vm)+ 0.3, [f_stim(3:end) ' Hz DBS']);

        % Determine limits based on what is plotted
        %if strcmp(f_region, 'r_M1') == 1
        %    Multi_func.plot_dbs_bar([0, 1], , [f_stim(3:end) 'Hz DBS']);
        %    y.Limits = [-5 10];
        %    Multi_func.set_spacing_axis(y, 5, 1);
        %elseif strcmp(f_region, 'r_combine') == 1
        %    y.Limits = [-2 8];
        %elseif strcmp(f_region, 'r_V1') == 1
        %    Multi_func.plot_dbs_bar([0, 1], 18, [f_stim(3:end) 'Hz DBS']);
        %    y.Limits = [-10 20];
        %    Multi_func.set_spacing_axis(y, 10, 1);
        %end

        %title(f_stim(3:end), 'Interpreter', 'none');
        ylabel('Normalized Vm');
        xlabel('Time from stim onset (s)');
    end

    %fontsize(7, "points")
    sgtitle([f_region(3:end) ' ' nr_pop ' Average Vm'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Average_Raw_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Average_Raw_Vm.pdf']);
    savefig(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Average_Raw_Vm.fig']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.eps']);
end


%% Violin plots Vm transient and sustained periods
% Determine which neuron population to use
nr_pop = 'all_mod';
%nr_pop = 'non_mod';
%nr_pop = 'vm_act';
%nr_pop = 'vm_sup';
%nr_pop = 'vm_non_mod';

% Struct to identify each group of data
sub_vm_stat_data = struct();
% Make Violin plots on Subthreshold Vm for transient and sustained
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    sub_vm_stat_data.(f_region) = struct();

    figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 8.096, 3]);
    for f_stim=stims'
        f_stim = f_stim{1};
        
        popul_data = data_bystim.(f_stim);
        
        % Determine specific neuron population to plot
        try
            % Store neurons based on their Vm modulation
            nr_trans_act = find([popul_data.Vm_trans_mod_stats.mod] > 0); 
            nr_sus_act = find([popul_data.Vm_sus_mod_stats.mod] > 0);

            nr_trans_sup = find([popul_data.Vm_trans_mod_stats.mod] < 0); 
            nr_sus_sup = find([popul_data.Vm_sus_mod_stats.mod] < 0);

            nr_non = find([popul_data.Vm_trans_mod_stats.mod] == 0 ...
                & [popul_data.Vm_sus_mod_stats.mod] == 0);

            % Consider activated neurons that have sustained activation, or transient activation that do not have sustained suppression
            nr_act = nr_sus_act;
            nr_act = unique([nr_act, setdiff(nr_trans_act, nr_sus_sup)]);
            
            % Consider suppressed neurons that have sustained suppression, or transient suppression without sustained activation
            nr_sup = nr_sus_sup;
            nr_sup = unique([nr_sup, setdiff(nr_trans_sup, nr_sus_act)]);
            
            % Remove the neurons that were non-modulated at all
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_non = nr_non(~ismember(nr_non, non_mod_nr));
 
            switch nr_pop
                case 'all_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);
                case 'non_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) == 0);
                case 'vm_act'
                    nr_idxs = nr_act;
                case 'vm_sup'
                    nr_idxs = nr_sup;
                case 'vm_non_mod'
                    nr_idxs = nr_non;
            end
        catch ME
            disp(ME.message);
        end
        
        % If there is no data, skip plot
        if isempty(nr_idxs)
            continue
        end

        nexttile;

        % Normalize the Vm
        %norm_vms = popul_data.neuron_RawVm(:, nr_idxs)./popul_data.neuron_spike_amp(:, nr_idxs);
        
        % %TODO need to change to use the adjusted baseline for transient and sustained
        %tim_stamp = popul_data.trace_timestamps(:, nr_idxs);
        %[trans_r, trans_c] = find(tim_stamp > Multi_func.trans_ped(1)/1000 & tim_stamp < Multi_func.trans_ped(2)/1000);
        %[sus_r, sus_c] = find(tim_stamp > Multi_func.sus_ped(1)/1000 & tim_stamp < Multi_func.sus_ped(2)/1000);
        %[base_r, base_c] = find(tim_stamp > Multi_func.base_ped(1)/1000 & tim_stamp < Multi_func.base_ped(2)/1000);
        %
        %base_r = unique(base_r);
        %base_c = unique(base_c);
        %trans_r = unique(trans_r);
        %trans_c = unique(trans_c);
        %sus_r = unique(sus_r);
        %sus_c = unique(sus_c);

        % % Find the base values and then reshape the array to maintain length
        %neuron_base_vm = nanmean(norm_vms(base_r, base_c), 1);
        %neuron_trans_vm = nanmean(norm_vms(trans_r, trans_c), 1) - neuron_base_vm;
        %neuron_sus_vm = nanmean(norm_vms(sus_r, sus_c), 1) - neuron_base_vm;

        %num_neurons = size(neuron_trans_vm, 2);
 
        % Store the number of neurons
        num_neurons = length(nr_idxs);
        
        % Grab the Vm from all neurons
        all_vm = popul_data.all_trial_rawVm(nr_idxs);
        all_sp_amp = num2cell(popul_data.neuron_spike_amp(:, nr_idxs));
        all_vm = cellfun(@(vm, sp_amp) vm./sp_amp, all_vm, all_sp_amp, 'UniformOutput', false);

        % Store the neuron timestamps into a cell array for compatibility with cellfun()
        tim_stamp = popul_data.trace_timestamps(:, nr_idxs);
        tim_stamp = mat2cell(tim_stamp, size(tim_stamp, 1), ones(1, size(tim_stamp, 2) ) );

        % Calculate the vm for the given period
        calc_vm = @(nr_vm, idxs) mean(nr_vm(idxs, :), 1);

        % Get the baseline idxs for transient period reflected over the onset into the baseline
        %get_trans_base_idx = @(times) find(times >= 0 - Multi_func.trans_ped(2)/1000 & times < 0 - Multi_func.trans_ped(1)/1000);
        %base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);
        
        % Baseline is 500 ms previous of onset
        get_trans_base_idx = @(times) find(times > Multi_func.base_ped(1)/1000 & times < Multi_func.base_ped(2)/1000);
        base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);
        
        % Get the stim idxs for the transient period
        get_trans_stim_idx = @(times) find(times >= 0 + Multi_func.trans_ped(1)/1000 & times < 0 + Multi_func.trans_ped(2)/1000);
        stim_trans_idx = cellfun(get_trans_stim_idx, tim_stamp, 'UniformOutput', false);

        % Calculate the vm for the transient period
        base_vm = cellfun(calc_vm, all_vm, base_trans_idx, 'UniformOutput', false);
        stim_vm = cellfun(calc_vm, all_vm, stim_trans_idx, 'UniformOutput', false);
        
        pop_trans_vms = cellfun(@(b_vm, s_vm) mean(s_vm - b_vm), base_vm, stim_vm, 'UniformOutput', false);
        pop_trans_vms = cell2mat(pop_trans_vms);


        % Get the baseline idxs for sustained period reflected over the onset into the baseline
        %get_sus_base_idx = @(times) find(times >= 0 - Multi_func.sus_ped(2)/1000 & times < 0 - Multi_func.sus_ped(1)/1000);
        %base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);
        
        % Baseline is the 500 ms previous of onset
        get_sus_base_idx = @(times) find(times > Multi_func.base_ped(1)/1000 & times < Multi_func.base_ped(2)/1000 );
        base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);

        % Get the stim idxs for the sustained period
        get_sus_stim_idx = @(times) find(times >= 0 + Multi_func.sus_ped(1)/1000 & times < 0 + Multi_func.sus_ped(2)/1000);
        stim_sus_idx = cellfun(get_sus_stim_idx, tim_stamp, 'UniformOutput', false);

        % Calculate the firing rate for the sustained period
        base_vm = cellfun(calc_vm, all_vm, base_sus_idx, 'UniformOutput', false);
        stim_vm = cellfun(calc_vm, all_vm, stim_sus_idx, 'UniformOutput', false);

        pop_sus_vms = cellfun(@(b_vm, s_vm) mean(s_vm - b_vm), base_vm, stim_vm, 'UniformOutput', false);
        pop_sus_vms = cell2mat(pop_sus_vms);
        
        % Initial p-value testing for each period
        disp([f_region ' ' f_stim ]);
        % Testing the sign stuff here
        trans_vm_p = signtest(pop_trans_vms)
        sus_vm_p = signtest(pop_sus_vms)

        % Plot violins
        labels = [];
        data = [];
        labels = [repmat({'Trans'}, num_neurons, 1); repmat({'Sus'}, num_neurons, 1)];
        data = [pop_trans_vms'; pop_sus_vms'];
 
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Trans', 'Sus'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.trans_color};
        violins(2).ViolinColor = {Multi_func.sus_color};

        hold on;

        % Plot zero horizontal line
        yline(0);

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        
        % Change axis cosmetics
        a = gca; y = a.YAxis;
        Multi_func.set_default_axis(gca);

        
        y.Limits = [-1.5, 2];
        Multi_func.set_spacing_axis(y, 1, 2);
        % Old way of setting axis limits
        % Determine parameters based on what is being plotted
        %if strcmp(f_region, 'r_combine') == 1
        %    y.Limits = [-2 8];
        %    Multi_func.set_spacing_axis(y, 4, 1);
        %elseif strcmp(f_region, 'r_M1') == 1
        %    y.Limits = [-1.5 2];
        %    %Multi_func.set_spacing_axis(y, 10, 1);
        %elseif strcmp(f_region, 'r_V1') == 1
        %    y.Limits = [-1 2];
        %    Multi_func.set_spacing_axis(y, 1, 1);
        %end

        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Normalized Vm Change');

        sub_vm_stat_data.(f_region).(f_stim).trans_vm.(nr_pop) = pop_trans_vms;
        sub_vm_stat_data.(f_region).(f_stim).sus_vm.(nr_pop) = pop_sus_vms;
    end
    
    sgtitle([f_region(3:end) ' Vm Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin_' nr_pop '.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin_' nr_pop '.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.eps'], 'epsc');
end

%% Violin plots firing rate transient and sustained periods
% Determine which neuron population to use

nr_pop = 'all_mod';
%nr_pop = 'non_mod';
%nr_pop = 'vm_act';
%nr_pop = 'vm_sup';
%nr_pop = 'vm_non_mod';

% Struct to identify each group of data
fr_stat_data = struct();

% Firing Rate transient and sustained
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    fr_stat_data.(f_region) = struct();

    figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 8.096, 3]);
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Determine specific neuron population to plot
        try
            % Store neurons based on their Vm modulation
            nr_trans_act = find([popul_data.Vm_trans_mod_stats.mod] > 0); 
            nr_sus_act = find([popul_data.Vm_sus_mod_stats.mod] > 0);

            nr_trans_sup = find([popul_data.Vm_trans_mod_stats.mod] < 0); 
            nr_sus_sup = find([popul_data.Vm_sus_mod_stats.mod] < 0);

            nr_non = find([popul_data.Vm_trans_mod_stats.mod] == 0 ...
                & [popul_data.Vm_sus_mod_stats.mod] == 0);

            % Consider activated neurons that have sustained activation, or transient activation that do not have sustained suppression
            nr_act = nr_sus_act;
            nr_act = unique([nr_act, setdiff(nr_trans_act, nr_sus_sup)]);
            
            % Consider suppressed neurons that have sustained suppression, or transient suppression without sustained activation
            nr_sup = nr_sus_sup;
            nr_sup = unique([nr_sup, setdiff(nr_trans_sup, nr_sus_act)]);
            
            % Remove the neurons that were non-modulated at all
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_non = nr_non(~ismember(nr_non, non_mod_nr));

            switch nr_pop
                case 'all_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);
                case 'non_mod'
                    nr_idxs = find(sum(popul_data.mod_matrix, 2) == 0);
                case 'vm_act'
                    nr_idxs = nr_act;
                case 'vm_sup'
                    nr_idxs = nr_sup;
                case 'vm_non_mod'
                    nr_idxs = nr_non;
            end
        catch ME
            disp(ME.message);
        end

        % If there is no data, skip plot
        if isempty(nr_idxs)
            continue
        end

        nexttile;

        % Calculate the spike rate stuff with the smoothed continuous rate, transient was different, but sustained was basically the same
        % probably because it averages across a large area
        
%        % Store the number of neurons
%        num_neurons = length(nr_idxs);
%        
%        % Grab the FR from all neurons
%        all_fr = popul_data.neuron_srate_50(:, nr_idxs);
%        all_fr = num2cell(all_fr, 1);
%
%        % Store the neuron timestamps into a cell array for compatibility with cellfun()
%        tim_stamp = popul_data.trace_timestamps(:, nr_idxs);
%        tim_stamp = mat2cell(tim_stamp, size(tim_stamp, 1), ones(1, size(tim_stamp, 2) ) );
%
%        % Calculate the fr for the given period
%        calc_fr = @(nr_fr, idxs) mean(nr_fr(idxs, :), 1);
%
%
%        % Get the baseline idxs for transient period
%        % This is the real one
%        get_trans_base_idx = @(times) find(times >= 0 - Multi_func.trans_ped(2)/1000 & times < 0 - Multi_func.trans_ped(1)/1000);
%        base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);
%        
%        % DEBUGGING checking with 500 ms previous of onset
%        %get_trans_base_idx = @(times) find(times > 0 - (500/1000) & times < 0);
%        %base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);
%        
%        % Get the stim idxs for the transient period
%        get_trans_stim_idx = @(times) find(times >= 0 + Multi_func.trans_ped(1)/1000 & times < 0 + Multi_func.trans_ped(2)/1000);
%        stim_trans_idx = cellfun(get_trans_stim_idx, tim_stamp, 'UniformOutput', false);
%
%        % Calculate the fr for the transient period
%        base_fr = cellfun(calc_fr, all_fr, base_trans_idx, 'UniformOutput', false);
%        stim_fr = cellfun(calc_fr, all_fr, stim_trans_idx, 'UniformOutput', false);
%        
%        pop_trans_frs = cellfun(@(b_fr, s_fr) s_fr - b_fr, base_fr, stim_fr, 'UniformOutput', false);
%        pop_trans_frs = cell2mat(pop_trans_frs);
%
%        % Get the baseline idxs for sustained period
%        %This is the real one
%        get_sus_base_idx = @(times) find(times >= 0 - Multi_func.sus_ped(2)/1000 & times < 0 - Multi_func.sus_ped(1)/1000);
%        base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);
%        
%        % DEBUGGING checking with the 500 ms previous of onset
%        %get_sus_base_idx = @(times) find(times > 0 - 500/1000 & times < 0 );
%        %base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);
%
%        % Get the stim idxs for the sustained period
%        get_sus_stim_idx = @(times) find(times >= 0 + Multi_func.sus_ped(1)/1000 & times < 0 + Multi_func.sus_ped(2)/1000);
%        stim_sus_idx = cellfun(get_sus_stim_idx, tim_stamp, 'UniformOutput', false);
%
%        % Calculate the firing rate for the sustained period
%        base_fr = cellfun(calc_fr, all_fr, base_sus_idx, 'UniformOutput', false);
%        stim_fr = cellfun(calc_fr, all_fr, stim_sus_idx, 'UniformOutput', false);
%
%        pop_sus_frs = cellfun(@(b_fr, s_fr) s_fr - b_fr, base_fr, stim_fr, 'UniformOutput', false);
%        pop_sus_frs = cell2mat(pop_sus_frs);
%        
%        % Initial p-value testing for each period
%        disp([f_region ' ' f_stim ]);
%        % Testing the sign stuff here
%        trans_fr_p = signtest(pop_trans_frs)
%        sus_fr_p = signtest(pop_sus_frs)




        % Calculate firing rate with the sum of spike events per the time of the area
        
        % Grab the baseline and period timestamp idxs
        %a Old way of calculating the transient and sustained period values
        %tim_stamp = popul_data.trace_timestamps(:, nr_idxs);
        %[trans_r, trans_c] = find(tim_stamp > Multi_func.trans_ped(1)/1000 & tim_stamp < Multi_func.trans_ped(2)/1000);
        %[sus_r, sus_c] = find(tim_stamp > Multi_func.sus_ped(1)/1000 & tim_stamp < Multi_func.sus_ped(2)/1000);
        %[base_r, base_c] = find(tim_stamp > Multi_func.base_ped(1)/1000 & tim_stamp < Multi_func.base_ped(2)/1000);
        %
        %base_r = unique(base_r);
        %base_c = unique(base_c);
        %trans_r = unique(trans_r);
        %trans_c = unique(trans_c);
        %sus_r = unique(sus_r);
        %sus_c = unique(sus_c);

        % Store the number of neurons
        num_neurons = length(nr_idxs);

        % Grab the spike raster
        all_sp = popul_data.all_trial_spike_rasters(nr_idxs);

        % Store the neuron timestamps into a cell array for compatibility with cellfun()
        tim_stamp = popul_data.trace_timestamps(:, nr_idxs);
        tim_stamp = mat2cell(tim_stamp, size(tim_stamp, 1), ones(1, size(tim_stamp, 2) ) );

        % Calculate the spike rate for the given period
        calc_fr = @(nr_sp, idxs, times) sum(nr_sp(idxs, :), 1)./( range(times(idxs), 1 ) ); 

        % Get the baseline idxs for transient period reflected over the onset into the baseline
        %get_trans_base_idx = @(times) find(times >= 0 - Multi_func.trans_ped(2)/1000 & times < 0 - Multi_func.trans_ped(1)/1000);
        %base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);
        
        % Get the baseline as the 500 ms right before onset
        get_trans_base_idx = @(times) find(times > Multi_func.base_ped(1)/1000 & times < Multi_func.base_ped(2)/1000);
        base_trans_idx = cellfun(get_trans_base_idx, tim_stamp, 'UniformOutput', false);

        % Get the stim idxs for the transient period
        get_trans_stim_idx = @(times) find(times >= 0 + Multi_func.trans_ped(1)/1000 & times < 0 + Multi_func.trans_ped(2)/1000);
        stim_trans_idx = cellfun(get_trans_stim_idx, tim_stamp, 'UniformOutput', false);

        % Calculate the firing rate for the transient period
        base_fr = cellfun(calc_fr, all_sp, base_trans_idx, tim_stamp, 'UniformOutput', false);
        stim_fr = cellfun(calc_fr, all_sp, stim_trans_idx, tim_stamp, 'UniformOutput', false);
        
        pop_trans_frs = cellfun(@(b_fr, s_fr) mean(s_fr - b_fr), base_fr, stim_fr, 'UniformOutput', false);
        pop_trans_frs = cell2mat(pop_trans_frs);
 

        % Get the baseline idxs for sustained period reflected over the onset into the baseline
        %get_sus_base_idx = @(times) find(times >= 0 - Multi_func.sus_ped(2)/1000 & times < 0 - Multi_func.sus_ped(1)/1000);
        %base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);
        
        % Check the 500 ms just before onset
        get_sus_base_idx = @(times) find(times > Multi_func.base_ped(1)/1000 & times < Multi_func.base_ped(2)/1000 );
        base_sus_idx = cellfun(get_sus_base_idx, tim_stamp, 'UniformOutput', false);

        % Get the stim idxs for the sustained period
        get_sus_stim_idx = @(times) find(times >= 0 + Multi_func.sus_ped(1)/1000 & times < 0 + Multi_func.sus_ped(2)/1000);
        stim_sus_idx = cellfun(get_sus_stim_idx, tim_stamp, 'UniformOutput', false);

        % Calculate the firing rate for the sustained period
        base_fr = cellfun(calc_fr, all_sp, base_sus_idx, tim_stamp, 'UniformOutput', false);
        stim_fr = cellfun(calc_fr, all_sp, stim_sus_idx, tim_stamp, 'UniformOutput', false);

        pop_sus_frs = cellfun(@(b_fr, s_fr) mean(s_fr - b_fr), base_fr, stim_fr, 'UniformOutput', false);
        pop_sus_frs = cell2mat(pop_sus_frs);

        % Testing the sign stuff here
        trans_fr_p = signtest(pop_trans_frs);
        sus_fr_p = signtest(pop_sus_frs);
        
        disp([f_region ' ' f_stim ]);
        disp(['Num neurons ' num2str(length(pop_trans_frs) )]);
        disp(['Trans fr p=' num2str(trans_fr_p)]);
        disp(['Sus fr p=' num2str(sus_fr_p)]);
        fprintf('\n\n');
        

        %DEBUG
        if strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_40') == 1
        %    figure;
        %    plot(ones(size(pop_trans_frs)), pop_trans_frs', '.');
        %    title('trans');
        %    
        %    figure;
        %    plot(ones(size(pop_sus_frs)), pop_sus_frs', '.');
        %    title('sus');
            
        %    figure;
        %    nrs_base = cell2mat(cellfun(@(x) mean(x, 'all'), base_fr, 'UniformOutput', false ) );
        %    nrs_stim = cell2mat(cellfun(@(x) mean(x, 'all'), stim_fr, 'UniformOutput', false ) );
        %    plot([1, 2], [nrs_base(:), nrs_stim(:)]);
        %    error('Done plot');;
        end


        % Plot violins
        labels = [];
        data = [];
        
        labels = [repmat({'Trans'}, num_neurons, 1); repmat({'Sus'}, num_neurons, 1)];
            
        data = [pop_trans_frs, pop_sus_frs];
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'MedianColor', [1 0 0], 'GroupOrder', {'Trans', 'Sus'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.trans_color};
        violins(2).ViolinColor = {Multi_func.sus_color};

        hold on;
        yline(0);

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Firing Rate Change (Hz)');
        
        % Change axis cosmetics
        a = gca; y = a.YAxis;
        Multi_func.set_default_axis(gca);
        
        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1
            y.Limits = [-10 40];
            Multi_func.set_spacing_axis(y, 4, 1);
        elseif strcmp(f_region, 'r_M1') == 1
            y.Limits = [-6.5 20];
            Multi_func.set_spacing_axis(y, 5, 1);
        elseif strcmp(f_region, 'r_V1') == 1
            y.Limits = [-19 70];
            Multi_func.set_spacing_axis(y, 20, 1);
        end
        
        %fr_stat_data.(f_region).(f_stim) = struct();
        fr_stat_data.(f_region).(f_stim).trans_fr.(nr_pop) = pop_trans_frs;
        fr_stat_data.(f_region).(f_stim).sus_fr.(nr_pop)= pop_sus_frs;

        % Just stop for V1 and 140 Hz stim
        %if strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_140') == 1
        %    error('pause');
        %end

    end
    
    sgtitle([f_region(3:end) ' Firing Rate Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin_' nr_pop '.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin_' nr_pop '.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.eps']);
end

%%
% Show violin plot on "stimulation" period for Firing Rate 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 5, 3]);
    for f_stim=stims'
        nexttile;
        % Plot violins
        labels = [];
        data = [];
        f_stim = f_stim{1};
        num_neurons = length(data_bystim.(f_stim).neuron_stim_FR);
        labels = [labels; repmat({'Stim'}, num_neurons, 1)];
        data = [data; data_bystim.(f_stim).neuron_stim_FR'];
        %TODO need to fix this, it is actually wrong because the +/- will change the exponential to either the fraction or what not, so really need to multiply the sign afterwards
        data = sign(data).*log10(1 + abs(data)/(10^1));
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'MedianColor', [1 0 0], 'GroupOrder', {'Stim'}, ViolinOpts);
        
        violins(1).ViolinColor = {Multi_func.stim_color};

        hold on;

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        y = gca; y = y.YAxis;
        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1
            y.Limits = [-10 40];
            Multi_func.set_spacing_axis(y, 4, 1);
        elseif strcmp(f_region, 'r_M1') == 1
            y.Limits = [-.1 0.4];
            Multi_func.set_spacing_axis(y, 0.2, 1);
        elseif strcmp(f_region, 'r_V1') == 1
            y.Limits = [-.1 0.8];
            Multi_func.set_spacing_axis(y, 0.4, 1);
        end

        
        %TODO resize the xaxis here before adding the 10 exnential value
        %TODO need to also show 0 tick, do this by setting proper limits and the spacing axis
        y_sign = arrayfun(@num2str, 10*sign(yticks), 'UniformOutput', 0);
        y_vals = arrayfun(@num2str, abs(yticks), 'UniformOutput', 0);

        ylab = strcat(y_sign, '^{', y_vals, '}');
        zero_idx = find(sign(yticks) == 0);
        if ~isempty(zero_idx)
            ylab{zero_idx} = '0';
        end
        yticklabels(ylab);
        
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Firing Rate Change (Hz)');

        fr_stat_data.(f_region).(f_stim).stim_fr = data_bystim.(f_stim).neuron_stim_FR;
    end
    
    sgtitle([f_region(3:end) ' Firing Rate Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_FR_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_FR_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.eps']);
end


%% All Vm showing all DBS pulses
% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

% Flag for removing non-modulated neurons from analysis
remove_nonmod_nrs = 1;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 9.5, 5]);
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
            
        % Filter out the neurons that were non-modulated at all
        if remove_nonmod_nrs == 1
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_idxs = nr_idxs(~ismember(nr_idxs, non_mod_nr));
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        norm_vms = popul_data.neuron_RawVm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_RawVm, 1);
        nexttile;

        % Plot the subthreshold Vm
        fill_h = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        hold on;
        
        % Plot the DBS stimulation time pulses
        stim_time = nanmedian(popul_data.stim_timestamps, 2);
        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;


        % Plot the timescale bar
        posx = -.100;
        posy = -0.5;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.1, '50ms');
        hold on;

        % Plot the Vm scale
        poxs = .2;
        posy = 0;
        Vm_scale = 0.5;
        plot([posx, posx], [posy, posy + Vm_scale], 'k', 'LineWidth', 2);
        text(posx - .01, posy, [num2str(Vm_scale) ' Norm Vm'], 'Rotation', 90);

        % Plot the zoom-in outlines for the onset and some points in the stimulation period
        min_val = min(cur_Vm - sem_Vm);
        max_val = max(cur_Vm + sem_Vm);
        plot([Multi_func.onset_ped(1), Multi_func.onset_ped(2), Multi_func.onset_ped(2), ...
            Multi_func.onset_ped(1), Multi_func.onset_ped(1) ]./1000, ...
            [min_val, min_val, max_val, max_val, min_val], ...
            'LineWidth', 0.5, 'Color', Multi_func.trans_color);
        hold on;

        min_val = min(cur_Vm - sem_Vm);
        max_val = max(cur_Vm + sem_Vm);
        plot([Multi_func.mid_stim_ped(1), Multi_func.mid_stim_ped(2), Multi_func.mid_stim_ped(2), ...
            Multi_func.mid_stim_ped(1), Multi_func.mid_stim_ped(1) ]./1000, ...
            [min_val, min_val, max_val, max_val, min_val], ...
            'LineWidth', 0.5, 'Color', Multi_func.sus_color);
        hold on;

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]); % Includes surrounding component
        %xlim([-0.80 2.05]); % Same time resolution as power spectra
        a = gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        set(gca, 'color', 'none')
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Average subthreshold Vm Showing all pulses ' nr_pop], 'Interpreter', 'none');
   
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_Vm.pdf']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_Vm.fig']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Vm.eps'], 'epsc');
end

%% Continuous firing rate showing all DBS pulses

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

% Flag for removing non-modulated neurons from analysis
remove_nonmod_nrs = 1;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 9.5, 5.0]);
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

        % Filter out the neurons that were non-modulated at all
        if remove_nonmod_nrs == 1
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_idxs = nr_idxs(~ismember(nr_idxs, non_mod_nr));
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        
        Fs = mean(popul_data.framerate);

        % -- Use a proper window for estimating spike rate
        %cur_srate = mean(popul_data.neuron_srate_3(:, nr_idxs), 2, 'omitnan');
        %std_srate = std(popul_data.neuron_srate_3(:, nr_idxs), 0, 2, 'omitnan');
        %num_neurons = size(popul_data.neuron_srate_3(:, nr_idxs), 2);
        %sem_srate = std_srate./sqrt(num_neurons);
        % -- End of proper window rate

        % -- Use instantaneous rate
        cur_srate = mean(popul_data.neuron_spikecounts_raster(:, nr_idxs)*Fs, 2, 'omitnan');
        std_srate = std(popul_data.neuron_spikecounts_raster(:, nr_idxs)*Fs, 0, 2, 'omitnan');
        num_neurons = size(popul_data.neuron_spikecounts_raster(:, nr_idxs), 2);
        sem_srate = std_srate./sqrt(num_neurons);
        % -- End of instantaneous rate

        %num_points = size(popul_data.neuron_srate, 1);
        
        nexttile;

        % Plot the subthreshold srate
        fill_h = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;
        
        % Plot the DBS stimulation time pulses
        stim_time = nanmean(popul_data.stim_timestamps, 2);
        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the timescale bar
        posx = -.100;
        posy = -8;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
        text(posx, posy - 0.5, '50ms', 'FontSize', 7);
        hold on;

        % Plot the srate scale
        poxs = -0.5;
        posy = 400;
        srate_scale = 100;
        plot([posx, posx], [posy, posy + srate_scale], 'k', 'LineWidth', 2);
        text(posx - .02, posy, [num2str(srate_scale) ' FR (Hz)'], 'Rotation', 90);

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]);
        a = gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        set(gca, 'color', 'none')
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Population Firing Rate Showing All Pulses ' nr_pop], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_FR.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.eps'], 'epsc');
end


%% Stats for Vm and Firing Rate of Transient and Sustained Period

nr_pop = 'all_mod';
%nr_pop = 'non_mod';
%nr_pop = 'vm_act';
%nr_pop = 'vm_sup';
%nr_pop = 'vm_non_mod';

% signtest for individual: trans, sus, stim
% signrank for paired, non-independent: 140 trans vs. 140 sus
% ranksum paired, independent distributions: 140 trans vs 40 trans
stats_log = [figure_path 'Average' f 'Population_Trans_and_Sustained_Vm_FR_stats_' nr_pop]
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

% Gonna use table instead
stats_t = table();

% Loop through each region and stimulation to perform statistical tests
for f_region = fieldnames(sub_vm_stat_data)'
    f_region = f_region{1};
        
    for f_stim = fieldnames(sub_vm_stat_data.(f_region))'
        f_stim = f_stim{1};
        
        % Log diary for conditions
        diary on;
        disp([f_region ' ' f_stim]);
        fprintf('\n');
        diary off;

        % Perform individual signtests on subvm
        for f_ped = fieldnames(sub_vm_stat_data.(f_region).(f_stim))'
            f_ped = f_ped{1};
            
            if contains(f_ped, 'stats')
                continue;
            end

            diary on;
            disp(['Period signtest from baseline ' f_ped]);
            [p, h, stats] = signtest(sub_vm_stat_data.(f_region).(f_stim).(f_ped).(nr_pop))
            diary off;
            
            % Add to stats table
            stats_t([f_region f_stim], ['vm_' f_ped ' p-val']) = {p};
            stats_t([f_region f_stim], ['vm_' f_ped ' statistic']) = {stats.sign};
            
            clear p, h, stats;
        end

        % Perform the trans and sus comparison (Vm)
        diary on;
        disp(['Trans Vm vs Sustained Vm Wilcoxon signed rank'])
        [p, h, stats] = signrank(sub_vm_stat_data.(f_region).(f_stim).trans_vm.(nr_pop), ...
            sub_vm_stat_data.(f_region).(f_stim).sus_vm.(nr_pop))
        diary off;

        % Add to stats table
        stats_t([f_region f_stim], 'vm trans vs sus p-value') = {p};
        stats_t([f_region f_stim], 'vm trans vs sus statistic') = {stats.signedrank};

        clear p, h, stats;

        % Perform individual signtests on firing rate
        for f_ped = fieldnames(fr_stat_data.(f_region).(f_stim))'
            f_ped = f_ped{1};

            if contains(f_ped, 'stats')
                continue;
            end

            diary on;
            disp(['Period signtest from baseline ' f_ped]);
            [p, h, stats] = signtest(fr_stat_data.(f_region).(f_stim).(f_ped).(nr_pop))
            diary off;
            
            % Add to stats table
            stats_t([f_region f_stim], ['fr_' f_ped ' p-value']) = {p};
            stats_t([f_region f_stim], ['fr_' f_ped ' statistic']) = {stats.sign};

            clear p, h, stats;
        end

        % Perform transient vs sustained firing rate comparison
        diary on;
        disp(['Trans FR vs Sustained FR Wilcoxon signed rank'])
        % Perform the trans and sus comparison (Firign Rate)
        [p, h, stats] = signrank(fr_stat_data.(f_region).(f_stim).trans_fr.(nr_pop), fr_stat_data.(f_region).(f_stim).sus_fr.(nr_pop))
        diary off;

        % Add to stats table
        stats_t([f_region f_stim], 'fr trans vs sus p-value') = {p};
        stats_t([f_region f_stim], 'fr trans vs sus statistic') = {stats.signedrank};
    end
end

disp(stats_t);
writetable(stats_t, [figure_path 'Average' f 'Population_pvals_trans_sus_Vm_FR_stats_' nr_pop '.csv'], 'WriteRowNames', true);


%% Make Violin plots on Subthreshold Vm for stimulation period

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 5, 3]);
    for f_stim=stims'
        nexttile;
        % Plot violins
        labels = [];
        data = [];
        f_stim = f_stim{1};
        num_neurons = length(data_bystim.(f_stim).neuron_stim_Vm);
        labels = [labels; repmat({'Stim'}, num_neurons, 1)];
        data = [data; data_bystim.(f_stim).neuron_stim_Vm'];
 
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Stim'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.stim_color};

        hold on;

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Vm Change');
        
        % Change axis cosmetics
        a = gca; y = a.YAxis;
        Multi_func.set_default_axis(gca);
        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1
            y.Limits = [-10 40];
            Multi_func.set_spacing_axis(y, 4, 1);
        elseif strcmp(f_region, 'r_M1') == 1
            y.Limits = [-8 20];
            Multi_func.set_spacing_axis(y, 10, 1);
        elseif strcmp(f_region, 'r_V1') == 1
            y.Limits = [-5 45];
            Multi_func.set_spacing_axis(y, 20, 1);
        end

        sub_vm_stat_data.(f_region).(f_stim).stim_vm = data_bystim.(f_stim).neuron_stim_Vm;
    end
    
    sgtitle([f_region(3:end) ' Sub Vm Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_Vm_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_Vm_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.eps'], 'epsc');
end
