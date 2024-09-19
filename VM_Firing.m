clear all;
close all;
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
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 17.56, 3]);
    for f_stim=stims'
        f_stim = f_stim{1};
        cur_win_srate = 50;
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(f_stim).neuron_srate_50, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim).neuron_srate_50, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim).neuron_srate_50, 2);
        %num_points = size(data_bystim.(f_stim).neurons_srate_50, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile;
        fill_h = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;
        
        % Some axis stuff
        a = gca; y = a.YAxis;
        Multi_func.set_default_axis(gca);
 
        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1
            Multi_func.plot_dbs_bar([0, 1], 7, [f_stim(3:end) 'Hz DBS']);
            y.Limits = [-2 8];
            Multi_func.set_spacing_axis(y, 2, 1);
        elseif strcmp(f_region, 'r_M1') == 1
            Multi_func.plot_dbs_bar([0, 1], 4, [f_stim(3:end) 'Hz DBS']);
            y.Limits = [-2 5];
            Multi_func.set_spacing_axis(y, 2, 1);
        elseif strcmp(f_region, 'r_V1') == 1
            Multi_func.plot_dbs_bar([0, 1], 18, [f_stim(3:end) 'Hz DBS']);
            y.Limits = [-5 20];
            Multi_func.set_spacing_axis(y, 5, 1);
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
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate' num2str(cur_win_srate) '.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate' num2str(cur_win_srate) '.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.eps']);
end

%% Compact population Subthreshold Vm

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    fontsize(gcf, 7, "points")
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
            end
        catch ME
            disp(ME.message);
        end

        timeline = nanmean(popul_data.trace_timestamps, 2);
        norm_vms = popul_data.neuron_Vm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;
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
        Multi_func.plot_dbs_bar([0, 1], max(cur_Vm)+ 0.3, [f_stim(3:end) 'Hz DBS']);

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
    sgtitle([f_region(3:end) ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Average_sub_thres.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Average_sub_thres.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.eps']);
end


%% Violin plots sub vm transient and sustained periods
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
        nexttile;

        % Normalize the Vm
        norm_vms = data_bystim.(f_stim).neuron_Vm./data_bystim.(f_stim).neuron_spike_amp;
        
        tim_stamp = data_bystim.(f_stim).trace_timestamps;
        [trans_r, trans_c] = find(tim_stamp > Multi_func.trans_ped(1)/1000 & tim_stamp < Multi_func.trans_ped(2)/1000);
        [sus_r, sus_c] = find(tim_stamp > Multi_func.sus_ped(1)/1000 & tim_stamp < Multi_func.sus_ped(2)/1000);
        [base_r, base_c] = find(tim_stamp > Multi_func.base_ped(1)/1000 & tim_stamp < Multi_func.base_ped(2)/1000);
        
        base_r = unique(base_r);
        base_c = unique(base_c);
        trans_r = unique(trans_r);
        trans_c = unique(trans_c);
        sus_r = unique(sus_r);
        sus_c = unique(sus_c);

        % Find the base values and then reshape the array to maintain length
        neuron_base_vm = nanmean(norm_vms(base_r, base_c), 1);
        neuron_trans_vm = nanmean(norm_vms(trans_r, trans_c), 1) - neuron_base_vm;
        neuron_sus_vm = nanmean(norm_vms(sus_r, sus_c), 1) - neuron_base_vm;

        num_neurons = size(neuron_trans_vm, 2);
 
        % Plot violins
        labels = [];
        data = [];
        labels = [labels; repmat({'Trans'}, num_neurons, 1); repmat({'Sus'}, num_neurons, 1)];
        data = [data; neuron_trans_vm', neuron_sus_vm'];
 
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
        y.Limits = [-1, 2];
        Multi_func.set_spacing_axis(y, 1, 2);

        % Old way of setting axis limits
        % Determine parameters based on what is being plotted
        %if strcmp(f_region, 'r_combine') == 1
        %    y.Limits = [-2 8];
        %    Multi_func.set_spacing_axis(y, 4, 1);
        %elseif strcmp(f_region, 'r_M1') == 1
        %    y.Limits = [-1 2];
        %    %Multi_func.set_spacing_axis(y, 10, 1);
        %elseif strcmp(f_region, 'r_V1') == 1
        %    y.Limits = [-5 45];
        %    Multi_func.set_spacing_axis(y, 20, 1);
        %end

        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Normalized Vm Change');

        sub_vm_stat_data.(f_region).(f_stim).trans_vm = neuron_trans_vm;
        sub_vm_stat_data.(f_region).(f_stim).sus_vm = neuron_sus_vm;
    end
    
    sgtitle([f_region(3:end) ' Sub Vm Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.eps'], 'epsc');
end

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

%% Violin plots firing rate transient and sustained periods

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
        nexttile;
        % Plot violins
        labels = [];
        data = [];
        f_stim = f_stim{1};
        num_neurons = length(data_bystim.(f_stim).neuron_trans_FR);
        labels = [labels; repmat({'Trans'}, num_neurons, 1); repmat({'Sus'}, num_neurons, 1)];
        data = [data; data_bystim.(f_stim).neuron_trans_FR', data_bystim.(f_stim).neuron_sus_FR'];
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
            y.Limits = [-2 15];
            Multi_func.set_spacing_axis(y, 5, 1);
        elseif strcmp(f_region, 'r_V1') == 1
            y.Limits = [-5 70];
            Multi_func.set_spacing_axis(y, 20, 1);
        end

        %fr_stat_data.(f_region).(f_stim) = struct();
        fr_stat_data.(f_region).(f_stim).trans_fr = data_bystim.(f_stim).neuron_trans_FR;
        fr_stat_data.(f_region).(f_stim).sus_fr = data_bystim.(f_stim).neuron_sus_FR;
    end
    
    sgtitle([f_region(3:end) ' Firing Rate Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.eps']);
end

% Show violin plot on stimulation period for Firing Rate 
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


%% Subthreshold Vm showing all DBS pulses

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
%nr_pop = 'etrain';
nr_pop = 'non';

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 9, 5.23]);
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
        norm_vms = popul_data.neuron_Vm(:, nr_idxs)./popul_data.neuron_spike_amp(nr_idxs);
        cur_Vm = mean(norm_vms, 2, 'omitnan');
        std_Vm = std(norm_vms, 0, 2, 'omitnan');
        num_neurons = size(norm_vms, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
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

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]);
        a = gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        set(gca, 'color', 'none')
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Average subthreshold Vm Showing all pulses ' nr_pop], 'Interpreter', 'none');
   
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_' nr_pop '_Display_All_Pulse_Vm.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Vm.eps'], 'epsc');
end

%% Continuous firing rate showing all DBS pulses

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
%nr_pop = 'etrain';
nr_pop = 'non';

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

        timeline = nanmean(popul_data.trace_timestamps, 2);
        cur_srate = mean(popul_data.neuron_srate_3(:, nr_idxs), 2, 'omitnan');
        std_srate = std(popul_data.neuron_srate_3(:, nr_idxs), 0, 2, 'omitnan');
        num_neurons = size(popul_data.neuron_srate_3(:, nr_idxs), 2);
        sem_srate = std_srate./sqrt(num_neurons);
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
        posy = -3;
        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 0.5);
        text(posx, posy - 0.5, '50ms', 'FontSize', 7);
        hold on;

        % Plot the srate scale
        poxs = .2;
        posy = 5;
        srate_scale = 10;
        plot([posx, posx], [posy, posy + srate_scale], 'k', 'LineWidth', 0.5);
        text(posx - .01, posy, [num2str(srate_scale) ' FR (Hz)'], 'Rotation', 90, 'FontSize', 7);

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


%% Calculate statistics for Vm and Firing Rate of Transient and Sustained Period
% signtest for individual: trans, sus, stim
% signrank for paired, non-independent: 140 trans vs. 140 sus
% ranksum paired, independent distributions: 140 trans vs 40 trans
stats_log = [figure_path 'Average' f 'Population_Trans_and_Sustained_Vm_FR_stats']
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

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
            disp(['Period ' f_ped]);
            [p, h, stats] = signtest(sub_vm_stat_data.(f_region).(f_stim).(f_ped))
            diary off;

            clear p, h, stats;
        end

        % Perform individual signtests on firing rate
        for f_ped = fieldnames(fr_stat_data.(f_region).(f_stim))'
            f_ped = f_ped{1};

            if contains(f_ped, 'stats')
                continue;
            end

            diary on;
            disp(['Period ' f_ped]);
            [p, h, stats] = signtest(fr_stat_data.(f_region).(f_stim).(f_ped))
            diary off;

            clear p, h, stats;
        end

        % Perform the trans and sus comparison (Vm)
        diary on;
        disp(['Trans Vm vs Sustained Vm'])
        [p, h, stats] = signrank(sub_vm_stat_data.(f_region).(f_stim).trans_vm, sub_vm_stat_data.(f_region).(f_stim).sus_vm)
        diary off;

        clear p, h, stats;

        diary on;
        disp(['Trans FR vs Sustained FR'])
        % Perform the trans and sus comparison (Firign Rate)
        [p, h, stats] = signrank(fr_stat_data.(f_region).(f_stim).trans_fr, fr_stat_data.(f_region).(f_stim).sus_fr)
        diary off;
    end
end

%% Functin to calculate the power spectra
% Calculate cwt for input signal and 
function [wt, f] = get_power_spec(signal, samp_freq)
    freqLimits = [0 150];
    fb = cwtfilterbank(SignalLength=length(signal),...
                       SamplingFrequency=samp_freq,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(signal, FilterBank=fb);
end
