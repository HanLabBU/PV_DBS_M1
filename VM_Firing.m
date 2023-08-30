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

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

srate_win = 100;

% Time periods for comparison of firing rate and sub Vm
trans_ped = [0, 150];
sus_ped = [150, 1000];
stim_ped = [0, 1000];

%%% END Modification

%% Check that the server path exists
%if ~isfolder(local_root_path)
%    disp('Server rootpath does not exist!!!');
%    return;
%end
%
%ses = dir([pv_data_path '*.mat']);
%
%all_matfiles = {ses.name};
%
%% Select matfiles by brain region
%[region_matfiles] = Multi_func.find_region(all_matfiles);
%region_data = struct();
%all_Fs = [];
%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%    %% Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    %% Select matfiles by stim condition for given region
%    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);
%
%    %% Loop through each field of the struct and concatenate everything together
%    % Store trace aspect data by each individual stimulation condition
%    data_bystim = struct();
%    % Store all of the calculated sampling frequencies
%
%    % Loop through each stimulation condition
%    for f_stim = fieldnames(matfile_stim)'
%        f_stim = f_stim{1};
%        matfiles = matfile_stim.(f_stim).names;    
%    
%        % Initialize field subthreshold array
%        data_bystim.(f_stim) = struct();
%        data_bystim.(f_stim).neuron_Vm = [];
%        data_bystim.(f_stim).neuron_srate = [];
%        data_bystim.(f_stim).stim_timestamps = [];
%        data_bystim.(f_stim).trace_timestamps = [];
%
%        % Store transient information    
%        data_bystim.(f_stim).neuron_trans_Vm = [];      
%        data_bystim.(f_stim).neuron_sus_Vm = [];        
%        data_bystim.(f_stim).neuron_trans_FR = [];     
%        data_bystim.(f_stim).neuron_sus_FR = [];
%
%        % Loop through each matfile of the current stimulation condition
%        for matfile = matfiles
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile{1}]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%            cur_fov_Fs = [];
%            cur_fov_subVm = [];
%            cur_fov_srate = [];
%            cur_fov_raster = [];
%            cur_fov_stim_time = [];
%            cur_fov_trace_time = [];
%         
%            cur_fov_trans_Vm = [];      
%            cur_fov_sus_Vm = [];        
%            cur_fov_trans_FR = [];     
%            cur_fov_sus_FR = [];
%            cur_fov_base_Vm = [];
%            cur_fov_base_FR = [];
%
%            % Loop through each ROI
%            for roi_idx=1:size(trial_data.detrend_traces, 2)
%                %Determine whether this roi is to be ignored for this particular trial
%                ri = strsplit(matfile{1}, '_');
%                try
%                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
%                catch
%                    trial_ignr_list = [];
%                end
%
%                % Remove trials from trial idx list
%                trial_idxs = setdiff(trial_idxs, trial_ignr_list);
%
%                % Skip ROI if there are at most 2 trials
%                if length(trial_idxs) <= 2
%                    continue;
%                end
%
%                % Loop through each trial                
%                for tr_idx=trial_idxs        
%                    trial_data = data.align.trial{tr_idx};
%                    raw_trial_data = data.raw.trial{tr_idx};
%
%                    % Store the camera framerate
%                    all_Fs(end+1) = trial_data.camera_framerate;
%                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;
%
%                    % Store all of the timestamp info
%                    stim_start = raw_trial_data.raw_stimulation_time(1);
%                    cur_stim_time = raw_trial_data.raw_stimulation_time - stim_start;
%                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
%                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, cur_stim_time);
%                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, cur_trace_time);
%                    
%                    % DEBUG
%                    if max(cur_fov_trace_time(:, end)) < 2
%                        stim_trace_sz = size(cur_fov_stim_time)
%                        trace_time_sz = size(cur_fov_trace_time)
%                        pause;
%                        matfile{1}
%                    end
%
%                    % Grab the subthreshold Vm
%                    % Chop the respective frames
%                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
%                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
%                    detrend_subVm = cur_trace_ws - baseline;
%                    % This normalization did not seem to change the plot much
%                    %detrend_subVm = detrend_subVm./mean(trial_data.spike_info375.spike_amplitude{1}, 'omitnan');
%                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');
%                    
%                    % Calculate the spike rate
%                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
%                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
%                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win, 1, 1);
%                    % Baseline subtract the mean
%                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
%                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
%
%                    cur_fov_srate = horzcat_pad(cur_fov_srate, cur_spikerate');            
%                    % Store the raster plot
%                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');
%
%                    % Grab transient and sustained period for spikes and Vm
%                    base_idx = find(cur_trace_time < cur_stim_time(1));
%                    trans_idx = find(cur_trace_time > 0 & cur_trace_time <= trans_ped(2)./1000);
%                    sus_idx = find(cur_trace_time > trans_ped(2)./1000 & cur_trace_time <= cur_stim_time(end));
%                    
%                    cur_fov_base_Vm = horzcat_pad(cur_fov_base_Vm, mean(detrend_subVm(base_idx)', 'omitnan'));
%                    cur_fov_trans_Vm = horzcat_pad(cur_fov_trans_Vm, mean(detrend_subVm(trans_idx)', 'omitnan'));
%                    cur_fov_sus_Vm = horzcat_pad(cur_fov_sus_Vm, mean(detrend_subVm(sus_idx)', 'omitnan'));
%
%                    cur_fov_base_FR = horzcat_pad(cur_fov_base_FR, sum(cur_raster(base_idx))');
%                    cur_fov_trans_FR = horzcat_pad(cur_fov_trans_FR, sum(cur_raster(trans_idx))./(diff(trans_ped)./1000)');
%                    cur_fov_sus_FR = horzcat_pad(cur_fov_sus_FR, sum(cur_raster(sus_idx))./(diff(sus_ped)./1000)');
%
%                end % End looping through each neuron
%            end
%            
%            % Skip rest of the calculations if the subthreshold Vm is nan
%            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
%                continue;
%            end
%
%            % Average for each neuron and save the subthreshold Vm
%            temp = data_bystim.(f_stim).neuron_Vm;
%            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
%            % Store average spike rate for each neuron
%            temp = data_bystim.(f_stim).neuron_srate;
%            data_bystim.(f_stim).neuron_srate = horzcat_pad(temp, nanmean(cur_fov_srate, 2));
%            % Store the timestamp data
%            temp = data_bystim.(f_stim).stim_timestamps;
%            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
%            temp = data_bystim.(f_stim).trace_timestamps;
%            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
%            
%            temp = data_bystim.(f_stim).neuron_trans_Vm;
%            data_bystim.(f_stim).neuron_trans_Vm = horzcat_pad(temp, mean(cur_fov_trans_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
%            temp = data_bystim.(f_stim).neuron_sus_Vm;
%            data_bystim.(f_stim).neuron_sus_Vm = horzcat_pad(temp, mean(cur_fov_sus_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
%            temp = data_bystim.(f_stim).neuron_trans_FR;
%            data_bystim.(f_stim).neuron_trans_FR = horzcat_pad(temp, mean(cur_fov_trans_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));
%            temp = data_bystim.(f_stim).neuron_sus_FR;
%            data_bystim.(f_stim).neuron_sus_FR = horzcat_pad(temp, mean(cur_fov_sus_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));
%
%        end % End looping through FOVs of a condition
%    end
%
%    % Save the VM to the specific region
%    region_data.(f_region).data_bystim = data_bystim;
%end
%---- End of data collection ---

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
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
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

% Set figure off currently
set(0,'DefaultFigureVisible','on');

%%Compact full collective spike rate over time
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 10.68, 3.0886]);
    for f_stim=stims'
        f_stim = f_stim{1};
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(f_stim).neuron_srate_large, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim).neuron_srate_large, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim).neuron_srate_large, 2);
        %num_points = size(data_bystim.(f_stim).neuron_srate_large, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;

        % Determine parameters based on what is being plotted
        if strcmp(f_region, 'r_combine') == 1
            Multi_func.plot_dbs_bar([0, 1], 7, [f_stim(3:end) 'Hz DBS']);
            y.Limits = [-2 8];
            Multi_func.set_spacing_axis(y, 6, 1);
        elseif strcmp(f_region, 'r_M1') == 1
            Multi_func.plot_dbs_bar([0, 1], 2, [f_stim(3:end) 'Hz DBS']);
            y.Limits = [-2 4];
            Multi_func.set_spacing_axis(y, 4, 1);
        end

        xlabel('Time from stim onset (s)');
        
        % Plot timeline depending on which chopping was done
        if exclude_200ms
            xlim([-0.85 2.05]);
        else
            xlim([-1 2.05]);
        end
        ylabel('Firing Rate (Hz)');
        a = gca; y = a.YAxis;
        
        Multi_func.set_default_axis(gca);
        %title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average Spike rate'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.eps']);
end

%%Compact population Subthreshold Vm
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    fontsize(gcf, 7, "points")
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 10.68, 3.0886]);
    for f_stim=stims'
        f_stim = f_stim{1};
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
        cur_Vm = mean(data_bystim.(f_stim).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5], 'linewidth', 0.2);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_Vm, 'k', 'LineWidth', 0.3);
        hold on;
        % Plot DBS pulse bar
        Multi_func.plot_dbs_bar([0, 1], 11, [f_stim(3:end) 'Hz DBS']);
        xlabel('Time from stim onset (s)');
        if exclude_200ms
            xlim([-0.85 2.05]);
        else
            xlim([-1 2.05]);
        end
        ylabel('Vm (A.U.)');
        Multi_func.set_default_axis(gca);
        x.Limits = [-4 12];

        % Determine limits based on what is plotted
        if strcmp(f_region, 'r_M1') == 1
            y.Limits = [-2 4];
        elseif strcmp(f_region, 'r_combine') == 1
            y.Limits = [-2 8];
        end

        Multi_func.set_spacing_axis(x, 5, 1);
        %title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.eps']);
end

% Make Violin plots on Subthreshold Vm for transient and sustained
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Struct to identify each group of data
    sub_vm_stat_data = struct();

    figure('visible', 'off', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 5, 3]);
    for f_stim=stims'
        nexttile;
        % Plot violins
        labels = [];
        data = [];
        f_stim = f_stim{1};
        num_neurons = length(data_bystim.(f_stim).neuron_trans_Vm);
        labels = [labels; repmat({'Trans'}, num_neurons, 1); repmat({'Sus'}, num_neurons, 1)];
        data = [data; data_bystim.(f_stim).neuron_trans_Vm', data_bystim.(f_stim).neuron_sus_Vm'];
 
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Trans', 'Sus'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.trans_color};
        violins(2).ViolinColor = {Multi_func.sus_color};

        hold on;

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Vm Change');
        ylim([-10 40]);

        %sub_vm_stat_data.(f_stim) = struct();
        sub_vm_stat_data.(f_stim).trans_vm = data_bystim.(f_stim).neuron_trans_Vm;
        sub_vm_stat_data.(f_stim).sus_vm = data_bystim.(f_stim).neuron_sus_Vm;
    end
    
    sgtitle([f_region ' Sub Vm Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.eps'], 'epsc');
end

% transient compared to baseline statistical test 
[p, h, stats] = signtest(sub_vm_stat_data.f_40.trans_vm)
SubVm_40_trans_stats = struct();
SubVm_40_trans_stats.p = p
SubVm_40_trans_stats.h = h;
SubVm_40_trans_stats.stats = stats;

[p, h, stats] = signtest(sub_vm_stat_data.f_140.trans_vm)
SubVm_140_trans_stats = struct();
SubVm_140_trans_stats.p = p
SubVm_140_trans_stats.h = h;
SubVm_140_trans_stats.stats = stats;

% transient compared to baseline statistical test 
[p, h, stats] = signtest(sub_vm_stat_data.f_40.sus_vm)
SubVm_40_sus_stats = struct();
SubVm_40_sus_stats.p = p
SubVm_40_sus_stats.h = h;
SubVm_40_sus_stats.stats = stats;

[p, h, stats] = signtest(sub_vm_stat_data.f_140.sus_vm)
SubVm_140_sus_stats = struct();
SubVm_140_sus_stats.p = p
SubVm_140_sus_stats.h = h;
SubVm_140_sus_stats.stats = stats;

% Trans and Sus SubVm statistical test
[p, h, stats] = signrank(sub_vm_stat_data.f_40.trans_vm, sub_vm_stat_data.f_40.sus_vm)
signrank_SubVm_40_trans_sus_stats = struct();
signrank_SubVm_40_trans_sus_stats.p = p
signrank_SubVm_40_trans_sus_stats.h = h;
signrank_SubVm_40_trans_sus_stats.stats = stats;

[p, h, stats] = signrank(sub_vm_stat_data.f_140.trans_vm, sub_vm_stat_data.f_140.sus_vm)
signrank_Subvm_140_trans_sus_stats = struct();
signrank_Subvm_140_trans_sus_stats.p = p
signrank_Subvm_140_trans_sus_stats.h = h;
signrank_Subvm_140_trans_sus_stats.stats = stats;

% Trans to Trans Subvm statistical test
[p, h, stats] = ranksum(sub_vm_stat_data.f_40.trans_vm, sub_vm_stat_data.f_140.trans_vm)
ranksum_Subvm_trans_trans_stats = struct();
ranksum_Subvm_trans_trans_stats.p = p
ranksum_Subvm_trans_trans_stats.h = h;
ranksum_Subvm_trans_trans_stats.stats = stats;

% Sus to sus Subvm statistical test
[p, h, stats] = ranksum(sub_vm_stat_data.f_40.sus_vm, sub_vm_stat_data.f_140.sus_vm)
ranksum_Subvm_sus_sus_stats = struct();
ranksum_Subvm_sus_sus_stats.p = p
ranksum_Subvm_sus_sus_stats.h = h;
ranksum_Subvm_sus_sus_stats.stats = stats;

% Make Violin plots on Subthreshold Vm for stimulation period
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Struct to identify each group of data
    sub_vm_stat_data = struct();

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
        
        if strcmp(f_region, 'r_M1') == 1
            ylim([-8 20]);
        elseif strcmp(f_region, 'r_combine') == 1
            ylim([-10 40]);
        end
        %sub_vm_stat_data.(f_stim) = struct();
        sub_vm_stat_data.(f_stim).stim_vm = data_bystim.(f_stim).neuron_stim_Vm;
        sub_vm_stat_data.(f_stim).stim_vm = data_bystim.(f_stim).neuron_stim_Vm;
    end
    
    sgtitle([f_region ' Sub Vm Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_Vm_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_Vm_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_Vm_violin.eps'], 'epsc');
end

[p, h, stats] = signtest(sub_vm_stat_data.f_140.stim_vm)
SubVm_140_stim_stats = struct();
SubVm_140_stim_stats.p = p
SubVm_140_stim_stats.h = h;
SubVm_140_stim_stats.stats = stats;

% Perform Statistical tests on Firing Rate 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Struct to identify each group of data
    fr_stat_data = struct();

    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 5, 3]);
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

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Firing Rate Change (Hz)');
        
        if strcmp(f_region, 'r_M1') == 1
            ylim([-5 20]);
        elseif strcmp(f_region, 'r_combine') == 1
            y.Limits = [-10 40];
        end

        %fr_stat_data.(f_stim) = struct();
        fr_stat_data.(f_stim).trans_fr = data_bystim.(f_stim).neuron_trans_FR;
        fr_stat_data.(f_stim).sus_fr = data_bystim.(f_stim).neuron_sus_FR;
    end
    
    sgtitle([f_region ' Firing Rate Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.eps']);
end

% transient compared to baseline statistical test 
[p, h, stats] = signtest(fr_stat_data.f_40.trans_fr)
FR_40_trans_stats = struct();
FR_40_trans_stats.p = p;
FR_40_trans_stats.h = h;
FR_40_trans_stats.stats = stats;

[p, h, stats] = signtest(fr_stat_data.f_140.trans_fr)
FR_140_trans_stats = struct();
FR_140_trans_stats.p = p;
FR_140_trans_stats.h = h;
FR_140_trans_stats.stats = stats;

% transient compared to baseline statistical test 
[p, h, stats] = signtest(fr_stat_data.f_40.sus_fr)
FR_40_sus_stats = struct();
FR_40_sus_stats.p = p;
FR_40_sus_stats.h = h;
FR_40_sus_stats.stats = stats;

[p, h, stats] = signtest(fr_stat_data.f_140.sus_fr)
FR_140_sus_stats = struct();
FR_140_sus_stats.p = p;
FR_140_sus_stats.h = h;
FR_140_sus_stats.stats = stats;

% Trans and Sus SubVm statistical test
[p, h, stats] = signrank(fr_stat_data.f_40.trans_fr, fr_stat_data.f_40.sus_fr);
signrank_FR_40_trans_sus_stats = struct();
signrank_FR_40_trans_sus_stats.p = p
signrank_FR_40_trans_sus_stats.h = h;
signrank_FR_40_trans_sus_stats.stats = stats;

[p, h, stats] = signrank(fr_stat_data.f_140.trans_fr, fr_stat_data.f_140.sus_fr);
signrank_FR_140_trans_sus_stats = struct();
signrank_FR_140_trans_sus_stats.p = p
signrank_FR_140_trans_sus_stats.h = h;
signrank_FR_140_trans_sus_stats.stats = stats;


% Trans to Trans Subfr statistical test
[p, h, stats] = ranksum(fr_stat_data.f_40.trans_fr, fr_stat_data.f_140.trans_fr);
ranksum_FR_trans_trans_stats = struct();
ranksum_FR_trans_trans_stats.p = p
ranksum_FR_trans_trans_stats.h = h;
ranksum_FR_trans_trans_stats.stats = stats;

% Sus to sus Subfr statistical test
[p, h, stats] = ranksum(fr_stat_data.f_40.sus_fr, fr_stat_data.f_140.sus_fr);
ranksum_FR_sus_sus_stats = struct();
ranksum_FR_sus_sus_stats.p = p
ranksum_FR_sus_sus_stats.h = h;
ranksum_FR_sus_sus_stats.stats = stats;

% Show violin plot on stimulation period for Firing Rate 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Struct to identify each group of data
    fr_stat_data = struct();

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
        data = sign(data).*log10(abs(data));
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'MedianColor', [1 0 0], 'GroupOrder', {'Stim'}, ViolinOpts);
        
        violins(1).ViolinColor = {Multi_func.stim_color};

        hold on;

        % Plot individual lines between violins
        %plot([1, 2], data, '-', 'Color', [0 0 0 0.2]);
        %hold on;
        Multi_func.set_default_axis(gca);
        title([f_stim(3:end) ' Hz DBS']);
        ylabel('Firing Rate Change (Hz)');
        %TODO resize the xaxis here before adding the 10 exnential value
        y = gca; y = y.YAxis;
        Multi_func.set_spacing_axis(y, 4, 1);
        y_sign = arrayfun(@num2str, 10*sign(yticks), 'UniformOutput', 0);
        y_vals = arrayfun(@num2str, abs(yticks), 'UniformOutput', 0);

        ylab = strcat(y_sign, '^{', y_vals, '}');
        zero_idx = find(sign(yticks) == 0);
        if ~isempty(zero_idx)
            ylab{zero_idx} = '0';
        end
        yticklabels(ylab);

        if strcmp(f_region, 'r_M1') == 1
            ylim([-2 2]);
        elseif strcmp(f_region, 'r_combine') == 1
            ylim([-10 70]);
        end

        %fr_stat_data.(f_stim) = struct();
        fr_stat_data.(f_stim).stim_fr = data_bystim.(f_stim).neuron_stim_FR;
    end
    
    sgtitle([f_region ' Firing Rate Violins'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_FR_violin.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_population_comp_stim_FR_violin.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_population_comp_FR_violin.eps']);
end
return;


% Subthreshold Vm showing all DBS pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_Vm = mean(data_bystim.(f_stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;

        % Plot the subthreshold Vm
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
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
        hold on;

        % Plot the Vm scale
        poxs = .2;
        posy = 5;
        Vm_scale = 2;
        plot([posx, posx], [posy, posy + Vm_scale], 'k', 'LineWidth', 2);
        text(posx - .01, posy, [num2str(Vm_scale) ' VM'], 'Rotation', 90);

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]);
        a = gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        set(gca, 'color', 'none')
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average subthreshold Vm Showing all pulses'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_Vm.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_Vm.eps'], 'epsc');
end

% Continuous firing rate showing all DBS pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 9.5, 5.0]);
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate_small, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim{1}).neuron_srate_small, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate_small, 2);
        sem_srate = std_srate./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
        nexttile;

        % Plot the subthreshold srate
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;
        
        % Plot the DBS stimulation time pulses
        stim_time = nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2);
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
        srate_scale = 2;
        plot([posx, posx], [posy, posy + srate_scale], 'k', 'LineWidth', 0.5);
        text(posx - .01, posy, [num2str(srate_scale) ' FR (Hz)'], 'Rotation', 90, 'FontSize', 7);

        % Increase timescale resolution
        xlim([0 - .100, max(stim_time) + 0.100]);
        a = gca;
        a.XAxis.Visible = 'off';
        a.YAxis.Visible = 'off';
        set(gca, 'color', 'none')
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Population Firing Rate Showing All Pulses'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.eps'], 'epsc');
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
