clear all;
close all;
clc;

%%User Modification
f = filesep;
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

% Flag to use low frequency oscillation
%TODO what was I trying to do here??
use_low_freq = 0;
use_plv_low = 1;

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
end
%Load the data
load(save_all_data_file);


if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

field1 = {'r_M1'};
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');

%% Loop through and calculate spike-Vm PLV values for all region and conditions
freqs = [1:200];

% Flag to use low frequency oscillation
use_low_freq = 0;
%use_plv_low = 1;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        base_plvs = [];
        stim_plvs = [];
        base_phase_vectors = [];
        stim_phase_vectors = [];
        base_plvs_adjusted = [];
        stim_plvs_adjusted = [];
        
        tic;
        
        % Parallel array that indicates if it has large 2-10Hz plv in baseline
        has_plv_low = [];
        vm_phases = {};

        for nr = 1:length(popul_data.all_trial_SubVm)
            
            % Skip if neuron does not have low frequency
            if use_low_freq && data_bystim.(f_stim).has_low_freq(nr) == 0
                disp('Skipping neuron ');
                continue;
            end


            %TODO this filt_trial for frequency and timeline has the arrays backwards
            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            partial_apply = applyFunToColi(filt_trial, popul_data.all_trial_SubVm{nr});
            
            % The rows are the frequency and columns are the time idx, 3rd are individual trials
            vm_phases{nr} = arrayfun(partial_apply, [1:size(popul_data.all_trial_SubVm{nr}, 2)]' , 'UniformOutput', false); 
            % 
            vm_phases{nr} = cat(3, vm_phases{nr}{:});
        
            time = popul_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(vm_phases{nr}, 3));

            base_idx = find(time < Multi_func.base_ped(2)/1000);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            
            spike_rasters = popul_data.all_trial_spike_rasters{nr};
            base_rasters = zeros(size(spike_rasters));
            stim_rasters = zeros(size(spike_rasters));
            
            base_rasters(base_idx) = spike_rasters(base_idx);
            stim_rasters(stim_idx) = spike_rasters(stim_idx);
            
            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, base_rasters, 7, 10);             
            base_plvs(nr, :) = PLV;
            base_plvs_adjusted(nr, :) = PLV2;
            base_phase_vectors = [base_phase_vectors; norm_vecs];

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, stim_rasters, 7, 10);           
            stim_plvs(nr, :) = PLV;
            stim_plvs_adjusted(nr, :) = PLV2;
            stim_phase_vectors = [stim_phase_vectors; norm_vecs];

            % TODO perform if it is important 
            %Determine if the neuron has baseline PLV to 2-10Hz
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).base_spikevm_plvs_adj = base_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_spikevm_plvs_adj = stim_plvs_adjusted;
        region_data.(f_region).(f_stim).base_spikevm_plvs = base_plvs;
        region_data.(f_region).(f_stim).stim_spikevm_plvs = stim_plvs;
        region_data.(f_region).(f_stim).base_spikevm_phase_vecs = base_phase_vectors;
        region_data.(f_region).(f_stim).stim_spikevm_phase_vecs = stim_phase_vectors;
        toc
    end
end

%TODO make a way of identifying neurons based on spike-Vm PLV baseline at 2-10Hz
% ensure that use_low_freq == 0 so that other neurons are not skipped

%% Loop and plot all of the spike-vm stuff

% Flag to determine which populations to plot
% The variable must be set from 'single_cell_mod'
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

%low_freqs = [1:50]; % Deprecated??

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
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

        % Setup figure with set dimensions
        figure('Renderer', 'painters', 'Position', [5, 5, 1000, 1000]);

        ax = gca;
        ax.Units = 'centimeters';
        ax.InnerPosition = [2 2 5 4.0];
        
        % Plot base data with error bars
        base_plvs_mean = nanmean(popul_data.base_spikevm_plvs_adj(nr_idxs, :), 1);
        base_plvs_std = nanstd(popul_data.base_spikevm_plvs_adj(nr_idxs, :), 1);
        num_base_plvs = size(popul_data.base_spikevm_plvs_adj(nr_idxs, :), 1);
        base_plvs_sem = base_plvs_std./sqrt(num_base_plvs);
 
        fill_h = fill([freqs, flip(freqs)], [[base_plvs_mean + base_plvs_sem], flip(base_plvs_mean - base_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(freqs, base_plvs_mean, 'color', Multi_func.base_color);
        hold on;
        fill_h = fill([freqs, flip(freqs)], [[base_plvs_mean + base_plvs_sem], flip(base_plvs_mean - base_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;

        % Plot stim data with error bars
        stim_plvs_mean = nanmean(popul_data.stim_spikevm_plvs_adj(nr_idxs, :), 1);
        stim_plvs_std = nanstd(popul_data.stim_spikevm_plvs_adj(nr_idxs, :), 1);
        num_stim_plvs = size(popul_data.stim_spikevm_plvs_adj(nr_idxs, :), 1);
        stim_plvs_sem = stim_plvs_std./sqrt(num_stim_plvs);
        
        fill_h = fill([freqs, flip(freqs)], [[stim_plvs_mean + stim_plvs_sem], flip(stim_plvs_mean - stim_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, stim_plvs_mean, 'color', Multi_func.stim_color);
        hold on;
        
        % Log frequency scale
        %set(gca,'Xscale','log');
        %xlabel('Frequency Log (Hz)');
        
        % Linear frequency scale
        xlabel('Frequency (Hz)');

        ylabel('Spike-Vm PLV^2');
        ylim([-0.07 0.5]);
        title([f_region(3:end) ' ' f_stim(3:end) ' ' nr_pop], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'PLV' f 'PLV_spike_vm_' f_region '_' nr_pop '_' f_stim '.png']);
        saveas(gcf, [figure_path 'PLV' f 'PLV_spike_vm_' f_region '_' nr_pop '_' f_stim '.pdf']);
        
        % Zoom in on the theta range and plot that
        %xlim([2, 10]);
        %
        %saveas(gcf, [figure_path 'PLV' f 'Theta_PLV_spike_vm_' f_region '_' nr_pop '_' f_stim '.png']);
        %saveas(gcf, [figure_path 'PLV' f 'Theta_PLV_spike_vm_' f_region '_' nr_pop '_' f_stim '.pdf']);
        
        % Grab single cell
        %num_neurons = size(popul_data.base_plvs_adjusted, 1);
        %tiledlayout(num_neurons, 1);
        %for nr = 1:num_neurons
        %    figure;
        %    %nexttile;
        %    plot(popul_data.base_plvs_adjusted(nr, :)', 'b');
        %    hold on;
        %    plot(popul_data.stim_plvs_adjusted(nr, :)', 'g');
        %    title(popul_data.neuron_name{nr}, 'Interpreter', 'none');
        %end 
    end
end

%% Save significant spike-Vm PLV of stim 140 Hz and 40 Hz compared to baseline
% Use Wilcoxon's sign rank test 

freqs = [1:200];

% Determine the neuron population
%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non';

stats_log = [figure_path 'PLV' f 'spikevm_plv_stats' nr_pop];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
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

        
        % Grab the current stimulation into a numer variable
        cur_stim = str2num(f_stim(3:end));

        % Store the baseline neuron PLVs
        base_plv = popul_data.base_spikevm_plvs_adj(nr_idxs, cur_stim);
        stim_plv = popul_data.stim_spikevm_plvs_adj(nr_idxs, cur_stim);
        
        [p_sign, ~, stats] = signrank(base_plv, stim_plv);
    
        % Print the stats data to the file
        diary on;
        fprintf('\n\n');
        disp([f_region ' ' f_stim]);
        disp(['p_value: ' num2str(p_sign)]);
        disp(['Statistic: ' num2str(stats.signedrank)]);
        diary off;
        
        % DEBUG
        figure;
        plot(ones(size(base_plv)), base_plv, '.');
        hold on;
        plot(2*ones(size(stim_plv)), stim_plv, '.');
        title([f_region ' ' f_stim], 'Interpreter', 'none');
        xlim([0 3]);
        ylim([0 1]);

        if strcmp(f_region, 'r_M1') && strcmp(f_stim, 'f_140')
            disp('op');
            error('asd');
        end
    end
end

%% Bootstrapping all of the spike-Vm PLVs across all narrow-band frequencies
% Try the boostrapping method here
num_iter = 1000; % Number of estimates for boostrapping
nspikes = 10; % This parameter indicates how many spikes to pool for a phase locking sample
num_w = 12;

freqs = [1:200];

% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        boot_base_plvs = [];
        boot_stim_plvs = [];
        
        %First test for 140Hz
        orig_base_plvs = popul_data.base_spikevm_phase_vecs;
        orig_stim_plvs = popul_data.stim_spikevm_phase_vecs;
        for freq = 1:length(freqs)
            % Check to make sure the window is less than the total number of spikes
            num_spikes = size(orig_base_plvs, 1);
            for id=1:num_iter
                spike_shuf_idx = randperm(num_spikes);
                spikes_select = spike_shuf_idx(1:nspikes);
                all_freq_plvs = abs(nanmean(orig_base_plvs(spikes_select, :), 1) );
                boot_base_plvs(id, freq) = nanmean(all_freq_plvs(freq));
            end
        end
        
        for freq = 1:length(freqs)
            % Check to make sure the window is less than the total number of spikes
            num_spikes = size(orig_stim_plvs, 1);
            for id=1:num_iter
                spike_shuf_idx = randperm(num_spikes);
                spikes_select = spike_shuf_idx(1:nspikes);
                all_freq_plvs = abs(nanmean(orig_stim_plvs(spikes_select, :), 1) );
                boot_stim_plvs(id, freq) = nanmean(all_freq_plvs(freq));
            end
        end
        region_data.(f_region).(f_stim).boot_stim_spikevm_plvs = boot_stim_plvs;
        region_data.(f_region).(f_stim).boot_base_spikevm_plvs = boot_base_plvs;
    end
end

%% Plot THETA quantifications from bootstrapped narrow band frequencies, essentially a PLV value for each frequency
stats_log = [figure_path 'PLV' f 'Theta_narrowband_avg_spike_vm_plv_violin_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;

theta_freqs = [2:10];
% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Violin plots for the average theta
        boot_base_theta_plvs = nanmean(popul_data.boot_base_spikevm_plvs(:, theta_freqs), 2);
        boot_stim_theta_plvs = nanmean(popul_data.boot_stim_spikevm_plvs(:, theta_freqs), 2);

        figure;
        ax = gca;
        ax.Units = 'Centimeters';
        ax.InnerPosition = [5 5 4 5.387];
        num_base = length(boot_base_theta_plvs)
        num_stim = length(boot_stim_theta_plvs)
        data = [boot_base_theta_plvs; boot_stim_theta_plvs];
        labels = [repmat({'Base'}, num_base, 1); repmat({'Stim'}, num_stim, 1)];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        Multi_func.set_default_axis(gca);
        title(['Theta spike-Vm PLVs ' f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'PLV' f 'Theta_spike_vm_boot_narrowband_violin ' f_region(3:end) ' ' f_stim(3:end) '.pdf']);
        saveas(gcf, [figure_path 'PLV' f 'Theta_spike_vm_boot_narrowband_violin ' f_region(3:end) ' ' f_stim(3:end) '.png']);

        diary on;
        disp(['Theta spike-Vm PLVs ' f_region(3:end) ' ' f_stim(3:end)]);
        [p, h, stats] = signtest(boot_base_theta_plvs, boot_stim_theta_plvs)
        diary off;
    end
end


%% Loop through and calculate broad-band THETA spike-Vm PLV values for all region and conditions
freqs = [2, 10]; 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        base_plvs = [];
        stim_plvs = [];
        base_plvs_adjusted = [];
        stim_plvs_adjusted = [];
        base_theta_phase_vectors = [];
        stim_theta_phase_vectors = [];
        
        tic;
        for nr = 1:length(popul_data.all_trial_SubVm)
            % Skip if neuron does not have low frequency
            if use_low_freq && data_bystim.(f_stim).has_low_freq(nr) == 0
                disp('Skippped-------------------');
                continue;
            end
            
            filt_trial = @(trial) (angle(Multi_func.filt_range(trial, freqs, avg_Fs)));
            
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            partial_apply = applyFunToColi(filt_trial, popul_data.all_trial_SubVm{nr});
            
            vm_phases{nr} = arrayfun(partial_apply, [1:size(popul_data.all_trial_SubVm{nr}, 2)]' , 'UniformOutput', false); 
            % 
            vm_phases{nr} = cat(3, vm_phases{nr}{:});
        
            time = popul_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(vm_phases{nr}, 3));

            base_idx = find(time < Multi_func.base_ped(2)/1000);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            
            spike_rasters = popul_data.all_trial_spike_rasters{nr};
            base_rasters = zeros(size(spike_rasters));
            stim_rasters = zeros(size(spike_rasters));
            
            base_rasters(base_idx) = spike_rasters(base_idx);
            stim_rasters(stim_idx) = spike_rasters(stim_idx);
            
            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, base_rasters, 7, 10);             
            base_plvs(nr, :) = PLV;
            base_plvs_adjusted(nr, :) = PLV2;
            base_theta_phase_vectors = [base_theta_phase_vectors; norm_vecs];

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, stim_rasters, 7, 10);           
            stim_plvs(nr, :) = PLV;
            stim_plvs_adjusted(nr, :) = PLV2;
            stim_theta_phase_vectors = [stim_theta_phase_vectors; norm_vecs];
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).base_spikevm_theta_plvs_adj = base_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_spikevm_theta_plvs_adj = stim_plvs_adjusted;

        region_data.(f_region).(f_stim).base_spikevm_theta_plvs = base_plvs;
        region_data.(f_region).(f_stim).stim_spikevm_theta_plvs = stim_plvs;
        region_data.(f_region).(f_stim).base_spikevm_theta_phase_vecs = base_theta_phase_vectors;
        region_data.(f_region).(f_stim).stim_spikevm_theta_phase_vecs = stim_theta_phase_vectors;
        toc
    end
end

%% Loop and plot all of the broadband THETA spike-vm stuff (non-bootstrapped)
stats_log = [figure_path 'PLV' f 'Theta_broadband_spike_vm_plv_nonboot_violin_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        figure;
        num_base = length(popul_data.base_spikevm_theta_plvs_adj)
        num_stim = length(popul_data.stim_spikevm_theta_plvs_adj)
        data = [popul_data.base_spikevm_theta_plvs_adj; popul_data.stim_spikevm_theta_plvs_adj];
        labels = [repmat({'Base'}, num_base, 1); repmat({'Stim'}, num_stim, 1)];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');

        % Grab single cell
        %num_neurons = size(popul_data.base_plvs_adjusted, 1);
        %tiledlayout(num_neurons, 1);
        %for nr = 1:num_neurons
        %    figure;
        %    %nexttile;
        %    plot(popul_data.base_plvs_adjusted(nr, :)', 'b');
        %    hold on;
        %    plot(popul_data.stim_plvs_adjusted(nr, :)', 'g');
        %    title(popul_data.neuron_name{nr}, 'Interpreter', 'none');
        %end 
    end
end

%% Bootstrapping broadband THETA spike-Vm PLVs 
% Try the boostrapping method here
num_iter = 500; % Number of estimates for boostrapping
nspikes = 10; % This parameter indicates how many spikes to pool for a phase locking sample
num_w = 12;

% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        boot_base_plvs = [];
        boot_stim_plvs = [];
        

        orig_base_plvs = popul_data.base_spikevm_theta_phase_vecs;
        orig_stim_plvs = popul_data.stim_spikevm_theta_phase_vecs;

        num_spikes = size(orig_base_plvs, 1);
        for id=1:num_iter
            spike_shuf_idx = randperm(num_spikes);
            spikes_select = spike_shuf_idx(1:nspikes);
            all_freq_plvs = abs(nanmean(orig_base_plvs(spikes_select, :), 1) );
            boot_base_plvs(id) = nanmean(all_freq_plvs);
        end
        
        % Check to make sure the window is less than the total number of spikes
        num_spikes = size(orig_stim_plvs, 1);
        for id=1:num_iter
            spike_shuf_idx = randperm(num_spikes);
            spikes_select = spike_shuf_idx(1:nspikes);
            all_freq_plvs = abs(nanmean(orig_stim_plvs(spikes_select, :), 1) );
            boot_stim_plvs(id) = nanmean(all_freq_plvs);
        end

        region_data.(f_region).(f_stim).boot_stim_theta_spikevm_plvs = boot_stim_plvs;
        region_data.(f_region).(f_stim).boot_base_theta_spikevm_plvs = boot_base_plvs;
    end
end

%% Plot THETA quantifications from bootstrapped broadband spike-vm PLV
stats_log = [figure_path 'PLV' f 'Theta_broadband_spike_vm_plv_boot_violin_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;

% Loop through each region
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Violin plots of the bootstrapped theta
        %boot_base_theta_plvs = nanmean(popul_data.boot_base_spikevm_plvs(:, theta_freqs), 2);
        %boot_stim_theta_plvs = nanmean(popul_data.boot_stim_spikevm_plvs(:, theta_freqs), 2);
        boot_base_theta_plvs = popul_data.boot_base_theta_spikevm_plvs;
        boot_stim_theta_plvs = popul_data.boot_stim_theta_spikevm_plvs;

        figure;
        ax = gca;
        ax.Units = 'Centimeters';
        ax.InnerPosition = [5 5 4 5.387];
        num_base = length(boot_base_theta_plvs)
        num_stim = length(boot_stim_theta_plvs)
        data = [boot_base_theta_plvs; boot_stim_theta_plvs];
        labels = [repmat({'Base'}, num_base, 1); repmat({'Stim'}, num_stim, 1)];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        Multi_func.set_default_axis(gca);
        title(['Theta spike-Vm PLVs ' f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'PLV' f 'Theta_spike_vm_boot_broadband_violin ' f_region(3:end) ' ' f_stim(3:end) '.pdf']);
        saveas(gcf, [figure_path 'PLV' f 'Theta_spike_vm_boot_broadband_violin ' f_region(3:end) ' ' f_stim(3:end) '.png']);

        diary on;
        disp(['Theta spike-Vm PLVs ' f_region(3:end) ' ' f_stim(3:end)]);
        [p, h, stats] = signtest(boot_base_theta_plvs, boot_stim_theta_plvs)
        diary off;
    end
end


%% Loop through and calculate spike-DBS PLV values for all region and conditions
freqs = [1:200];
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        stim_plvs = [];
        stim_plvs_adjusted = [];
        stim_phase_vecs = [];
        
        disp(f_stim);
        tic;
        for nr = 1:length(popul_data.all_trial_SubVm)
            % Creates a dbs signal just using a one at pulse times
            dbs_signal = zeros(size(popul_data.trace_timestamps(:, nr)));
            pulse_time = popul_data.stim_timestamps(:, nr);
            dbs_start = popul_data.stim_timestamps(1, nr) - popul_data.trace_timestamps(1, nr);
            stim_idx =  ceil(((dbs_start + pulse_time)*avg_Fs) + 1); 
            dbs_signal(stim_idx) = 1;
            dbs_signal = zscore(dbs_signal);
            dbs_signal = repmat(dbs_signal, 1, size(popul_data.all_trial_SubVm{nr}, 2));
            
            %pulse_time = popul_data.stim_timestamps(:, nr);
            %period = 2*pi/(mean(diff(pulse_time)));
            %dbs_signal = zeros(size(popul_data.trace_timestamps(:, nr)));
            %time = popul_data.trace_timestamps(:, nr);
            %stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            %dbs_signal(stim_idx) = cos(period*time(stim_idx));
            %dbs_signal = zscore(dbs_signal);
            %dbs_signal = repmat(dbs_signal, 1, size(popul_data.all_trial_SubVm{nr}, 2));


            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            partial_apply = applyFunToColi(filt_trial, dbs_signal);
            

            dbs_phases{nr} = arrayfun(partial_apply, [1:size(dbs_signal, 2)]' , 'UniformOutput', false); 
            % 
            dbs_phases{nr} = cat(3, dbs_phases{nr}{:});
        
            time = popul_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(dbs_phases{nr}, 3));
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);

            spike_rasters = popul_data.all_trial_spike_rasters{nr};
            stim_rasters = zeros(size(spike_rasters));
            stim_rasters(stim_idx) = spike_rasters(stim_idx);

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(dbs_phases{nr}, stim_rasters, 7, 10);             

            stim_plvs(nr, :) = PLV;
            stim_plvs_adjusted(nr, :) = PLV2;
            stim_phase_vecs = [stim_phase_vecs; norm_vecs];

            %DEBUG
            %figure;
            %plot(PLV2);
            %title(f_stim, 'Interpreter', 'none');

            %DEBUG
            %figure;
            %spike_idx = find(stim_rasters(:, 1) == 1);
            %plot(spike_idx, stim_rasters(spike_idx, 1), '|');
            %hold on;
            %plot(dbs_signal(:, 1));
            %title(f_stim, 'Interpreter', 'none');
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).stim_spikedbs_plvs_adj = stim_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_spikedbs_plvs = stim_plvs;
        region_data.(f_region).(f_stim).stim_spikedbs_phase_vecs = stim_phase_vecs;
        toc
    end
end

%% Loop and plot all of the spike-dbs stuff
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        figure;
        plvs_mean = nanmean(popul_data.stim_spikedbs_plvs_adj, 1);
        plvs_std = nanstd(popul_data.stim_spikedbs_plvs_adj, 1);
        num_plvs = size(popul_data.stim_spikedbs_plvs_adj, 1);
        plvs_sem = plvs_std./sqrt(num_plvs);
        
        yline(0);
        hold on;
        fill_h = fill([freqs, flip(freqs)], [[plvs_mean + plvs_sem], flip(plvs_mean - plvs_sem)], [0.5 0.5 0.5]);
        hold on;
        plot(freqs, plvs_mean, 'color', Multi_func.stim_color);
        Multi_func.set_fill_properties(fill_h);
        Multi_func.set_default_axis(gca);

        legend({'stim'});
        set(gca,'Xscale','log');

        %ylim([-0.07 0.5]);
        %xlim([2, 10]);
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'PLV' f 'PLV_DBSStim_' f_region '_' f_stim '.png']);
        saveas(gcf, [figure_path 'PLV' f 'PLV_DBSStim_' f_region '_' f_stim '.pdf']);
        


        % Plot all of the PLVs
        figure;
        plot(freqs, popul_data.stim_spikedbs_plvs_adj);
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
        num_neurons = sum(~isnan(popul_data.stim_spikedbs_plvs_adj(:, 1)));
        legend(num2str(num_neurons));


        % Plot the spike population wide PLV
        figure;
        num_spikes = size(popul_data.stim_spikedbs_phase_vecs, 1);
        PLV = (1/num_spikes)*abs(sum(popul_data.stim_spikedbs_phase_vecs, 1));
        PLV_adj = (1/((num_spikes - 1))).*((PLV.^2).*num_spikes - 1); 
        plot(freqs, PLV_adj);
        ylim([-0.1 0.8]);
        title([f_region(3:end) ' ' f_stim(3:end) ' from pop spikes'], 'Interpreter', 'none');
        legend(num2str(num_spikes));

        % Grab single cell
        %num_neurons = size(popul_data.stim_spikedbs_plvs_adj, 1);
        %for nr = 1:num_neurons
        %    figure;
        %    plot(popul_data.stim_spikedbs_plvs_adj(nr, :)', 'g');
        %    title(popul_data.neuron_name{nr}, 'Interpreter', 'none');
        %end 
    end
end

%% Bootstrapping all of the spike-dbs PLVs across all frequencies
% Try the boostrapping method here
num_iter = 1000; % Number of estimates for boostrapping
nspikes = 10; % This parameter indicates how many spikes to pool for a phase locking sample
num_w = 12;

freqs = [1:200];

% Loop through each region
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        boot_stim_plvs = [];
        
        orig_stim_plvs = popul_data.stim_spikedbs_phase_vecs;
        
        for freq = 1:length(freqs)
            % Check to make sure the window is less than the total number of spikes
            num_spikes = size(orig_stim_plvs, 1);
            for id=1:num_iter
                spike_shuf_idx = randperm(num_spikes);
                spikes_select = spike_shuf_idx(1:nspikes);
                all_freq_plvs = abs(nanmean(orig_stim_plvs(spikes_select, :), 1) );
                boot_stim_plvs(id, freq) = nanmean(all_freq_plvs(freq));
            end
        end
        region_data.(f_region).(f_stim).boot_stim_spikedbs_plvs = boot_stim_plvs;
    end
end


%% DEBUG compare boostrapped with original data
figure;
plot(region_data.(f_region).(f_stim).base_spikevm_plvs_adj');
title('Original Base Data');

figure;
plot(boot_base_plvs');
title('Boot base Data');

figure;
plot(region_data.(f_region).(f_stim).stim_spikevm_plvs_adj');
title('Original Stim Data');

figure;
plot(boot_stim_plvs');
title('Boot Stim Data');


figure; 
plot(nanmean(boot_stim_plvs, 1), 'g');
hold on;
plot(nanmean(boot_base_plvs, 1), 'b');

%% Tried to vectorize with funciton handles to make PLV calculation easier
%
%%Store PLVs
%spike_vm_plvs = zeros();
%parfor
%% Single trial call
%%[plv, plv_adjusted] = Multi_func.plv_calc(phase_signals, freqs, spikes);
%calc_plv = @(trial) (Multi_func.plv_calc(Multi_func.filt_data(trial, frs, Fs)));
%% Anonymous function to apply 
%applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
%partial_apply = applyFunToCols(calc_plv, nr);
%neuron_plv = @(nr) arrafun(partial_apply, 1:size(nr, 2)));
%
%cellfun(neuron_plv, popul_data.all_trial_SubVm);
