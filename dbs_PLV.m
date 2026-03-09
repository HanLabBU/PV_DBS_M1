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

field1 = fieldnames(region_data);
field1 = field1(1);
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');


%% Loop through and calculate dbs-Vm PLV values for all region and conditions
% Needed for ;single_cell_mod.m'

% The base data here I do not believe will have any relevant data since
% pulses do not occur during the baseline period

freqs = Multi_func.entr_freqs;
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        
        avg_Fs = mean(stim_data.framerate, 'omitnan');

        % Loop through each neuron
        base_plvs = [];
        stim_plvs = [];
        base_phase_vectors = {};
        stim_phase_vectors = {};
        base_plvs_adjusted = [];
        stim_plvs_adjusted = [];
        
        vm_phases = {};

        tic;
        %parfor (nr = 1:length(stim_data.all_trial_SubVm), 0)
        for nr = 1:length(stim_data.all_trial_SubVm)
            % Function that waits for a whole signal to be passed then perform filtering 
            %with constants freq and samp_freq
            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            % Function that waits for a function and matrix and then returns a function waiting for a column value
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            % Aply the filtering function with the whole trial matrix, so all is left is waiting for a column number
            partial_apply = applyFunToColi(filt_trial, stim_data.all_trial_rawVm{nr});
            % Iterate through each column trial and concatenate all o fthe results
            vm_phases{nr} = arrayfun(partial_apply, [1:size(stim_data.all_trial_rawVm{nr}, 2)]' , 'UniformOutput', false); 
            % 
            vm_phases{nr} = cat(3, vm_phases{nr}{:});
        
            time = stim_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(vm_phases{nr}, 3));

            base_idx = find(time < Multi_func.base_ped(2)/1000);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            
            % Calculate stimulation raster points in trace timepoints
            nr_stim_time = reshape(stim_data.stim_timestamps(:, nr), 1, []);
            nr_frame_time = stim_data.trace_timestamps(:, nr);
            diffs = abs(nr_stim_time - nr_frame_time);
            [~, stim_idx_i] = min(diffs, [], 1);
            
            dbs_raster = zeros(size(nr_frame_time));
            dbs_raster(stim_idx_i) = 1;
            dbs_rasters = repmat(dbs_raster, 1, size(vm_phases{nr}, 3));

            base_rasters = zeros(size(dbs_rasters));
            stim_rasters = zeros(size(dbs_rasters));
            
            base_rasters(base_idx) = dbs_rasters(base_idx);
            stim_rasters(stim_idx) = dbs_rasters(stim_idx);
            
            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, base_rasters, 0, 10);   
            base_plvs(:, nr) = PLV(:);
            base_plvs_adjusted(:, nr) = PLV2(:);
            base_phase_vectors{nr} = norm_vecs;

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, stim_rasters, 0, 10);           
            stim_plvs(:, nr) = PLV(:);
            stim_plvs_adjusted(:, nr) = PLV2(:); %TODO I changed this to transpose the dimensions
            stim_phase_vectors{nr} = norm_vecs;
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).base_dbsvm_plvs_adj = base_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_dbsvm_plvs_adj = stim_plvs_adjusted;
        region_data.(f_region).(f_stim).base_dbsvm_plvs = base_plvs;
        region_data.(f_region).(f_stim).stim_dbsvm_plvs = stim_plvs;
        region_data.(f_region).(f_stim).base_dbsvm_phase_vecs = base_phase_vectors;
        region_data.(f_region).(f_stim).stim_dbsvm_phase_vecs = stim_phase_vectors;
        toc
    end
end

%% Calculate each cell's shuffled DBS-Vm PLV for the frequency sweep

% Specify randomization parameters
rng(123);
num_iter = 500;

freqs = Multi_func.entr_freqs;

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        avg_fs = mean(popul_data.framerate, 'omitnan');

        % Struct for all of the shuffled data across frequencies
        all_shuf_plv = struct;
        
        % Loop through each neuron
        tic;
        for nr=1:length(popul_data.neuron_name)
            % Try to do vectorization for the frequency
            vm_phases = {};

            % Function that waits for a whole signal to be passed then perform filtering 
            %with constants freq and samp_freq
            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            % Function that waits for a function and matrix and then returns a function waiting for a column value
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            % Aply the filtering function with the whole trial matrix, so all is left is waiting for a column number
            partial_apply = applyFunToColi(filt_trial, popul_data.all_trial_rawVm{nr});
            % Iterate through each column trial and concatenate all o fthe results
            vm_phases = arrayfun(partial_apply, [1:size(popul_data.all_trial_rawVm{nr}, 2)]' , 'UniformOutput', false);  
            vm_phases = cat(3, vm_phases{:});

            % Create stimulation raster points in trace timepoints
            nr_stim_time = reshape(popul_data.stim_timestamps(:, nr), 1, []);
            nr_frame_time = popul_data.trace_timestamps(:, nr);
            diffs = abs(nr_stim_time - nr_frame_time);
            [~, stim_idx_i] = min(diffs, [], 1);
           
            dbs_raster = zeros(size(nr_frame_time));
            dbs_raster(stim_idx_i) = 1;
            dbs_rasters = repmat(dbs_raster, 1, size(vm_phases, 3));

            shuf_plv_wfreqs = zeros(length(freqs), num_iter);
            shuf_plv_adj_wfreqs = [];

            % Shuffle the DBS pulses and calculate the PLV for each
            % iteration
            parfor(i=1:num_iter, 5)
                % Need to figure out how to randomly place DBS trains
                % within each frequency for each trial

                % Randomly select stim onset time
                onset_rand = randi(last_index, 1, size(vm_phases, 3));
                
                % Reset the starting stim_idx for each trial and keep the spacing between indices the same
                rand_stim_idx = repmat(stim_idx_i' - stim_idx_i(1), 1, size(vm_phases, 3)) + onset_rand;
                    
                col_ind_mat = repmat(1:size(rand_stim_idx, 2), size(rand_stim_idx, 1), 1);
                stim_rasters_wfreqs = zeros(size(dbs_rasters));   
                stim_rasters_wfreqs(sub2ind(size(stim_rasters_wfreqs),...
                                     rand_stim_idx(:),...
                                     col_ind_mat(:))) = 1;

                [sh_PLV_wfreqs, sh_PLV2_wfreqs, ~] = Multi_func.spike_field_PLV(vm_phases, stim_rasters_wfreqs, 0, 10);
                % Store shuffled plvs for all frequencies
                % Rows are the frequencies and columns are the individual
                shuf_plv_wfreqs(:, i) = sh_PLV_wfreqs;
                shuf_plv_adj_wfreqs(:, i) = sh_PLV2_wfreqs;
            
            end
            toc
            % Save the shuffled data
            all_shuf_plv(nr).shuf_plv_wfreqs = shuf_plv_wfreqs;
            all_shuf_plv(nr).shuf_plv_adj_wfreqs = shuf_plv_adj_wfreqs;
            
        end
        
        region_data.(f_region).(f_stim).shuf_plv_data = all_shuf_plv;
    end
end

%% Plot the shuffled distribution with each frequency to check that values are shuffled
freqs = Multi_func.entr_freqs;
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        stim_freq = str2num(f_stim(3:end));
        figure('Position', [100 100 1500 2000]);
        tiledlayout(length(popul_data.neuron_name), 2, 'TileSpacing','tight', 'Padding','loose');
            
        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            nexttile;
            
            histogram(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs([1:stim_freq-1, stim_freq+1:end], :), 1000);
            
            nexttile;
            histogram(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs(stim_freq, :), 1000);
            
%             % Just plot the shuffled in one plot with Y being the frequency
%             f_x = min(freqs):5:max(freqs);
% 
%             tic
%             for i= f_x
%                 nexttile;
%                 histogram(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs(i, :), 1000);
%                 title(num2str(i));
%             end
%             toc
%             % Dont do all 200, just some linspace of them
%             %TODO create histograms maybe?
%             xlabel('PLV');
        end
    end
end


%% Loop and plot all of the dbs-vm stuff
freqs = Multi_func.entr_freqs;
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        stim_freq = str2num(f_stim(3:end));

        figure;
        
        % Baseline PLV does not make sense for DBS
        % Plot base data with error bars
        %base_plvs_mean = nanmean(stim_data.base_dbsvm_plvs_adj, 1);
        %base_plvs_std = nanstd(stim_data.base_dbsvm_plvs_adj, 1);
        %num_base_plvs = size(stim_data.base_dbsvm_plvs_adj, 1);
        %base_plvs_sem = base_plvs_std./sqrt(num_base_plvs);
 
        %fill_h = fill([freqs, flip(freqs)], [[base_plvs_mean + base_plvs_sem], flip(base_plvs_mean - base_plvs_sem)], [0.5 0.5 0.5]);
        %Multi_func.set_fill_properties(fill_h);
        %hold on;
        %plot(freqs, base_plvs_mean, 'color', Multi_func.base_color);
        %hold on;
        %fill_h = fill([freqs, flip(freqs)], [[base_plvs_mean + base_plvs_sem], flip(base_plvs_mean - base_plvs_sem)], [0.5 0.5 0.5]);
        %Multi_func.set_fill_properties(fill_h);
        %hold on;

        % Plot the shuffled distribution for each frequency
        %shuf_plvs = arrayfun(@(x) disp(x.shuf_plv_adj_wfreqs), stim_data.shuf_plv_data, ...
        %    'UniformOutput',false);


        shuf_plvs = arrayfun(@(x) mean(x.shuf_plv_adj_wfreqs, 2), stim_data.shuf_plv_data, ...
            'UniformOutput',false);
        shuf_plvs = cat(2, shuf_plvs{:});
        
        shuf_plvs_mean = nanmean(shuf_plvs, 2)';
        shuf_plvs_std = nanstd(shuf_plvs, [], 2)';
        num_shuf_plvs = size(shuf_plvs, 2);
        shuf_plvs_sem = shuf_plvs_std./sqrt(num_shuf_plvs);
 
        fill_h = fill([freqs, flip(freqs)], [[shuf_plvs_mean + shuf_plvs_sem], flip(shuf_plvs_mean - shuf_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(freqs, shuf_plvs_mean, 'color', Multi_func.shuf_color);
        hold on;
        fill_h = fill([freqs, flip(freqs)], [[shuf_plvs_mean + shuf_plvs_sem], flip(shuf_plvs_mean - shuf_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;

        % Plot stim data with error bars
        stim_plvs_mean = nanmean(stim_data.stim_dbsvm_plvs_adj, 2)';
        stim_plvs_std = nanstd(stim_data.stim_dbsvm_plvs_adj, [], 2)';
        num_stim_plvs = size(stim_data.stim_dbsvm_plvs_adj, 2);
        stim_plvs_sem = stim_plvs_std./sqrt(num_stim_plvs);
        
        fill_h = fill([freqs, flip(freqs)], [[stim_plvs_mean + stim_plvs_sem], flip(stim_plvs_mean - stim_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, stim_plvs_mean, 'color', Multi_func.stim_color);
        hold on;
        
        ax = gca;
        
        %set(ax,'Xscale','log');
        xlabel('Frequency (Hz)');
        ylabel('Pulse-Vm PLV^2');
        
        %legend({'base', 'stim'}, 'Location', 'west');
        ax.Units = 'centimeters';
        %ax.InnerPosition = [2 2 3.91 3.24];% Old dimensions
        ax.InnerPosition = [2 2 3 3];
        
        
        % Test significance between shuffled and observed for stimulation
        % frequency
        p = signrank(stim_data.stim_dbsvm_plvs_adj(stim_freq, :), shuf_plvs(stim_freq, :));

        ylim([-0.07 1]);
        title([f_stim(3:end) ' ' f_region(3:end) ' Hz p=' num2str(p)], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'PLV' f 'PLV_dbsvm_' f_region '_' f_stim '.png']);
        saveas(gcf, [figure_path 'PLV' f 'PLV_dbsvm_' f_region '_' f_stim '.pdf']);
        
        % Grab single cell
        %num_neurons = size(stim_data.base_plvs_adjusted, 1);
        %tiledlayout(num_neurons, 1);
        %for nr = 1:num_neurons
        %    figure;
        %    %nexttile;
        %    plot(stim_data.base_plvs_adjusted(nr, :)', 'b');
        %    hold on;
        %    plot(stim_data.stim_plvs_adjusted(nr, :)', 'g');
        %    title(stim_data.neuron_name{nr}, 'Interpreter', 'none');
        %end 
    end
end

%% Plot population roseplot of the phases vectors from DBS-Vm and the shuffled distribution
% Note: requires 'single cell DBS-PLV' section executed from 'single_cell_mod.m'

%TODO need to finish implementing this function

%nr_pop = 'all';
nr_pop = 'etrain';
%nr_pop = 'non_entr';

remove_nonmod_nrs = 1;

polar_edges = linspace(0, 2*pi, 24);

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Setup figure to show frequencies side-by-side
    figure;
    tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);


    % Loop through stim frequencies
    for f_stim = stims'
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

        
        nexttile;
        % Plot the stim phase-vectors
        %stim_phases = 
        polarhistogram(base_phases, polar_edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.base_color);
        title([f_stim]);
    end
end

%% Plot individual neuron roseplot of the phases vectors from DBS-Vm and the shuffled distribution
% Note: requires 'single cell DBS-PLV' section executed from 'single_cell_mod.m'

polar_edges = linspace(0, 2*pi, 24);

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Skip CA1 region
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Setup figure to show frequencies side-by-side
    %figure;
    %tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);


    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Loop through each neuron name
        for nr=1:length(popul_data.neuron_name)
            
            figure;
            
            % Grab the stim-phases of the stimulation frequency
            stim_phases = popul_data.stim_dbsvm_phase_vecs{nr};
            stim_phases = stim_phases(:, find(Multi_func.entr_freqs == str2num(f_stim(3:end))));
            
            % Grab the shuffled phase vectors
            shuf_phases = popul_data.plv_mod_stats(nr).shuf_phase_vecs;
            shuf_phases = cat(2, shuf_phases{:});
            
            % Save the PLV mod stats
            nr_plv_mod = popul_data.plv_mod_stats(nr).mod;

            nexttile;
            polarhistogram(angle(stim_phases), polar_edges, 'Normalization', 'probability', 'EdgeColor', Multi_func.stim_color);
            title('Stim Phases');

            nexttile;
            polarhistogram(angle(shuf_phases), polar_edges, 'Normalization', 'probability', 'EdgeColor', Multi_func.base_color)
            title('Shuffled Phases');

            sgtitle([f_region(3:end) f_stim(3:end) ' mod: ' num2str(nr_plv_mod)], 'Interpreter', 'none');

            saveas(gcf, [figure_path 'Phase' f 'Neuronwise' f 'DBS_Vm_Phase_roseplot_' popul_data.neuron_name{nr} '.png']);
            saveas(gcf, [figure_path 'Phase' f 'Neuronwise' f 'DBS_Vm_Phase_roseplot_' popul_data.neuron_name{nr} '.pdf']);        
        end
    end
end

% Just want to close all figures
close all;

%% Combine the DBS-PLV for both CA1 and M1 stim data
test_regions = {'r_CA1', 'r_M1'};
stims = fieldnames(region_data.r_M1)';
freqs = [1:200];
for f_stim=stims'
    f_stim = f_stim{1};
    
    figure();

    for f_region = test_regions
        f_region = f_region{1};
        stim_data = region_data.(f_region).(f_stim);
        
        stim_plvs_mean = nanmean(stim_data.stim_dbsvm_plvs_adj, 1);
        stim_plvs_std = nanstd(stim_data.stim_dbsvm_plvs_adj, 1);
        num_stim_plvs = size(stim_data.stim_dbsvm_plvs_adj, 1);
        stim_plvs_sem = stim_plvs_std./sqrt(num_stim_plvs);
        
        fill_h = fill([freqs, flip(freqs)], [[stim_plvs_mean + stim_plvs_sem], flip(stim_plvs_mean - stim_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        Multi_func.set_default_axis(gca);
        hold on;

        if strcmp(f_region, 'r_CA1') == 1
            cur_color = Multi_func.ca1_color;
            label = 'CA1';
        elseif strcmp(f_region, 'r_M1') == 1
            cur_color = Multi_func.m1_color;
            label = 'M1';
        end
        plot(freqs, stim_plvs_mean, 'color', cur_color, 'DisplayName', label);
        hold on;

    end
    legend();  
    ax = gca;
    set(ax,'Xscale','log');
    title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
end
