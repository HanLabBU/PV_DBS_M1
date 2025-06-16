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

% Be careful if you current amplitude data saved
%%% Only use when other variables are saved to the region data
%save(save_all_data_file, 'region_data', '-v7.3');


%% Setup the color variables
[red_blue_color_cmap] = (cbrewer('div', 'RdBu',500));
red_blue_color_cmap(red_blue_color_cmap > 1) = 1;
red_blue_color_cmap(red_blue_color_cmap < 0) = 0;
red_blue_color_cmap = flipud(red_blue_color_cmap);

%% Calculate the Vm modulation stats from onset, transient period
% Creates 'Vm_trans_mod_stats' field
ped_start = Multi_func.trans_ped(1)./1000; % Convert to sec
ped_stop = Multi_func.trans_ped(2)./1000;

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
            
            % Old way
            %base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            %stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            
            base_idx = find(nr_timeline > 0 - ped_stop & ...
                            nr_timeline < 0 - ped_start );

            stim_idx = find(nr_timeline > 0 + ped_start & ...
                            nr_timeline < 0 + ped_stop );

            % Calculate the average amplitude Vm for each trial
            Vm_mean_pre = mean(popul_data.all_trial_rawVm{nr}(base_idx, :), 1);
            Vm_mean_stim = mean(popul_data.all_trial_rawVm{nr}(stim_idx, :), 1);

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

%% Calculate the Vm modulation stats for sustained period
% Creates 'Vm_sus_mod_stats' field

ped_start = Multi_func.sus_ped(1)/1000; % convert to sec
ped_stop = Multi_func.sus_ped(2)/1000;

% For debugging purposes
show_figures = 1;

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

            % Old way
            %base_idx = find(nr_timeline > 0 - range([per_start, per_stop]) & nr_timeline < 0);
            %stim_idx = find(nr_timeline > per_start & nr_timeline < per_stop);
            
            % Reflect the stimulation time period across the onset
            base_idx = find(nr_timeline > 0 - ped_stop & ...
                            nr_timeline < 0 - ped_start )

            stim_idx = find(nr_timeline > 0 + ped_start & ...
                            nr_timeline < 0 + ped_stop );

            % Calculate the average amplitude Vm for each trial
            Vm_mean_pre = mean(popul_data.all_trial_rawVm{nr}(base_idx, :), 1);
            Vm_mean_stim = mean(popul_data.all_trial_rawVm{nr}(stim_idx, :), 1);
            
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

            % Debug the figures
            if show_figures == 1 && strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_140') == 1 && nr == 24
                figure;
                plot(popul_data.all_trial_rawVm{nr} + 50*[1:size(popul_data.all_trial_rawVm{nr}, 2)]);
                legend(['p=' num2str(p_sign)]);
                sgtitle(popul_data.neuron_name{nr}, 'Interpreter', 'none');
                
                %TODO something with the base_idx??????
                %figure;
                % %plot(nr_timeline, '.');
                %plot(base_idx, '.')

                %figure;
                %plot(stim_idx, '.')
                %popul_data.all_trial_rawVm{nr}(base_idx, :)
            end
        end
        
        %DEBUG
        if strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_140') == 1
            %disp(sum([Vm_sus_mod_stats.mod] == -1, "all"));
            disp(find([Vm_sus_mod_stats.mod] == -1));
            %e rror('Pause');
        end
        
        region_data.(f_region).(f_stim).Vm_sus_mod_stats = Vm_sus_mod_stats;
    end
end

%% Calculate modulation of spike rate in the transient period
% Creates 'fr_trans_mod_stats' field

ped_start = Multi_func.trans_ped(1)./1000; % Convert to sec
ped_stop = Multi_func.trans_ped(2)./1000;

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
           
            % Old way
            %base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            %stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            
            base_idx = find(nr_timeline > 0 - ped_stop & ...
                            nr_timeline < 0 - ped_start);

            stim_idx = find(nr_timeline > 0 + ped_start & ...
                            nr_timeline < 0 + ped_stop);

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

%% Calculate modulation of spike rate in the sustained period
% Creates 'fr_sus_mod_stats' field

ped_start = Multi_func.sus_ped(1)./1000; % Concert to Sec
ped_stop = Multi_func.sus_ped(2)./1000;

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
            
            % Old way but it is wrong!!!!!!!
            %base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            %stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);
            % -- End of wrong way
                
            base_idx = find(nr_timeline > 0 - ped_stop & ...
                            nr_timeline < 0 - ped_start );

            stim_idx = find(nr_timeline > 0 + ped_start & ...
                            nr_timeline < 0 + ped_stop );
            

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
            fr_sus_mod_stats(nr).p_sign = p_sign;
            fr_sus_mod_stats(nr).sign_stats = stats;
            fr_sus_mod_stats(nr).mod = mod;

        end
        region_data.(f_region).(f_stim).fr_sus_mod_stats = fr_sus_mod_stats;
    end
end

%% Calculate single cell DBS-PLV as well as shuffled distribution
% This is used to determine between entrained and non-entrained neurons
% Creates 'plv_mod_stats' field

num_iter = 500; %change this back to 500
wind_dist = 1000/1000; %ms

% Determine whether or not to show figures
show_figures = 0;

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
        plv_freq = str2num(erase(f_stim, 'f_'));
        avg_Fs = mean(popul_data.framerate, 'omitnan');
 
        
            disp('Fig created');

        % Check if plv_mod_stats already exists
        %try
        %    plv_mod_stats = popul_data.plv_mod_stats;
        %    continue;
        %catch ME
            plv_mod_stats = struct;
        %end

        if show_figures == 1
            figure('Position', [0, 0, 800, 1800]);
            tiledlayout(length(popul_data.neuron_name), 2, 'TileSpacing', 'none', 'Padding', 'none');
        end

        for nr=1:length(popul_data.neuron_name)

            nr_timeline = popul_data.trace_timestamps(:, nr);
            base_idx = find(nr_timeline > 0 - wind_dist & nr_timeline < 0);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 0 + wind_dist);

            % Get the vector of all phases for all trials at the stimulation frequency
            filt_trial = @(trial) (angle(Multi_func.filt_single_freq(trial, plv_freq, avg_Fs)));
            
            % Function that waits for a function and matrix and then returns a function waiting for a column value
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            % Aply the filtering function with the whole trial matrix, so all is left is waiting for a column number
            partial_apply = applyFunToColi(filt_trial, popul_data.all_trial_rawVm{nr});
            % Iterate through each column trial and concatenate all o fthe results
            vm_phases = arrayfun(partial_apply, [1:size(popul_data.all_trial_rawVm{nr}, 2)]' , 'UniformOutput', false); 

            vm_phases = cat(3, vm_phases{:});

            % Calculate stimulation raster with the index of the whole trace
            nr_stim_time = reshape(popul_data.stim_timestamps(:, nr), 1, []);
            nr_frame_time = popul_data.trace_timestamps(:, nr);
            diffs = abs(nr_stim_time - nr_frame_time);
            [~, stim_idx_i] = min(diffs, [], 1);
            
            dbs_raster = zeros(size(nr_frame_time));
            dbs_raster(stim_idx_i) = 1;
            dbs_rasters = repmat(dbs_raster, 1, size(vm_phases, 3));
                  
            % DEBUG ensuring that the DBS pulses are properly shown (Its COrrect!)
            %nexttile;
            %plot(dbs_rasters + [1:size(dbs_rasters, 2)]);

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, dbs_rasters, 0, 10);             
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
            
            if show_figures == 1
                % Plot the trial-averaged Vm
                %nexttile;
                %plot(nr_timeline, popul_data.neuron_RawVm(:, nr));
        
                % plot the shuffled distribution
                nexttile;
                histogram(shuf_plv_adj, 1000);
                hold on;
                xline(high_prc, '-b');
                hold on;
                Multi_func.set_default_axis(gca);
            end

            % Check if significant
            if obs_PLV2 > high_prc
                if show_figures == 1
                    xline(obs_PLV2, '-g');
                end
                plv_mod_stats(nr).mod = 1;
            else
                if show_figures == 1
                    xline(obs_PLV2, '-r');
                end
                plv_mod_stats(nr).mod = -1;
            end

            % DEBUG plot the last randomized stim pulses
            %nexttile;
            %plot(dbs_raster);
            %hold on;
            %plot(stim_rasters + [1:size(stim_rasters, 2)])
    
    
            if show_figures == 1
                %DEBUG Plot the raw trials with spike phase information
                nexttile;
                all_trs = popul_data.all_trial_rawVm{nr};
                all_trs = Multi_func.raw_filt(all_trs, plv_freq, avg_Fs);
                norm_trs = Multi_func.norm_signals(all_trs);
                plot(norm_trs + repmat(1:size(norm_trs, 2), size(norm_trs, 1), 1));
                hold on;
                xline(find(dbs_raster == 1));
                 %plot(find(dbs_raster == 1), 1.2.*size(norm_trs, 2).*ones(sum(dbs_raster), 1), '|');
                xlim([700 800]);

                sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
            end

            % Store the PLV data into the neuron structure
            plv_mod_stats(nr).obs_PLV = obs_PLV;
            plv_mod_stats(nr).obs_PLV2 = obs_PLV2;
            plv_mod_stats(nr).shuf_PLV = sh_PLV;
            plv_mod_stats(nr).shuf_PLV2 = sh_PLV2;

        end %neuron loop
        % Add plv mod stats to structure
        region_data.(f_region).(f_stim).plv_mod_stats = plv_mod_stats;
        
        if show_figures == 1
            % Save the figures
            saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_stimPLV_shuf.png']);
            saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_stimPLV_shuf.pdf']);
        end
    end % Stim freq loop
end % Region loop

%% Determine neurons by Vm modulation, firing rate modulation, and entrainment
% Creates the 'mod_matrix'
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
    
        % Need to convert the non-entrained label of -1 to 0. 
        %Easier for intersection purposes
        nr_entr_mod = [popul_data.plv_mod_stats.mod]';
        nr_entr_mod(nr_entr_mod == -1) = 0;

        % Create matrix for all neuron modulation feature
        % I am taking the absolute value here to combine activated and suppressed
        mod_matrix = [nr_entr_mod, ...
                abs([popul_data.Vm_trans_mod_stats.mod])', ...
                abs([popul_data.Vm_sus_mod_stats.mod])', ...
                abs([popul_data.fr_trans_mod_stats.mod])', ...
                abs([popul_data.fr_sus_mod_stats.mod])', ...
                ];
        
        % Save modulation to region_data structure
        region_data.(f_region).(f_stim).mod_matrix = mod_matrix;

        % Save matrix to csv file
        headers = {['Ent ' num2str(100*sum(mod_matrix(:, 1))/size(mod_matrix, 1))], ...
            ['Vm trans ' num2str(100*sum(mod_matrix(:, 2))/size(mod_matrix, 1))], ...
            ['Vm sus ' num2str(100*sum(mod_matrix(:, 3))/size(mod_matrix, 1))], ...
            ['FR trans ' num2str(100*sum(mod_matrix(:, 4))/size(mod_matrix, 1))], ...
            ['FR sus ' num2str(100*sum(mod_matrix(:, 5))/size(mod_matrix, 1))]};

        T = array2table(mod_matrix, 'VariableNames', headers);
        writetable(T, [figure_path 'Neuronwise' f 'mod_matrix_' f_region '_' f_stim '.csv'], ...
            'WriteRowNames', true);

        % Find neurons that do not have any modulatory property
        disp([f_region ' ' f_stim]);
        non_mod = find(sum(mod_matrix, 2) == 0)

        % Find neurons that were sustained activated, but not in the transient
        sus_mod = find(mod_matrix(:, 3) == 1);
        trans_non_mod = find(mod_matrix(:, 2) == 0);
        disp('Sus only');
        sus_only = intersect(sus_mod, trans_non_mod)

        if isempty(non_mod)
            continue
        end

        % Display neuron trials of neurons without modulatory effect
        figure;
        tiledlayout(length(non_mod), 1, 'TileSpacing', 'none', 'Padding', 'none');
        for nr=non_mod'
            nexttile;
            plot(Multi_func.norm_signals(popul_data.all_trial_rawVm{nr}) ...
                + [1:size(popul_data.all_trial_rawVm{nr}, 2)]);
            hold on;
            stim = popul_data.stim_timestamps(:, nr);
            frame_time = popul_data.trace_timestamps(:, nr);
            stim_onset = find(frame_time < stim(1));
            stim_onset = stim_onset(end);

            stim_offset = find(frame_time > stim(end));
            stim_offset = stim_offset(1);

            xline([stim_onset, stim_offset]);
            title(['Non-Feature ' popul_data.neuron_name{nr}], 'Interpreter', 'none');
        end
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        legend(['Num Neuron: ' num2str(length(popul_data.neuron_name))]);

    end
end

%% Save the number of neurons per modulation

% Note: to get percentages, I have to copy the top part of this csv into the top of formula.xlsx file

% a flag here for whether or not to include non-modulated neurons or not
% If 1, then keep the non-modulated neurons in the counts calculations
% if 0, then remove the non-modulated neurons from the count calculations
keep_nonmod_nrs = 0;

stats_t = table();
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'

        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
       
        %--- Count neurons that have any modulation vs no modulation
        mod_nr = find(sum(popul_data.mod_matrix, 2) > 0);
        non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);

        % Grab the neuron number of each Vm trans positive modulation
        vm_trans_act_nr = find([popul_data.Vm_trans_mod_stats.mod] > 0);
        vm_trans_non_nr = find([popul_data.Vm_trans_mod_stats.mod] == 0);
        vm_trans_sup_nr = find([popul_data.Vm_trans_mod_stats.mod] < 0);

        % Grab the firing rate modulated trans data
        fr_trans_act_nr = find([popul_data.fr_trans_mod_stats.mod] > 0);
        fr_trans_non_nr = find([popul_data.fr_trans_mod_stats.mod] == 0);
        fr_trans_sup_nr = find([popul_data.fr_trans_mod_stats.mod] < 0);

        %----- Modulated data for sustained period
        % Grab the neuron number of each Vm sus positive modulation
        vm_sus_act_nr = find([popul_data.Vm_sus_mod_stats.mod] > 0);
        vm_sus_non_nr = find([popul_data.Vm_sus_mod_stats.mod] == 0);
        vm_sus_sup_nr = find([popul_data.Vm_sus_mod_stats.mod] < 0);

        % Grab the firing rate modulated sus data
        fr_sus_act_nr = find([popul_data.fr_sus_mod_stats.mod] > 0);
        fr_sus_non_nr = find([popul_data.fr_sus_mod_stats.mod] == 0);
        fr_sus_sup_nr = find([popul_data.fr_sus_mod_stats.mod] < 0);

        % Count the unique neurons that are Vm or firing rate modulated,
        % regardless of transient or sustained
        vm_act_tot_nrs = union(vm_trans_act_nr, vm_sus_act_nr);
        vm_sup_tot_nrs = union(vm_trans_sup_nr, vm_sus_sup_nr);
        vm_non_tot_nrs = intersect(vm_trans_non_nr, vm_sus_non_nr);

        fr_act_tot_nrs = union(fr_trans_act_nr, fr_sus_act_nr);
        fr_sup_tot_nrs = union(fr_trans_sup_nr, fr_sus_sup_nr);    
        fr_non_tot_nrs = intersect(fr_trans_non_nr, fr_sus_non_nr);

        % Count the total neurons that are Vm modulated for the given period
        vm_trans_tot_nrs = union(vm_trans_act_nr, vm_trans_sup_nr);
        vm_sus_tot_nrs = union(vm_sus_act_nr, vm_sus_sup_nr);

        %--- Determine the entrainment PLV mod stuff
        etrain_nr = find([popul_data.plv_mod_stats.mod] > 0);
        non_etrain_nr = find([popul_data.plv_mod_stats.mod] < 0);

        % Count the total number of neurons
        total_nrs = 1:length([popul_data.plv_mod_stats.mod]);

        % Determine whether to keep the non-modulated neurons
        if keep_nonmod_nrs == 0
            vm_trans_non_nr = setdiff(vm_trans_non_nr, non_mod_nr);
            fr_trans_non_nr = setdiff(fr_trans_non_nr, non_mod_nr);
            vm_sus_non_nr = setdiff(vm_sus_non_nr, non_mod_nr);
            fr_sus_non_nr = setdiff(fr_sus_non_nr, non_mod_nr);
            vm_non_tot_nrs = setdiff(vm_non_tot_nrs, non_mod_nr);
            fr_non_tot_nrs = setdiff(fr_non_tot_nrs, non_mod_nr);
            non_etrain_nr = setdiff(non_etrain_nr, non_mod_nr);
            total_nrs = setdiff(total_nrs, non_mod_nr);
        end

        % Add all of the counts to table
        stats_t([f_region f_stim], 'Vm Trans Act') = {length(vm_trans_act_nr)};
        stats_t([f_region f_stim], 'Vm sus Act') = {length(vm_sus_act_nr)};

        stats_t([f_region f_stim], 'Vm Trans Sup') = {length(vm_trans_sup_nr)};
        stats_t([f_region f_stim], 'Vm sus Sup') = {length(vm_sus_sup_nr)};

        stats_t([f_region f_stim], 'FR Trans Act') = {length(fr_trans_act_nr)};
        stats_t([f_region f_stim], 'FR sus Act') = {length(fr_sus_act_nr)};
        
        stats_t([f_region f_stim], 'PLV Etrain') = {length(etrain_nr)};
        stats_t([f_region f_stim], 'PLV Non Etrain') = {length(non_etrain_nr)};

        stats_t([f_region f_stim], 'Modulated Neurons') = {length(mod_nr)};

        %--- Add total number of neurons
        stats_t([f_region f_stim], 'Total Neurons') = {length(total_nrs)};

        % Determine how many transiently modulated neurons are sustainly modulated
        disp(['---' f_region f_stim '---']);
        disp(['Different temporally modulated nrs ' ...
            num2str(setdiff(vm_trans_tot_nrs, vm_sus_tot_nrs))]);
        vm_trans_tot_nrs
        vm_sus_tot_nrs
        fprintf("\n\n");

        % %--- Add total number of neurons
        %stats_t([f_region f_stim], 'Total Neurons') = {length(total_nrs)};

        % Total Vm/Firing rate modulated regardless of transient or sustained
        stats_t([f_region f_stim], 'Vm Total Activated') = {length(vm_act_tot_nrs)}; 
        stats_t([f_region f_stim], 'Vm Total Supressed') = {length(vm_sup_tot_nrs)};
        stats_t([f_region f_stim], 'Vm Total Unchanged') = {length(vm_non_tot_nrs)};

        stats_t([f_region f_stim], 'FR Total Activated') = {length(fr_act_tot_nrs)}; 
        stats_t([f_region f_stim], 'FR Total Supressed') = {length(fr_sup_tot_nrs)}; 
        stats_t([f_region f_stim], 'FR Total Unchanged') = {length(fr_non_tot_nrs)};

        % Total Vm modulated for each transient or sustained
        stats_t([f_region f_stim], 'Vm Total Transient') = {length(vm_trans_tot_nrs)};
        stats_t([f_region f_stim], 'Vm Total Sustained') = {length(vm_sus_tot_nrs)};

        % Add all of the counts to table (-- Old format)
        % stats_t([f_region f_stim], 'Vm Trans Act') = {length(vm_trans_act_nr)};
        % stats_t([f_region f_stim], 'Vm Trans Sup') = {length(vm_trans_sup_nr)};
        % stats_t([f_region f_stim], 'Vm Trans Non') = {length(vm_trans_non_nr)};
    
        % stats_t([f_region f_stim], 'FR Trans Act') = {length(fr_trans_act_nr)};
        % stats_t([f_region f_stim], 'FR Trans Sup') = {length(fr_trans_sup_nr)};
        % stats_t([f_region f_stim], 'FR Trans Non') = {length(fr_trans_non_nr)};

        % stats_t([f_region f_stim], 'Vm sus Act') = {length(vm_sus_act_nr)};
        % stats_t([f_region f_stim], 'Vm sus Sup') = {length(vm_sus_sup_nr)};
        % stats_t([f_region f_stim], 'Vm sus Non') = {length(vm_sus_non_nr)};
    
        % stats_t([f_region f_stim], 'FR sus Act') = {length(fr_sus_act_nr)};
        % stats_t([f_region f_stim], 'FR sus Sup') = {length(fr_sus_sup_nr)};
        % stats_t([f_region f_stim], 'FR sus Non') = {length(fr_sus_non_nr)};
        % 
        % stats_t([f_region f_stim], 'PLV Etrain') = {length(etrain_nr)};
        % stats_t([f_region f_stim], 'PLV non-etrain') = {length(non_etrain_nr)};

        % stats_t([f_region f_stim], 'Modulated Neurons') = {length(mod_nr)};
        % stats_t([f_region f_stim], 'Non-Modulated Neurons') = {length(non_mod_nr)};

        % %--- Add total number of neurons
        % stats_t([f_region f_stim], 'Total Neurons') = {length(total_nrs)};

        % Loop and create percentages
        %for row = stats_t.Properties.RowNames'
        %    stats_t{row{1}, :} = stats_t{row{1}, :}*100 / stats_t{row{1}, 'Total Neurons'};
        %end
    end
end
disp(stats_t);
writetable(stats_t, [figure_path 'Neuronwise' f 'modulation_stats_keepNonMod' num2str(keep_nonmod_nrs) '.csv'], 'WriteRowNames', true);

%% Plot the Vm modulation heatmap and modulated averages based on the transient OR sustained phase
display_names = 0;
remove_nonmod_nrs = 1;
stim_wind_sort = 100/1000;

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Get the neuron idx of neurons that were modulated
        %Just for neurons from the transient phase
        %nr_act = find([popul_data.Vm_trans_mod_stats.mod] > 0);
        %nr_sup = find([popul_data.Vm_trans_mod_stats.mod] < 0);
        %nr_non = find([popul_data.Vm_trans_mod_stats.mod] == 0);
 
        nr_trans_act = find([popul_data.Vm_trans_mod_stats.mod] > 0); 
        nr_sus_act = find([popul_data.Vm_sus_mod_stats.mod] > 0);

        nr_trans_sup = find([popul_data.Vm_trans_mod_stats.mod] < 0); 
        nr_sus_sup = find([popul_data.Vm_sus_mod_stats.mod] < 0);

        nr_non = find([popul_data.Vm_trans_mod_stats.mod] == 0 ...
            & [popul_data.Vm_sus_mod_stats.mod] == 0);
        
        % These set operations will include the sustained modulation and neurons
        % that have the same transient modulation AND do not have the sustained opposite modulation
        % This way how ever the neuron is modulated during the sustained period will take priority for plotting
        % than its transient modulation

        % Consider activated neurons that have sustained activation, or transient activation that do not have sustained suppression
        nr_act = nr_sus_act;
        nr_act = unique([nr_act, setdiff(nr_trans_act, nr_sus_sup)]);
        
        % Consider suppressed neurons that have sustained suppression, or transient suppression without sustained activation
        nr_sup = nr_sus_sup;
        nr_sup = unique([nr_sup, setdiff(nr_trans_sup, nr_sus_act)]);

        % Check to ensure there is no overlap between nr_act and nr_sup
        if sum(ismember(nr_act, nr_sup)) > 0
            error('Overlapping activation and suppression');
        end

        % Filter out the neurons that were non-modulated at all
        if remove_nonmod_nrs == 1
            non_mod_nr = find(sum(popul_data.mod_matrix, 2) == 0);
            nr_non = nr_non(~ismember(nr_non, non_mod_nr));
        end

        %DEBUG
        %disp([f_region ' ' f_stim]);
        %pause();

        % Grab idxs for different parts of the signal
        base_idx = find(mean(popul_data.trace_timestamps, 2) < 0);

        stim_idx = find(mean(popul_data.trace_timestamps, 2) > 0 &...
                    mean(popul_data.trace_timestamps, 2) < 0 + stim_wind_sort);
        
        % Perform spike-amplitude normalization and baseline subtraction
        norm_vms = popul_data.neuron_RawVm./popul_data.neuron_spike_amp;
        pop_base = mean(norm_vms(base_idx, :), 1);
        norm_vms = norm_vms - pop_base;

        % Sorts groups from the RawVm
        [~, act_i] = sort(mean(norm_vms(stim_idx, nr_act), 1), 'descend');
        [~, sup_i] = sort(mean(norm_vms(stim_idx, nr_sup), 1), 'descend');
        [~, non_i] = sort(mean(norm_vms(stim_idx, nr_non), 1), 'descend');
        act_i = nr_act(act_i);
        sup_i = nr_sup(sup_i);
        non_i = nr_non(non_i);

        % Save the trial-averaged neuron Vm
        pop_act_vm = norm_vms(:, act_i);
        pop_non_vm = norm_vms(:, non_i);
        pop_sup_vm = norm_vms(:, sup_i);

        plot_act_vm = pop_act_vm;
        plot_non_vm = pop_non_vm;
        plot_sup_vm = pop_sup_vm;       
        
        % Zscore all of the trial-averaged neuron Vm
        %plot_act_vm = zscore(plot_act_vm, [], 1);
        %plot_non_vm = zscore(plot_non_vm, [], 1);
        %plot_sup_vm = zscore(plot_sup_vm, [], 1);

        %Do I need to sort again after zscoring the Vm??

        % Baseline subtract the vm data
        %plot_act_vm = plot_act_vm - mean(plot_act_vm(base_idx, :), 1);
        %plot_non_vm = plot_non_vm - mean(plot_non_vm(base_idx, :), 1);
        %plot_sup_vm = plot_sup_vm - mean(plot_sup_vm(base_idx, :), 1);

        % Zscore the trial-averaged neuron Vm using baseline mean and standard deviation
        plot_act_vm = (plot_act_vm - mean(plot_act_vm(base_idx, :), 1))./std(plot_act_vm(base_idx, :), 0, 1);
        plot_non_vm = (plot_non_vm - mean(plot_non_vm(base_idx, :), 1))./std(plot_non_vm(base_idx, :), 0, 1);
        plot_sup_vm = (plot_sup_vm - mean(plot_sup_vm(base_idx, :), 1))./std(plot_sup_vm(base_idx, :), 0, 1);
        
        % Just divide by the max value of the Vm
        % This DOES NOT work well
        %plot_act_vm = plot_act_vm./max(plot_act_vm, 1);
        %plot_non_vm = plot_non_vm./max(plot_non_vm, 1);
        %plot_sup_vm = plot_sup_vm./max(plot_sup_vm, 1);

        % Construct Vm heatmap
        vm_heatmap = plot_act_vm;
        vm_heatmap = horzcat_pad(vm_heatmap, plot_non_vm);
        vm_heatmap = horzcat_pad(vm_heatmap, plot_sup_vm)';

        % Calculate the population average of each groups
        % This is used for the average traces right below the heatmaps
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
        %surface(timeline, 1:size(vm_heatmap, 1), vm_heatmap, ...
        %    'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        
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
                text(2.2, i, [num2str(popul_data.Vm_trans_mod_stats(nr).mod) ' '...
                            popul_data.neuron_name{nr}], ...
                            'Interpreter', 'none');
                i = i+1;

                %TODO need to figure out how to print values here for DEBUGGING why the plots get different, honestly
                % First printing the different neuron activation and suppression may be the way to go
                if  strcmp(f_region, 'r_V1') == 1 && strcmp(f_stim, 'f_140') == 1
                    disp(['Sustained neurons list ' num2str(sup_i)]);
                    disp(find([popul_data.Vm_sus_mod_stats.mod] == -1) );
                    %error('Pause for V1 140 Hz');
                end

            end
            xlim([-1, 4.5]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end

        %TODO maybe plot the DBS bar on top??
        ylim([0.5 size(vm_heatmap, 1) + 0.5]);
        xlabel('Time from onset (S)');

        %colormap(red_blue_color_cmap);
        colormap(Multi_func.red_purple_blue_color);
        c = colorbar;
        c.Label.String = 'Vm';
        max_abs = max(abs(c.Limits));
        caxis([-max_abs, max_abs]);

        % Plot the population average for each group
        nexttile([1 1]);

        fill_h = fill([timeline; flip(timeline)], ...
                [non_Vm_avg + non_Vm_sem; flipud(non_Vm_avg - non_Vm_sem)], ...
                Multi_func.non_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);

        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, non_Vm_avg, 'DisplayName', 'Non-Responsive Vm', 'Color', Multi_func.non_color);
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
        if display_names == 1
            xlim([-.7 2.05]);
        else
            xlim([min(timeline) - .01, max(timeline)]);
        end

        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_Vm_mod_plots.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_' num2str(wind_dist) '_Vm_mod_plots.pdf']);
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


        colormap(red_blue_color_cmap);
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

%% 
% Debugging script to determine the sustained mods for each region
T = array2table([region_data.r_M1.f_140.fr_sus_mod_stats.mod]);
writetable(T, 'M1_140');


T = array2table([region_data.r_M1.f_40.fr_sus_mod_stats.mod]);
writetable(T, 'M1_40');


T = array2table([region_data.r_V1.f_140.fr_sus_mod_stats.mod]);
writetable(T, 'V1_140');


T = array2table([region_data.r_V1.f_40.fr_sus_mod_stats.mod]);
writetable(T, 'V1_40');



%% Calculate significant stim-Vm PLV for both halves of the stimulation period
num_iter = 500;
end_point = 1000/1000; %ms
start_point = 0/1000; %ms
mid_point = 500/1000 %ms

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
        
        % Check if plv_mod_stats already exists
        try
            plv_mod_stats = popul_data.plv_mod_stats;
        catch ME
            plv_mod_stats = struct;
        end

        figure('Position', [0, 0, 800, 1000]);
        tiledlayout(length(popul_data.neuron_name), 2, 'TileSpacing', 'none', 'Padding', 'none');

        for nr=1:length(popul_data.neuron_name)
            nr_timeline = popul_data.trace_timestamps(:, nr);
            stim_f_idx = find(nr_timeline > start_point & nr_timeline < mid_point);
            stim_l_idx = find(nr_timeline > mid_point & nr_timeline < end_point);

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
            
            % Calculate PLVs for the first half of stimulation
            %TODO I need to double check this is doing the proper thing
            stim_rasters = zeros(size(dbs_rasters));
            stim_f_idx = repmat(stim_f_idx, 1, size(vm_phases, 3));
            col_ind_mat = repmat(1:size(stim_f_idx, 2), size(stim_f_idx, 1), 1);
            stim_rasters(sub2ind(size(stim_rasters), ...
                        stim_f_idx(:), ...
                        col_ind_mat(:))) = ...
            dbs_rasters(sub2ind(size(stim_rasters), ...
                        stim_f_idx(:), ...
                        col_ind_mat(:)));
            
            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, stim_rasters, 0, 10);             
            first_obs_PLV = PLV;
            first_obs_PLV2 = PLV2;
            first_obs_norm_vecs= norm_vecs;

            % Calculate PLVs for the second half of stimulation
            stim_rasters = zeros(size(dbs_rasters));
            stim_l_idx = repmat(stim_l_idx, 1, size(vm_phases, 3));
            col_ind_mat = repmat(1:size(stim_l_idx, 2), size(stim_l_idx, 1), 1);
            stim_rasters(sub2ind(size(stim_rasters), ...
                        stim_l_idx(:), ...
                        col_ind_mat(:))) = ...
            dbs_rasters(sub2ind(size(stim_rasters), ...
                        stim_l_idx(:), ...
                        col_ind_mat(:)));

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, stim_rasters, 0, 10);             
            last_obs_PLV = PLV;
            last_obs_PLV2 = PLV2;
            last_obs_norm_vecs= norm_vecs;

            % Find the upper limit for randomization
            filt_time = nr_frame_time(mid_point > nr_frame_time);
            % Compute absolute differences
            diffs = abs(filt_time - mid_point);
            % Find index of minimum difference in filtered_array
            [minDiff, minIndex] = min(diffs);
            % Find index of corresponding element in original array
            last_index = find(nr_frame_time == filt_time(minIndex), 1);

            % Shuffled PLV data
            shuf_plv = [];
            shuf_plv_adj = [];

            % Reset the starting stim_idx for each trial and keep the spacing between indices the same
            stim_idx_i = stim_idx_i(ceil(length(stim_idx_i)/2):end);

            % Shuffle the start of the dbs timepoints and recalculate the PLV
            % Need to randomize where the DBS points are along the rasters
            for i=1:num_iter
                % Randomly select stim onset time
                onset_rand = randi(last_index, 1, size(vm_phases, 3));
                
                
                rand_stim_idx = repmat(stim_idx_i' - stim_idx_i(1), 1, size(vm_phases, 3)) + onset_rand;
                col_ind_mat = repmat(1:size(rand_stim_idx, 2), size(rand_stim_idx, 1), 1);
                stim_rasters = zeros(size(dbs_rasters));
                stim_rasters(sub2ind(size(stim_rasters), ...
                                     rand_stim_idx(:), ...
                                     col_ind_mat(:))) = 1;
            
                [sh_PLV, sh_PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, stim_rasters, 0, 10);            
                shuf_plv(end + 1) = sh_PLV;
                shuf_plv_adj(end + 1) = sh_PLV2;
            end
        
            % Calculate the 95th percentile
            high_prc = prctile(shuf_plv_adj, 95);

            % Plot the trial-averaged Vm
            %nexttile;
            %plot(nr_timeline, popul_data.neuron_RawVm(:, nr));

            % plot the shuffled distribution
            nexttile;
            histogram(shuf_plv_adj, 1000);
            hold on;
            xline(high_prc, '-b');
            hold on;
            Multi_func.set_default_axis(gca);
            
            % Check if significant for the first half of stimulation
            if first_obs_PLV2 > high_prc
                xline(first_obs_PLV2, '-g');
                plv_mod_stats(nr).first_mod = 1;
            else
                xline(first_obs_PLV2, '-r');
                plv_mod_stats(nr).first_mod = -1;
            end
        
            % Check if significant for the second half of stimulation
            if last_obs_PLV2 > high_prc
                xline(last_obs_PLV2, '--g');
                plv_mod_stats(nr).last_mod = 1;
            else
                xline(last_obs_PLV2, '--r');
                plv_mod_stats(nr).last_mod = -1;
            end

            % DEBUG plotting randomized stim within vectors
            nexttile;
            plot(dbs_raster);
            hold on;
            plot(stim_rasters + [1:size(stim_rasters, 2)]); % The original DBS pulses
        
            %DEBUG plotting the power spectra for each neuron
            %nexttile;

            %get_base_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp < stim_tmstp(1));
            %get_stim_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp >= stim_tmstp(1) & ...
            %                                           tr_tmstmp <= stim_tmstp(end));

            %calc_time_pow = @(trial_spec, time_idxs) mean(trial_spec(:, time_idxs, :), 2);
            %
            %calc_trial_spec = @(trial_spec, base_pow, stim_pow) (trial_spec - base_pow)./(base_pow + stim_pow);
            %
            %calc_nr_spec = @(trial_spec, tr_tmstmp, stim_tmstp) mean(calc_trial_spec(trial_spec, ...
            %        calc_time_pow(trial_spec, get_base_idxs(tr_tmstmp, stim_tmstp)),  ...
            %        calc_time_pow(trial_spec, get_stim_idxs(tr_tmstmp, stim_tmstp))), 3, 'omitnan');
            %
            %norm_spec = calc_nr_spec(popul_data.all_trial_power_spec{nr}, ...
            %                        popul_data.trace_timestamps(:, nr), ...
            %                        popul_data.stim_timestamps(:, nr));

            % %avg_spec = mean(popul_data.all_trial_power_spec{nr}, 3);
            %avg_freqs = mean(popul_data.all_trial_spec_freq{nr}, 3);
            %surface(nr_timeline, avg_freqs, norm_spec, ...
            %'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            %xlim([-0.75 2.02]);
            %hold on;
            %xline([0 0.5]);

            sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');

            % Store the PLV data into the neuron structure
            plv_mod_stats(nr).first_obs_PLV = first_obs_PLV;
            plv_mod_stats(nr).first_obs_PLV2 = first_obs_PLV2;
            plv_mod_stats(nr).last_obs_PLV = last_obs_PLV;
            plv_mod_stats(nr).last_obs_PLV2 = last_obs_PLV2;
            plv_mod_stats(nr).half_shuf_PLV = sh_PLV;
            plv_mod_stats(nr).half_shuf_PLV2 = sh_PLV2;

        end %neuron loop
        
        % Add plv mod stats to structure
        region_data.(f_region).(f_stim).plv_mod_stats = plv_mod_stats;
            
        % Save the figures
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_stimPLV_half_shuf.png']);
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_stimPLV_half_shuf.pdf']);

    end % Stim freq loop
end % Region loop

%% Plot single cell pulse triggered average across transient and sustained based on stim-Vm PLV entrainement
% Need to run All the options for the variable 'nr_pop' 
%Vm pulse-triggered for both 'etrain' and 'non_entr' in 'pulse_small_res.m' for plotting the heatmaps
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
        

        % Set up figure for plotting
        % Each column is for: 'all', 'trans', and 'sus' 
        figure('Position', [0, 0, 1500, 1000]);
        %tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.9, 6]);
        tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 11.5, 6]);
        %tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        N_col = 3;

        % Loop through all, sus, and transient pulses
        tilenum = 1;

        % Reset the sorting for the entrained and non-entrained neurons
        etrain_i = [];
        non_etrain_i = [];
        for f_ped = {'all', 'trans', 'sus'}
            f_ped = f_ped{1};
            
            % Create heatmap for each neuron
            %cur_pulse_vms = popul_data.etrain.([f_ped '_pulse_trig_Vm']);
            num_stim_pulses = str2num(erase(f_stim, 'f_'));
            num_pulse_wind = size(popul_data.etrain.all_pulse_trig_Vm, 1); 
            base_sub = @(mat) mat - mat(extra_trace + 1, :);
            
            etrain_avg_pulse_mat = cellfun(@(x) mean(base_sub(reshape(x', num_pulse_wind, [])), 2, 'omitnan'), ...
                popul_data.etrain.([f_ped '_pulse_vm_by_neuron']), 'UniformOutput', false)';
            etrain_avg_pulse_mat = cat(2, etrain_avg_pulse_mat{:})';

            % Non-entrained pulse window
            non_etrain_avg_pulse_mat = cellfun(@(x) mean(base_sub(reshape(x', num_pulse_wind, [])), 2, 'omitnan'), ...
                popul_data.non_entr.([f_ped '_pulse_vm_by_neuron']), 'UniformOutput', false)';
            non_etrain_avg_pulse_mat = cat(2, non_etrain_avg_pulse_mat{:})';

            % Sort by the max Vm and concatenate the Vms
            if isempty(etrain_i) && isempty(non_etrain_i)
                [~, etrain_i] = sort(max(etrain_avg_pulse_mat, [], 2), 'descend');
                [~, non_etrain_i] = sort(max(non_etrain_avg_pulse_mat, [], 2), 'descend');
            end
            pulse_heatmap = cat(1, etrain_avg_pulse_mat(etrain_i, :), ...
                                non_etrain_avg_pulse_mat(non_etrain_i, :) );

            %-- Plot the heatmap with the pulse-triggered averages--
            nexttile(tilenum);
            
            timeline = ([0:size(pulse_heatmap, 2) - 1]- extra_trace)*1000./avg_Fs;
            imagesc('XData', timeline, 'YData', 1:size(pulse_heatmap, 1), 'CData', pulse_heatmap);
            hold on;
            yline(size(etrain_avg_pulse_mat, 1) + 0.5);
            hold on;
            Multi_func.set_default_axis(gca);
            xlim([min(timeline) - .01, max(timeline)]);
            
            ylim([0.5 size(pulse_heatmap, 1) + 0.5]);
            xlabel('Time from onset (ms)');
            
            if strcmp(f_ped, 'sus') == 1
                colormap(Multi_func.red_purple_blue_color);
                c = colorbar;
                c.Label.String = 'Vm';
            end
            
            clim = caxis;
            caxis([-1 1]*range(clim)/2);
            
             title([f_stim ' ' f_ped], 'Interpreter', 'none');
            %-- End plotting the heatmap with the pulse-triggered averages--
            
            % -- Plotting the averages of the entrained and non-entrained population
            nexttile(tilenum + N_col);
            cur_etrain_pulse_vms = popul_data.etrain.([f_ped '_pulse_trig_Vm']);
            cur_etrain_avg_Vm = mean(cur_etrain_pulse_vms, 2, 'omitnan');
            cur_etrain_sem_Vm = std(cur_etrain_pulse_vms, 0, 2, 'omitnan')./sqrt(size(cur_etrain_pulse_vms, 2));
            
            fill_h = fill([timeline, flip(timeline)], [cur_etrain_avg_Vm + cur_etrain_sem_Vm; flipud(cur_etrain_avg_Vm - cur_etrain_sem_Vm)], [0.5 0.5 0.5]);
            if ~isempty(fill_h)
                Multi_func.set_fill_properties(fill_h);
            end
            hold on;
            plot(timeline, cur_etrain_avg_Vm, '-r', 'LineWidth', 1, 'DisplayName', 'Entrained');
            hold on;
            
            %TODO plot the red points on the actual line
            plot(timeline(sig_idx), max(cur_etrain_avg_Vm + cur_etrain_sem_Vm)*ones(size(sig_idx)), ...
                'b*', 'MarkerSize', 4)
            %plot(timeline(sig_idx), cur_etrain_avg_Vm(sig_idx), ...
            %    'b.', 'MarkerSize', 6)
            hold on;

            % Try to plot non-entrained neurons if it exists
            try
                cur_nonetrain_pulse_vms = popul_data.non_entr.([f_ped '_pulse_trig_Vm']);
                cur_nonetrain_avg_Vm = mean(cur_nonetrain_pulse_vms, 2, 'omitnan');
                cur_nonetrain_sem_Vm = std(cur_nonetrain_pulse_vms, 0, 2, 'omitnan')./sqrt(size(cur_nonetrain_pulse_vms, 2));
                
                fill_h = fill([timeline, flip(timeline)], [cur_nonetrain_avg_Vm + cur_nonetrain_sem_Vm; flipud(cur_nonetrain_avg_Vm - cur_nonetrain_sem_Vm)], [0.5 0.5 0.5]);
                if ~isempty(fill_h)
                    Multi_func.set_fill_properties(fill_h);
                end
                hold on;
                plot(timeline, cur_nonetrain_avg_Vm, '-k', 'LineWidth', 1, 'DisplayName', 'Non-Entrained');
                hold on;
            
            catch ME
                
            end
            
            % Plot the DBS stimulation time pulses
            nr_avg_pulse_width_time = popul_data.nr_avg_pulse_width_time;
            xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], 'Color', Multi_func.pulse_color, 'LineWidth', 2);
            hold on;

            % Increase timescale resolution
            %xlim([0 - .100, 0 + .100]);
            Multi_func.set_default_axis(gca);
            ylabel('Normalized Vm');
            xlabel('Time from pulse(ms)');

            % -- End plotting population averages

            % Increment to the next period pulse-triggered
            tilenum = tilenum + 1;
        end

        % Give the whole title
        sgtitle([f_region], 'Interpreter', 'none');

        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_entrain_pulse.png']);
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_entrain_pulse.pdf']);
    end
end

%% Plot different aspects of neuron's data based on stim-Vm PLV entrainment
% This is doing 'all_pulse' all stimulation stuff
% Note:To use 'stim_dbsvm_plvs_adj', need to run 'dbs_PLV.m' and the section on calculating dbs-Vm PLVs
% Note:To use ''

%plot_mode = 'pow'; % Use stimulation frequency power spectra for each neuron
%plot_mode = 'Vm'; % Use the trial averaged Vm
plot_mode = 'pulse'; % use Pulse-triggered average
freqs = Multi_func.entr_freqs;
extra_trace = 3;

% This is mostly for debugging purposes
display_names = 0;

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
                pop_etrain_vm = popul_data.neuron_RawVm(:, etrain_i);
                pop_non_vm = popul_data.neuron_RawVm(:, non_i);

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
        %tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 10, 20]);
        tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 3.25, 10]);

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

                colormap(Multi_func.red_purple_blue_color);
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

                colormap(Multi_func.red_purple_blue_color);
                c = colorbar;
                c.Label.String = 'Vm (z-scored)';
                %---------- End Vm heatmap plotting
            case 'pulse'
                nexttile;
                
                % Just take from first neuron
                nr_avg_pulse_width_time = mean(diff(popul_data.stim_timestamps(:, 1) ), 'omitnan');
                
                timeline = ([0:size(pulse_heatmap, 2) - 1]- extra_trace )*1000./avg_Fs;

                imagesc('XData', timeline, 'YData', 1:size(pulse_heatmap, 1), 'CData', pulse_heatmap);
                hold on;
                yline(0.5 + [length(etrain_i), length(etrain_i) + length(non_i)]);
                hold on;

                xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], ...
                    'Color', 'k', 'LineWidth', 2);
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

                colormap(Multi_func.red_purple_blue_color);
                c = colorbar;
                clim = caxis;
                caxis([-1 1]*range(clim)/2);
                c.Label.String = 'Vm';

        end

        % Plot the population average PLV for each group entrainment
        %nexttile;
        %
        %fill_h = fill([freqs'; flip(freqs)'], ...
        %        [non_plv_avg + non_plv_sem; flipud(non_plv_avg - non_plv_sem)], ...
        %        Multi_func.non_color, 'HandleVisibility', 'off');
        %Multi_func.set_fill_properties(fill_h);

        %fill_h.EdgeAlpha = 0;
        %hold on;
        %plot(freqs', non_plv_avg, 'DisplayName', 'Non-entrained', 'Color', Multi_func.non_color);

        %hold on;
        %fill_h = fill([freqs'; flip(freqs)'], ...
        %        [etrain_plv_avg + etrain_plv_sem; flipud(etrain_plv_avg - etrain_plv_sem)], ...
        %        Multi_func.act_color, 'HandleVisibility', 'off');
        %Multi_func.set_fill_properties(fill_h);
        %fill_h.EdgeAlpha = 0;
        %hold on;
        %plot(freqs', etrain_plv_avg, 'DisplayName', 'Entrained', 'Color', Multi_func.act_color);
        %hold on;

        %set(gca, 'XScale', 'log');
        %Multi_func.set_default_axis(gca);
        %legend();
        %ylabel('DBS-Vm PLV^2');
        %xlabel('Frequency (Hz)');
        % -- Done plotting the population average for each group entrainment

        % Plot the pulse-triggered average for both entrained and non-entrained
        nexttile;
        etrain_pulse_avg = mean(pop_etrain_pulse, 2);
        etrain_pulse_sem = std(pop_etrain_pulse, 0, 2)./sqrt(size(pop_etrain_pulse, 2));
        non_pulse_avg = mean(pop_non_pulse, 2);
        non_pulse_sem = std(pop_non_pulse, 0, 2)./sqrt(size(pop_non_pulse, 2));
 
        % Plot the Non-Entrained pulse-triggered average
        fill_h = fill([timeline, flip(timeline)], ...
            [non_pulse_avg + non_pulse_sem; flipud(non_pulse_avg - non_pulse_sem)], ...
            Multi_func.non_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, non_pulse_avg, 'DisplayName', 'Non-entrained', 'Color', Multi_func.non_color);
        hold on;

        % Plot the Entrained pulse-triggered average
        fill_h = fill([timeline, flip(timeline)], ...
            [etrain_pulse_avg + etrain_pulse_sem; flipud(etrain_pulse_avg - etrain_pulse_sem)], ...
            Multi_func.act_color, 'HandleVisibility', 'off');
        Multi_func.set_fill_properties(fill_h);
        fill_h.EdgeAlpha = 0;
        hold on;
        plot(timeline, etrain_pulse_avg, 'DisplayName', 'Entrained', 'Color', Multi_func.act_color);
        hold on;
        Multi_func.set_default_axis(gca);

        xlim([min(timeline) - .01, max(timeline)]);
        xline([0:nr_avg_pulse_width_time*1000:nr_avg_pulse_width_time*1000], ...
            'Color', Multi_func.pulse_color, 'LineWidth', 2, ...
            'HandleVisibility', 'off');
        hold on;
        legend();
        xlabel('Time from pulse onset (S)');
        ylabel('Normalized Vm');

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

        colormap(Multi_func.red_purple_blue_color);
        c = colorbar;
        caxis([-1 1]);
        xlim([0 200]);
        ylim([1 size(plv_heatmap, 1) ]);
        xlabel('Frequency (Hz)');
        set(gca, 'XScale', 'log');
        c.Label.String = 'DBS-Vm PLV';

        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        
        % Save figure stuff
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_' num2str(wind_dist) '_stimVm_plv_mod_plots_mode' plot_mode '.png']);
        saveas(gcf, [figure_path 'Neuronwise' f f_region '_' f_stim '_' num2str(wind_dist) '_stimVm_plv_mod_plots_mode' plot_mode '.pdf']);
    end
end

%% Calculate the cross-correlation of the stimulation period for all neurons

for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Get average framerate
        avg_Fs = mean(popul_data.framerate, 'omitnan');

        % Figure for each condition
        figure('Position', [0, 0, 800, 2000]);
        tiledlayout(length(popul_data.neuron_name), 2, 'TileSpacing', 'none', 'Padding', 'none');

        % Iterate through neurons
        for nr=1:length(popul_data.neuron_name)
            
            % Grab the stimulation idxs
            nr_timeline = popul_data.trace_timestamps(:, nr);
            stim_idx = find(nr_timeline > 0 & nr_timeline < 1);

            % De-mean and calculate the auto-correlation
            de_mean = @(col_idx) xcorr(popul_data.all_trial_SubVm{nr}(stim_idx, col_idx) ...
                                    - mean(popul_data.all_trial_SubVm{nr}(stim_idx, col_idx)));
            

            [corrs, lags] = arrayfun(de_mean, [1:size(popul_data.all_trial_SubVm{nr}, 2)], 'UniformOutput', false);
            
            % Calculate the average correlations and lags
            avg_corrs = mean(cat(2, corrs{:}), 2);
            avg_lags = mean(cat(1, lags{:})', 2);
            sem_corrs = std(cat(2, corrs{:}), [], 2)./sqrt(length(corrs));

            % Chop the data in half because it is symmetrical
            avg_corrs = avg_corrs(ceil(length(avg_corrs)/2):end);
            avg_lags = avg_lags(ceil(length(avg_lags))./2:end);
            sem_corrs = sem_corrs(ceil(length(sem_corrs))./2:end);

            % Calculate the fft of the average cross-correlation
            fft_y = fft(avg_corrs);
            sig_length = length(avg_corrs);
            freq_domain = (avg_Fs/sig_length)*(0:sig_length - 1);

            % Plot the auto-correlation
            nexttile;
            fill_h = fill([avg_lags; flip(avg_lags)], ...
                [avg_corrs + sem_corrs; flipud(avg_corrs - sem_corrs)], ...
                [0 0 1] , 'HandleVisibility', 'off');
            Multi_func.set_fill_properties(fill_h);
            hold on;
            plot(avg_lags, avg_corrs);
            xlim([-200 200]);

            % Plot the fourier-transform
            nexttile;
            
            % Check if the PLV mod is high or low
            if popul_data.plv_mod_stats(nr).mod > 0
                cur_color = 'g';
            else
                cur_color = 'r';
            end
            plot(freq_domain, abs(fft_y), 'Color', cur_color);
            xlim([5 200]);
        end

        % Title of condition
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
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
        
        disp([f_region ' ' f_stim ' ' num2str(length([popul_data.Vm_trans_mod_stats.mod]))]);

        % Display the modulation of Vm
        disp('Vm');
        disp(['Activated: ' num2str(length(find([popul_data.Vm_trans_mod_stats.mod] > 0 )))]);
        disp(['Non-modulated: ' num2str(length(find([popul_data.Vm_trans_mod_stats.mod] == 0 )))]);
        disp(['Suppression: ' num2str(length(find([popul_data.Vm_trans_mod_stats.mod] < 0 )))]);
        fprintf('\n');
        
        % Display the modulation of Firing Rate
        disp('Firing Rate');
        disp(['Activated: ' num2str(length(find([popul_data.fr_trans_mod_stats.mod] > 0 )))]);
        disp(['Non-modulated: ' num2str(length(find([popul_data.fr_trans_mod_stats.mod] == 0 )))]);
        disp(['Suppression: ' num2str(length(find([popul_data.fr_trans_mod_stats.mod] < 0 )))]);
        fprintf('\n');
    
        % Display the modulation of Firing Rate
        disp('PLV');
        disp(['Entrained: ' num2str(length(find([popul_data.plv_mod_stats.mod] > 0 )))]);
        disp(['Non-entrained: ' num2str(length(find([popul_data.plv_mod_stats.mod] < 0 )))]);
        fprintf('\n');

        % Display the number of neurons that have first and second half PLV entrainment
        disp('Stimulation Split PLV');
        disp(['First half entrained: ' num2str(length(find([popul_data.plv_mod_stats.first_mod] > 0 )))]);
        disp(['First half not entrained: ' num2str(length(find([popul_data.plv_mod_stats.first_mod] < 0 )))]);
        
        disp(['Second half entrained: ' num2str(length(find([popul_data.plv_mod_stats.last_mod] > 0 )))]);
        disp(['Second half not entrained: ' num2str(length(find([popul_data.plv_mod_stats.last_mod] < 0 )))]);

        fprintf('\n');
        fprintf('\n');
    end
end

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
        vm_act_nr = find([popul_data.Vm_trans_mod_stats.mod] > 0);
        vm_non_nr = find([popul_data.Vm_trans_mod_stats.mod] == 0);
        vm_sup_nr = find([popul_data.Vm_trans_mod_stats.mod] < 0);

        % Grab the firing rate modulated data
        fr_act_nr = find([popul_data.fr_trans_mod_stats.mod] > 0);
        fr_non_nr = find([popul_data.fr_trans_mod_stats.mod] == 0);
        fr_sup_nr = find([popul_data.fr_trans_mod_stats.mod] < 0);

        % Grab the PLV modulation
        etrain_nr = find([popul_data.plv_mod_stats.mod] > 0);
        non_etrain_nr = find([popul_data.plv_mod_stats.mod] < 0);

        sets = {vm_act_nr, vm_non_nr, vm_sup_nr,...
            fr_act_nr, fr_non_nr, fr_sup_nr, etrain_nr, non_etrain_nr};

        setLabels = {'Vm activated', 'Vm non-modul', 'Vm suppr', 'FR activated',...
            'FR non-modul', 'FR suppr', 'Entrained', 'Non-entrained'};

        DrawVennDiag(length(sets), sets, setLabels, colors);
    end
end

%% Entrained full stim and each halves in a Venn Diagram
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Define the sets from each entrainment method
        nr_full_etrain = [popul_data.plv_mod_stats.mod] > 0;
        nr_full_nonetrain = [popul_data.plv_mod_stats.mod] < 0;
        
        nr_first_etrain = [popul_data.plv_mod_stats.first_mod] > 0;
        nr_first_nonetrain = [popul_data.plv_mod_stats.first_mod] < 0;

        nr_last_etrain = [popul_data.plv_mod_stats.last_mod] > 0;
        nr_last_nonetrain = [popul_data.plv_mod_stats.last_mod] < 0;

        sets = {nr_full_etrain, nr_first_etrain, ...
            nr_last_etrain, nr_full_nonetrain, ...
            nr_first_nonetrain, nr_last_nonetrain };
        
        labels = {'Full Etrain', 'First Etrain', 'Last etrain', 'Full Non-etrain',...
            'First non-etrain', 'Last non-etrain'};

        % Visualize using a table
        T = table();

        % Loop through each entrainment method
        for var=1:round(length(sets)/2) % Only use the entraind population
            T = addvars(T, string(transpose(sets{var})), 'NewVariableNames', labels{var});     
        end
        
        % Loop through each column without using each name
        for i=1:width(T)
            col = T.Properties.VariableNames{i};
            T.(col)(T.(col) == 'true') = 'X';
            T.(col)(T.(col) == 'false') = ' ';
        end
        
        % Display the table
        writetable(T, [figure_path 'Neuronwise' f 'Entrainment_Stats' f ...
            f_region '_' f_stim '_entrain_table.csv']); 
        
        % DEBUG for specific conditions
        %if strcmp(f_region, 'r_M1') && strcmp(f_stim, 'f_40')
        %    T(:, 1)
        %    return;
        %end

    end
end

%% DEBUG venn digram

sets = {[1, 2, 3, 5, 6, 7], [3, 5, 5, 7, 8, 9]}

labels = {'Set 1', 'Set 2'};
circle_size = [length(sets{1}), length(sets{2})];

cols = {[1, 0, 0], [1, 0, 1]};

DrawVennDiag(length(circle_size), sets, labels, circle_size, cols);
