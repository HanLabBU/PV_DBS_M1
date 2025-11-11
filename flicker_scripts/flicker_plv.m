clc;
clear all;
close all;

%%
f= filesep;

addpath('..');

% This structure is taken from the python export file
interm_data_path = [f 'home' f 'pierfier' f 'Projects' f ...
      'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f];

savefig_path = [Multi_func.save_plot  'Flicker' f]; 

% Index for 'nr_mod' array the description of each  is storing
mod_str = {'inc', 'dec', 'unc'};

%res =  24; % Change to 24 for higher res, 10 for low res
voices_per_octaves = {'vpo_10', 'vpo_17', 'vpo_24'};

%% Calculate parameters for recording

total_frames = 2500;
total_time = 5.0258; % in Sec
soft_Fs = total_frames/total_time;

front_frame_drop = 14;

%%
save_data_file = [interm_data_path 'v1_flicker.mat'];
%save_data_file = [interm_data_path 'v1_flicker_aligned_dat.mat'];

data = load(save_data_file);

%% Calculate the phase-locking value between flicker and Vm across frequencies with trial aggregate
% This will show the flicker-only periods as well athe population level
freqs = Multi_func.entr_freqs;

for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;
    avg_Fs = 1/mean(diff(popul_data.interp_time));

    % Grab the flicker idxs
    flicker_base_idxs = find(timeline > 0 & timeline < 1);
    estim_idxs = find(timeline > 1 & timeline < 2);
    flicker_offset_idxs = find(timeline > 2 & timeline < 3);
    
    % Arrays to store PLVs
    flicker_base_phase_vectors = [];
    flicker_base_plvs = [];
    flicker_base_plvs_adj = [];
    
    estim_phase_vectors = [];
    estim_plvs = [];
    estim_plvs_adj = [];
    
    flicker_offset_phase_vectors = [];
    flicker_offset_plvs = [];
    flicker_offset_plvs_adj = [];
    

    % Loop through all of the neurons
    for f_nr = fieldnames(popul_data.nr_name)'
        f_nr = f_nr{1};

        % Get the vector of all phases for all trials at the stimulation frequency
        filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
        applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
        partial_apply = applyFunToColi(filt_trial, popul_data.raw_vm.(f_nr));

        vm_phases = arrayfun(partial_apply, [1:size(popul_data.raw_vm.(f_nr), 2)]', 'UniformOutput', false);
        vm_phases = cat(3, vm_phases{:});

        % Create a raster of the flicker onset rising edge
        flicker_onset_raster = zeros(size(timeline));
        flicker_rise_times = Multi_func.get_ephys_rise_times(timeline', popul_data.flicker_raster');
        [~, idxs] = ismember(flicker_rise_times, timeline);
        flicker_onset_raster(idxs) = 1;
        flicker_onset_alltrial_rasters = repmat(flicker_onset_raster', 1, size(vm_phases, 3));
        
        % Setup rasters
        flicker_base_rasters = zeros(size(flicker_onset_alltrial_rasters));
        flicker_base_rasters(flicker_base_idxs, :) = flicker_onset_alltrial_rasters(flicker_base_idxs, :);

        estim_rasters = zeros(size(flicker_onset_alltrial_rasters));
        estim_rasters(estim_idxs, :) = flicker_onset_alltrial_rasters(estim_idxs, :);

        flicker_offset_rasters = zeros(size(flicker_onset_alltrial_rasters));
        flicker_offset_rasters(flicker_offset_idxs, :) = flicker_onset_alltrial_rasters(flicker_offset_idxs, :);

        % Calculate PLVs
        [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, flicker_base_rasters, 0, 10);             
        flicker_base_plvs = [flicker_base_plvs, PLV(:)];
        flicker_base_plvs_adj = [flicker_base_plvs_adj, PLV2(:)];
        flicker_base_phase_vectors = [flicker_base_phase_vectors; norm_vecs];

        [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, estim_rasters, 0, 10);             
        estim_plvs = [estim_plvs, PLV(:)];
        estim_plvs_adj = [estim_plvs_adj, PLV2(:)];
        estim_phase_vectors = [estim_phase_vectors; norm_vecs];

        [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases, flicker_offset_rasters, 0, 10);             
        flicker_offset_plvs = [flicker_offset_plvs, PLV(:)];
        flicker_offset_plvs_adj = [flicker_offset_plvs_adj, PLV2(:)];
        flicker_offset_phase_vectors = [flicker_offset_phase_vectors; norm_vecs];
    end

    % Save the PLVs
    data.(f_stim).flicker_base_plvs = flicker_base_plvs;
    data.(f_stim).flicker_base_plvs_adj = flicker_base_plvs_adj;
    data.(f_stim).flicker_base_phase_vectors = flicker_base_phase_vectors;
    
    data.(f_stim).estim_plvs = estim_plvs;
    data.(f_stim).estim_plvs_adj = estim_plvs_adj;
    data.(f_stim).estim_phase_vectors = estim_phase_vectors;

    data.(f_stim).flicker_offset_plvs = flicker_offset_plvs;
    data.(f_stim).flicker_offset_plvs_adj = flicker_offset_plvs_adj;
    data.(f_stim).flicker_offset_phase_vectors = flicker_offset_phase_vectors;
end

%% Calculate the phase-locking value between flicker and Vm across frequencies
% for individual neurons by calculating individual trial PLVs
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;
    avg_Fs = 1/mean(diff(popul_data.interp_time));

    % Grab the flicker idxs
    flicker_base_idxs = find(timeline > 0 & timeline < 1);
    estim_idxs = find(timeline > 1 & timeline < 2);
    flicker_offset_idxs = find(timeline > 2 & timeline < 3);

    % Loop through all of the neurons
    for f_nr = fieldnames(popul_data.nr_name)'
        f_nr = f_nr{1};
        
        % Store the neuron's PLVs in a structure
        % TODO

        % Get the vector of all phases for all trials at the stimulation frequency
        filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
        applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
        partial_apply = applyFunToColi(filt_trial, popul_data.raw_vm.(f_nr));

        vm_phases = arrayfun(partial_apply, [1:size(popul_data.raw_vm.(f_nr), 2)]', 'UniformOutput', false);
        vm_phases = cat(3, vm_phases{:});

         % Create a raster of the flicker onset rising edge
        flicker_onset_raster = zeros(size(timeline));
        flicker_rise_times = Multi_func.get_ephys_rise_times(timeline', popul_data.flicker_raster');
        [~, idxs] = ismember(flicker_rise_times, timeline);
        flicker_onset_raster(idxs) = 1;
        flicker_onset_alltrial_rasters = repmat(flicker_onset_raster', 1, size(vm_phases, 3));
        
        % Setup rasters
        flicker_base_rasters = zeros(size(flicker_onset_alltrial_rasters));
        flicker_base_rasters(flicker_base_idxs, :) = flicker_onset_alltrial_rasters(flicker_base_idxs, :);

        estim_rasters = zeros(size(flicker_onset_alltrial_rasters));
        estim_rasters(estim_idxs, :) = flicker_onset_alltrial_rasters(estim_idxs, :);

        flicker_offset_rasters = zeros(size(flicker_onset_alltrial_rasters));
        flicker_offset_rasters(flicker_offset_idxs, :) = flicker_onset_alltrial_rasters(flicker_offset_idxs, :);

        %-- PLVs for each trial that is stored within
        % 'vm_phases' with size [freq, timepoint, trial]
        
        
        % Calculate the base flicker PLVs
        partial_apply_PLV = @(tr) Multi_func.spike_field_PLV(vm_phases(:, :, tr), flicker_base_rasters(:, tr), 0, 6);

        [PLVs, PLV2s, norm_vecs] = arrayfun(partial_apply_PLV, [1:size(vm_phases, 3)], 'UniformOutput',false);
        PLV2s = cellfun(@(ar) ar(:), PLV2s, 'UniformOutput', false);
        PLV2s = [PLV2s{:}];

        data.(f_stim).nr_plvs.(f_nr).base_PLV2s = PLV2s;

        % Calculate the stim flicker PLVs
        partial_apply_PLV = @(tr) Multi_func.spike_field_PLV(vm_phases(:, :, tr), estim_rasters(:, tr), 0, 6);

        [PLVs, PLV2s, norm_vecs] = arrayfun(partial_apply_PLV, [1:size(vm_phases, 3)], 'UniformOutput',false);
        PLV2s = cellfun(@(ar) ar(:), PLV2s, 'UniformOutput', false);
        PLV2s = [PLV2s{:}];

        data.(f_stim).nr_plvs.(f_nr).estim_PLV2s = PLV2s;

        % Calculate the offset flicker PLVs
        partial_apply_PLV = @(tr) Multi_func.spike_field_PLV(vm_phases(:, :, tr), flicker_offset_rasters(:, tr), 0, 6);

        [PLVs, PLV2s, norm_vecs] = arrayfun(partial_apply_PLV, [1:size(vm_phases, 3)], 'UniformOutput',false);
        PLV2s = cellfun(@(ar) ar(:), PLV2s, 'UniformOutput', false);
        PLV2s = [PLV2s{:}];

        data.(f_stim).nr_plvs.(f_nr).offset_PLV2s = PLV2s;
    end
end

%% Plot the neuron population phase-locking values as mean+SEM
freqs = Multi_func.entr_freqs;
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;
    
    % Plot the figure with different colors between flicker base, estim,
    % and flicker offset
    
    figure;
    
    % Plotting the flicker baseline
    flicker_base_mean = mean(popul_data.flicker_base_plvs_adj, 2, 'omitnan');
    flicker_base_std = std(popul_data.flicker_base_plvs_adj, [], 2, 'omitnan');
    num_flick_base = size(popul_data.flicker_base_plvs_adj, 2);
    flicker_base_sem = flicker_base_std./sqrt(num_flick_base);
    
    fill_h = fill([freqs, flip(freqs)], [flicker_base_mean' + flicker_base_sem', ...
        flip(flicker_base_mean' - flicker_base_sem')], [0.5 0.5 0.5]);
    Multi_func.set_fill_properties(fill_h);

    Multi_func.set_default_axis(gca);
    hold on;
    plot(freqs, flicker_base_mean, 'color', Multi_func.base_color, DisplayName='Base');
    hold on;

    % Plotting the estim
    estim_mean = mean(popul_data.estim_plvs_adj, 2, 'omitnan');
    estim_std = std(popul_data.estim_plvs_adj, [], 2, 'omitnan');
    num_flick_base = size(popul_data.estim_plvs_adj, 2);
    estim_sem = estim_std./sqrt(num_flick_base);
    
    fill_h = fill([freqs, flip(freqs)], [estim_mean' + estim_sem', ...
        flip(estim_mean' - estim_sem')], [0.5 0.5 0.5]);
    Multi_func.set_fill_properties(fill_h);

    Multi_func.set_default_axis(gca);
    hold on;
    plot(freqs, estim_mean, 'color', Multi_func.stim_color, DisplayName='Estim');
    hold on;

    % Plotting the flicker offset
    flicker_offset_mean = mean(popul_data.flicker_offset_plvs_adj, 2, 'omitnan');
    flicker_offset_std = std(popul_data.flicker_offset_plvs_adj, [], 2, 'omitnan');
    num_flick_base = size(popul_data.flicker_offset_plvs_adj, 2);
    flicker_offset_sem = flicker_offset_std./sqrt(num_flick_base);
    
    fill_h = fill([freqs, flip(freqs)], [flicker_offset_mean' + flicker_offset_sem', ...
        flip(flicker_offset_mean' - flicker_offset_sem')], [0.5 0.5 0.5]);
    Multi_func.set_fill_properties(fill_h);

    Multi_func.set_default_axis(gca);
    hold on;
    plot(freqs, flicker_offset_mean, 'color', Multi_func.post_color, DisplayName='Offset');
    hold on;

    ax = gca;
        
    %set(ax,'Xscale','log');
    xlabel('Frequency (Hz)');
    ylabel('PLV');
    
    % Show only plots that have DisplayNames
    lines = findall(gca, 'Type', 'Line');
    lines = lines(~cellfun(@isempty, get(lines, 'DisplayName')));
    legend(lines, 'Location', 'east');
   
    title([f_stim(3:end)], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_pop_' f_stim '.png']);
    saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_pop_' f_stim '.pdf']);
        
end

%% Plot the population PLVs with individual neuron PLVs as points
freqs = Multi_func.entr_freqs;
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;
    
    figure;
    plot(freqs, popul_data.flicker_base_plvs_adj, '.', 'color', Multi_func.base_color, 'DisplayName', 'Base');
    hold on;
    plot(freqs, popul_data.estim_plvs_adj, '.', 'color', Multi_func.stim_color, 'DisplayName', 'Stim');
    hold on;
    plot(freqs, popul_data.flicker_offset_plvs_adj, '.', 'color', Multi_func.post_color, 'DisplayName', 'Post');
    hold on;
    
    ax = gca;
        
    %set(ax,'Xscale','log');
    xlabel('Frequency (Hz)');
    ylabel('PLV');

    
    title([f_stim(3:end)], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_Indi_' f_stim '.png']);
    saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_Indi_' f_stim '.pdf']);
end

%% Plot the individual neuron phase-locking values with trials as mean+SEM
% TODO change values for this

% Set variable x-axis limit
%set_xlim = [0, 200];
set_xlim = [0, 10];

freqs = Multi_func.entr_freqs;
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;

    % Plot the figure with different colors between flicker base, estim,
    % and flicker offset
    % Loop through all of the neurons
    for f_nr = fieldnames(popul_data.nr_name)'
        f_nr = f_nr{1};

        figure;

        % Plotting the flicker baseline
        flicker_base_mean = mean(popul_data.nr_plvs.(f_nr).base_PLV2s, 2, 'omitnan');
        flicker_base_std = std(popul_data.nr_plvs.(f_nr).base_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).base_PLV2s, 2);
        flicker_base_sem = flicker_base_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [flicker_base_mean' + flicker_base_sem', ...
            flip(flicker_base_mean' - flicker_base_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, flicker_base_mean, 'color', Multi_func.base_color, DisplayName='Base');
        hold on;

        % Plotting the estim
        estim_mean = mean(popul_data.nr_plvs.(f_nr).estim_PLV2s, 2, 'omitnan');
        estim_std = std(popul_data.nr_plvs.(f_nr).estim_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).estim_PLV2s, 2);
        estim_sem = estim_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [estim_mean' + estim_sem', ...
            flip(estim_mean' - estim_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, estim_mean, 'color', Multi_func.stim_color, DisplayName='Estim');
        hold on;

        % Plotting the flicker offset
        flicker_offset_mean = mean(popul_data.nr_plvs.(f_nr).offset_PLV2s, 2, 'omitnan');
        flicker_offset_std = std(popul_data.nr_plvs.(f_nr).offset_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).offset_PLV2s, 2);
        flicker_offset_sem = flicker_offset_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [flicker_offset_mean' + flicker_offset_sem', ...
            flip(flicker_offset_mean' - flicker_offset_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, flicker_offset_mean, 'color', Multi_func.post_color, DisplayName='Offset');
        hold on;

        ax = gca;

        %set(ax,'Xscale','log');
        xlabel('Frequency (Hz)');
        ylabel('PLV');
        xlim(set_xlim);


        % Show only plots that have DisplayNames
        lines = findall(gca, 'Type', 'Line');
        lines = lines(~cellfun(@isempty, get(lines, 'DisplayName')));
        legend(lines, 'Location', 'east');

        title([ f_nr ' ' f_stim(3:end)], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_pop_' f_stim '_' f_nr '_xlimUp' num2str(set_xlim) '.png']);
        saveas(gcf, [figure_path 'Flicker' f 'PLV' f 'PLV_pop_' f_stim '_' f_nr '_xlimUp' num2str(set_xlim) '.pdf']);
    end
end