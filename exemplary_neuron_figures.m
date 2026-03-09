% Housekeeping clear all stuff
clc;
clear all;
close all;
f = filesep;

% Other sample traces to try
%617100_M1_rec20211111_FOV3_140 trial #4 (has delta)
% 31556noeartag_M1_rec20221206_FOV1_140 (hyperpolarization) trial 3 or 9

% Specify which photobleaching detrending to use
% (2) for photobleaching only the baseline
% (1) for the photobleaching taking into account the stimulation bump
% (0) for the basic full trace expontential fit subtraction detrending

sophis_bleachdetrend = 1;

% Option for excluding or including first 200ms
exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;

% Maingear office computer
local_root_path = '/home/pierfier/Projects/';
server_root_path = '~/handata_server/';
% Local windows machine
%local_root_path = 'Z:\Local_Data\';
%server_root_path = 'X:\';
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
pv_data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Path to save the figures
%savefig_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f 'Exemplary' f];
% handata3 server path
%savefig_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'Plots' f];


% Check if the figure path exists
if ~exist(savefig_path)
    disp('Figure path not found');
    return;
end

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
%Load the data
load(save_all_data_file);

%% Setup color variables
[red_blue_color_cmap] = (cbrewer('div', 'RdBu',500));
red_blue_color_cmap(red_blue_color_cmap > 1) = 1;
red_blue_color_cmap(red_blue_color_cmap < 0) = 0;
red_blue_color_cmap = flipud(red_blue_color_cmap);

[spectral_cmap] = (cbrewer('div', 'Spectral',500));
spectral_cmap(spectral_cmap > 1) = 1;
spectral_cmap(spectral_cmap < 0) = 0;
spectral_cmap = flipud(spectral_cmap);



savefig_path = Multi_func.save_plot;

%% Read in data for M1
f_region = 'r_M1';
data_bystim = region_data.(f_region);
%% Info for a 140Hz neuron regular polarized <--------
neuron = '617100_M1_rec20211110_FOV4_140_60';
trial_num = 3;
nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
neuron_data = data_bystim.f_140;

%% Infor for a 140Hz neuron hyperpolarized
neuron = '31556noeartag_M1_rec20221206_FOV1_140';
trial_num = 1;
nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
neuron_data = data_bystim.f_140;

%% Info for a 40Hz example without delta
neuron = '50373_M1_rec20230801_FOV2_40_200';
trial_num = 8;
nr_idx = find(contains(data_bystim.f_40.neuron_name, neuron));
neuron_data = data_bystim.f_40;

%% Info for a 40Hz example with delta <--------
neuron = '617100_M1_rec20211110_FOV3_40_60'; 
trial_num = 4;
nr_idx = find(contains(data_bystim.f_40.neuron_name, neuron));
neuron_data = data_bystim.f_40;

%% Read in data for V1
f_region = 'r_V1';
data_bystim = region_data.(f_region);
%% Info for a 140Hz reguarly depolarized neuron <--------
neuron = '109558_V1_rec20240110_FOV1_140_90.mat_3';
trial_num = 1; %5;
nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
neuron_data = data_bystim.f_140;
%% Info for a 40Hz hyperpolarized neuron <--------
neuron = '96334_V1_rec20231120_FOV1_40_150.mat_1';
trial_num = 10;
nr_idx = find(contains(data_bystim.f_40.neuron_name, neuron));
neuron_data = data_bystim.f_40;

%% Plot lined example trial, all trial heatmaps, and trial-averaged
% This requires using neuron info from the above selected subsection
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
tiledlayout(4, 1, 'TileSpacing', 'tight', 'Padding', 'tight', 'Units', 'centimeters', 'InnerPosition', [4 5 7 8.5]);
timeline = neuron_data.trace_timestamps(:, nr_idx);

%-- Plot example trial
nexttile;
raw_tr = neuron_data.all_trial_rawVm{nr_idx}(:, trial_num);
tr_noise = neuron_data.all_trial_trace_noise{nr_idx}(trial_num);
plot(timeline, raw_tr./tr_noise, 'k');
hold on;

spike_idx = neuron_data.all_trial_spikeidx{nr_idx}(:, trial_num);
%spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
%spike_idx = spike_idx - front_frame_drop + 1;
%spike_idx = spike_idx + 1;
spike_idx(isnan(spike_idx)) = [];
plot(timeline(spike_idx), ones(size(spike_idx)).*(max(raw_tr(spike_idx)./tr_noise) + 1), '.r');
hold on;

% Plot the stim pulses
%stim_idx = neuron_data.stim_timestamps(:, nr_idx);
%plot(stim_idx, ones(size(stim_idx))*(max(raw_tr(spike_idx)./tr_noise) + 3.5), '|k');
%hold on;
%plot(timeline, ones(size(timeline))*(max(raw_tr(spike_idx)./tr_noise) + 3.5), '-k');
%hold on;

plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)
hold on;

xlim([min(timeline) max(timeline)]);

% Add DBS bar
hold on;
Multi_func.plot_dbs_bar([0, 1], (max(raw_tr(spike_idx)./tr_noise) + 3.5), '');

ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);

%-- Plot the fluorescence heatmap
nexttile;
vm_map = [];
spidx_map = [];
% Loop through each trial of given neuron
for i = 1:size(neuron_data.all_trial_rawVm{nr_idx}, 2)
    trace_noise = neuron_data.all_trial_trace_noise{nr_idx}(i);

    % Add the trace with SBR as value
    vm_map = [vm_map; neuron_data.all_trial_rawVm{nr_idx}(:, i)'./trace_noise];
    cur_spikeidx = neuron_data.all_trial_spikeidx{nr_idx}(:, i);
    spidx_map = [spidx_map; cur_spikeidx'];
end
imagesc('XData', timeline, 'YData', 1:size(vm_map, 1), 'CData', vm_map);
ylabel('Trial #');
cb  = colorbar;
cb.Label.String = 'SBR';
hold on;

%colormap(spectral_cmap);
colormap(Multi_func.red_purple_blue_color);
%colormap(Multi_func.get_red_blue_cmap());
%Old colorscheme
%colormap(jet*.8);

% Reset scale with the limits to the max value of whats plotted
cb = colorbar;
max_ax = max(abs(cb.Limits));
caxis([-max_ax max_ax]);

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(vm_map, 1)]);

ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

%-- Plot the population average trace
nexttile;

all_traces = neuron_data.all_trial_rawVm{nr_idx};
snr_traces = all_traces./neuron_data.all_trial_trace_noise{nr_idx};

% Calculate the average and SEM
avg_trace = nanmean(snr_traces, 2);
std_trace = nanstd(snr_traces, 0, 2);
num_trials = size(snr_traces, 2);
sem_trace = std_trace./sqrt(num_trials);

% Plot neuron's trial SEM
fill_h = fill([timeline; flip(timeline)], [avg_trace + sem_trace; flipud(avg_trace - sem_trace)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;

% Plot trial-averaged trace
plot(timeline, avg_trace, 'k');

hold on;

% Plotting the SNR scale
plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)

xlim([min(timeline) max(timeline)]);
ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);

%-- Plot the spectra for this single neuron
nexttile;
nr_pow = neuron_data.all_trial_rawvm_power_spec{nr_idx};
nr_idx

base_idxs = find(neuron_data.trace_timestamps(:, nr_idx) < neuron_data.stim_timestamps(1, nr_idx));
stim_idxs = find(neuron_data.trace_timestamps(:, nr_idx) >= neuron_data.stim_timestamps(1, nr_idx)& ...
neuron_data.trace_timestamps(:, nr_idx) <= neuron_data.stim_timestamps(end, nr_idx) );

base_pow = mean(nr_pow(:, base_idxs, :), 2);
stim_pow = mean(nr_pow(:, stim_idxs, :), 2);

nr_pow = (nr_pow - base_pow) ./(base_pow + stim_pow);

surface(timeline, ... 
        nanmean(neuron_data.neuron_rawvm_spec_freq(:, :, nr_idx), 2), ...
        mean(nr_pow, 3, 'omitnan'), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');

% Reset scale with the limits to the max value of whats plotted
cb = colorbar;
max_ax = max(abs(cb.Limits));
caxis([-max_ax max_ax]);

% Add colorbar
a = colorbar;
set(a, 'TickDirection', 'out');
a.Ticks = linspace(a.Limits(1), a.Limits(2), 3);
a.TickLabels = num2cell(round(linspace(a.Limits(1), a.Limits(2), 3), 1));
a.Label.String = 'Relative Power';

Multi_func.set_default_axis(gca);
xlabel('Time from Stim onset(s)');
xlim([-0.80 2.05]);
ylabel('Freq (Hz)');

% Save figure
saveas(gcf, [savefig_path 'Exemplary' f neuron '_heatmap.png']);
saveas(gcf, [savefig_path 'Exemplary' f neuron '_heatmap.pdf']);

%% -- Show single neuron example of phase locking to stimulation with summary pulse-triggered Vm average
% for that neuron

%% TODO find a neuron that was deemed entrained
f_region = 'r_V1';
data_bystim = region_data.(f_region);

%neuron = '23072_V1_rec20220407_FOV1_140_180_.mat_1'; % idx is 20
neuron = '109557_V1_rec20240110_FOV1_140_120.mat_7'; % idx 7

nr_idx = find(contains(data_bystim.f_140.neuron_name, neuron));
popul_data = data_bystim.f_140;

%% Plot single neuron's pulse-triggered Vm
timeline = [ [0:size(popul_data.nr_pulse_windows{nr_idx}, 1) - 1] - extra_trace]*1000./popul_data.framerate(nr_idx);

% Calculate the average Vm and the shuffled Vm
avg_vm = mean(popul_data.nr_pulse_windows{nr_idx}, 2, 'omitnan');
sem_vm = std(popul_data.nr_pulse_windows{nr_idx}, 0, 2, 'omitnan')./sqrt(size(popul_data.nr_pulse_windows{nr_idx}, 2));

% Calculate the shuffled distribution error bars
low_pulse_perc = prctile(popul_data.nr_shuf_pulses{nr_idx}, 2.5, 2);
high_pulse_perc = prctile(popul_data.nr_shuf_pulses{nr_idx}, 97.5, 2);
shuf_pulse_mean = mean(popul_data.nr_shuf_pulses{nr_idx}, 2);

% Find Vm that is significantly higher than the shuffled
sig_idx = find(avg_vm > high_pulse_perc);
% Find Vm that is significantly lower than the shuffled
sig_idx = [find(avg_vm < low_pulse_perc); sig_idx ];
sig_idx(sig_idx <= extra_trace + 1 | sig_idx > size(popul_data.nr_pulse_windows{nr_idx}, 1) - extra_trace - 1) = [];

figure;

% Plot pulse triggered Vm
fill_h = fill([timeline, flip(timeline)], [avg_vm + sem_vm; flipud(avg_vm - sem_vm)], [0.5 0.5 0.5]);
if ~isempty(fill_h)
    Multi_func.set_fill_properties(fill_h);
end
hold on;
plot(timeline, avg_vm, 'k', 'LineWidth', 1);
hold on;

% Plot the percentile Vm
plot(timeline, [low_pulse_perc high_pulse_perc], '--', 'Color', Multi_func.shuf_color);
hold on;
plot(timeline, shuf_pulse_mean, 'Color', Multi_func.shuf_color);
hold on;

% Plot siginificant points with larger markerSizes
plot(timeline(sig_idx), avg_vm(sig_idx), 'b.', 'MarkerSize', 10);
hold on;

% Plot the pulse bars
xline([timeline(extra_trace + 1), timeline(size(popul_data.nr_pulse_windows{nr_idx}, 1) - extra_trace)], ...
    'Color', Multi_func.pulse_color, 'LineWidth', 2);
hold on;
Multi_func.set_default_axis(gca);
ylabel('Normalized Vm');
xlabel('Time from pulse(ms)');

%% Plot single-cells PLV with frequency sweep
freqs = Multi_func.entr_freqs;

% Calculate the avg PLVs and SEMs
shuf_plvs_mean = mean(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs, 2, 'omitnan')';
num_plvs = size(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs ,2);
shuf_plvs_sem = std(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs, 0, 2, 'omitnan')'./sqrt(num_plvs);

% Calculate the percentiles
high_plv_prc = prctile(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs',97.5);
half_plv_prc = prctile(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs',50);
low_plv_prc = prctile(popul_data.shuf_plv_data(nr_idx).shuf_plv_adj_wfreqs',2.5);

figure;
% Plot the shuffled distribution for all frequencies
%fill_h = fill([freqs, flip(freqs)], [[shuf_plvs_mean + shuf_plvs_sem], flip(shuf_plvs_mean - shuf_plvs_sem)], [0.5 0.5 0.5]);
%        Multi_func.set_fill_properties(fill_h);
fill_h = fill([freqs, flip(freqs)], [[high_plv_prc], flip(low_plv_prc)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

hold on;
plot(freqs, half_plv_prc, 'color', Multi_func.shuf_color);
%plot(freqs, shuf_plvs_mean, 'color', Multi_func.shuf_color);
hold on;
plot(freqs, popul_data.stim_dbsvm_plvs_adj(:, nr_idx), 'color', Multi_func.stim_color);

Multi_func.set_default_axis(gca);
xlabel('Frequency (Hz)');
ylabel('Pulse-Vm PLV^2');

%% Plot each neuron their single cell pulse-triggered Vm and PLV sweep

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end
    
    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Loop through each neuron
        tic;
        for nr=1:length(popul_data.neuron_name)

            % Calculate the average Vm and the shuffled Vm
            avg_vm = mean(popul_data.nr_pulse_windows{nr}, 2, 'omitnan');
            sem_vm = std(popul_data.nr_pulse_windows{nr}, 0, 2, 'omitnan')./sqrt(size(popul_data.nr_pulse_windows{nr}, 2));

            % Calculate the shuffled distribution error bars
            low_pulse_perc = prctile(popul_data.nr_shuf_pulses{nr}, 2.5, 2);
            high_pulse_perc = prctile(popul_data.nr_shuf_pulses{nr}, 97.5, 2);
            shuf_pulse_mean = mean(popul_data.nr_shuf_pulses{nr}, 2);

            % Find Vm that is significantly higher than the shuffled
            sig_idx = find(avg_vm > high_pulse_perc);
            % Find Vm that is significantly lower than the shuffled
            sig_idx = [find(avg_vm < low_pulse_perc); sig_idx ];
            sig_idx(sig_idx <= extra_trace + 1 | sig_idx > size(popul_data.nr_pulse_windows{nr}, 1) - extra_trace - 1) = [];

            % Calculate the avg PLVs and SEMs
            shuf_plvs_mean = mean(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs, 2, 'omitnan')';
            num_plvs = size(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs ,2);
            shuf_plvs_sem = std(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs, 0, 2, 'omitnan')'./sqrt(num_plvs);

            % Calculate the percentiles
            high_plv_prc = prctile(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs',97.5);
            half_plv_prc = prctile(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs',50);
            low_plv_prc = prctile(popul_data.shuf_plv_data(nr).shuf_plv_adj_wfreqs',2.5);

            % -- Start figure plotting

            figure('Position', [200 200 1500 1000]);
            tiledlayout(3, 1, 'Units', 'centimeters', 'InnerPosition', [20 2 7.5 11]);
            
            % Plot the average Vm
            nexttile;
            timeline = popul_data.trace_timestamps(:, nr);

            plot(timeline, mean(popul_data.all_trial_rawVm{nr}./popul_data.neuron_spike_amp(nr), 2, 'omitnan'), 'k');
            %xlim([0 1]);
            Multi_func.set_default_axis(gca);
            xlabel('Time from onset (s)');
            ylabel('Normalized Vm');
            ax = gca;
            ax.XAxis.TickLabelGapOffset = -1;
            ax.YAxis.TickLabelGapOffset = -1;

            % Plot the PLV with frequency sweep
            nexttile;

            % Show the shuffled distribution for all frequencies
            %fill_h = fill([freqs, flip(freqs)], [[shuf_plvs_mean + shuf_plvs_sem], flip(shuf_plvs_mean - shuf_plvs_sem)], [0.5 0.5 0.5]);
            %        Multi_func.set_fill_properties(fill_h);
            fill_h = fill([freqs, flip(freqs)], [[high_plv_prc], flip(low_plv_prc)], [0.5 0.5 0.5]);
            Multi_func.set_fill_properties(fill_h);

            hold on;
            plot(freqs, half_plv_prc, 'color', Multi_func.shuf_color);
            %plot(freqs, shuf_plvs_mean, 'color', Multi_func.shuf_color);
            hold on;
            plot(freqs, popul_data.stim_dbsvm_plvs_adj(:, nr), 'color', Multi_func.stim_color);

            Multi_func.set_default_axis(gca);
            xlabel('Frequency (Hz)');
            ylabel('Pulse-Vm PLV^2');
            ax = gca;
            ax.XAxis.TickLabelGapOffset = -1;
            ax.YAxis.TickLabelGapOffset = -1;

            % Plot the pulse-triggered average
            nexttile;

            timeline = [ [0:size(popul_data.nr_pulse_windows{nr}, 1) - 1] - extra_trace]*1000./popul_data.framerate(nr);

            fill_h = fill([timeline, flip(timeline)], [avg_vm + sem_vm; flipud(avg_vm - sem_vm)], [0.5 0.5 0.5]);
            if ~isempty(fill_h)
                Multi_func.set_fill_properties(fill_h);
            end
            hold on;
            plot(timeline, avg_vm, 'k', 'LineWidth', 1);
            hold on;

            % Plot the percentile Vm
            plot(timeline, [ ...
                low_pulse_perc high_pulse_perc], '--', 'Color', Multi_func.shuf_color);
            hold on;
            plot(timeline, shuf_pulse_mean, 'Color', Multi_func.shuf_color);
            hold on;

            % Plot siginificant points with larger markerSizes
            plot(timeline(sig_idx), avg_vm(sig_idx), 'b.', 'MarkerSize', 10);
            hold on;

            % Plot the pulse bars
            xline([timeline(extra_trace + 1), timeline(size(popul_data.nr_pulse_windows{nr}, 1) - extra_trace)], ...
                'Color', Multi_func.pulse_color, 'LineWidth', 2);
            hold on;
            Multi_func.set_default_axis(gca);
            ylabel('Normalized Vm');
            xlabel('Time from pulse (ms)');
            ax = gca;
            ax.XAxis.TickLabelGapOffset = -1;
            ax.YAxis.TickLabelGapOffset = -1;
            
            sgtitle([num2str(popul_data.plv_mod_stats(nr).mod) ' ' ...
                f_stim(3:end) ' ' ...
                popul_data.neuron_name{nr}], 'Interpreter', 'none');

            % Set fontsize
            fontsize(gcf, 14, 'points');

            % Save the plots
            saveas(gcf, [figure_path 'Examplary' f 'Entrainment' f popul_data.neuron_name{nr} '_plv_pulse_trig.png']);
            saveas(gcf, [figure_path 'Examplary' f 'Entrainment' f popul_data.neuron_name{nr} '_plv_pulse_trig.pdf']);
        end
    end
end

% ----- Everything below here appears to be DEPRECATED code --------
%% Code to plot individual trials with raster plot and average
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
tiledlayout(3, 1, 'TileSpacing', 'loose', 'Padding', 'loose', 'Units', 'centimeters', 'InnerPosition', [4 5 8 15]);
timeline = neuron_data.trace_timestamps(:, nr_idx);


% Plot trial-averaged fluorescence
nexttile;
all_traces = neuron_data.all_trial_rawVm{nr_idx};
snr_traces = all_traces./neuron_data.all_trial_trace_noise{nr_idx};

% Calculate the average and SEM
avg_trace = nanmean(snr_traces, 2);
std_trace = nanstd(snr_traces, 0, 2);
num_trials = size(snr_traces, 2);
sem_trace = std_trace./sqrt(num_trials);

% Plot neuron's trial SEM
fill_h = fill([timeline; flip(timeline)], [avg_trace + sem_trace; flipud(avg_trace - sem_trace)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;

% Plot trial-averaged trace
plot(timeline, avg_trace, 'k');

hold on;

% Plot the stim pulses
stim_idx = neuron_data.stim_timestamps(:, nr_idx);
plot(stim_idx, ones(size(stim_idx))*(max(avg_trace + sem_trace) + 1), '|k');
hold on;
plot(timeline, ones(size(timeline))*(max(avg_trace + sem_trace) + 1), '-k');
hold on;

% Plotting the SNR scale
plot([min(timeline) min(timeline)], [-5 0], 'b', 'LineWidth', 2)

xlim([min(timeline) max(timeline)]);
ax = gca;
ax.YAxis.Visible = 'off';
Multi_func.set_default_axis(gca);
% Plot each individual trial trace
nexttile;
trace_map = [];
spidx_map = [];
% Loop through each trial of given neuron
for i = 1:size(neuron_data.all_trial_rawVm{nr_idx}, 2)
    tr_trace = neuron_data.all_trial_rawVm{nr_idx}(:, i)';
    norm_trace = (tr_trace - min(tr_trace)) ./(max(tr_trace) - min(tr_trace));
    trace_map = [trace_map; norm_trace + (i - 0.5)];
    cur_spikeidx = neuron_data.all_trial_spikeidx{nr_idx}(:, i);
    spidx_map = [spidx_map; cur_spikeidx'];
end
plot(timeline, trace_map, 'k');
hold on;
% loop through and plot the spikes detected
for i =1:size(spidx_map, 1)
    temp_spidx = spidx_map(i, :);
    temp_spidx(isnan(temp_spidx)) = [];
    plot(timeline(temp_spidx), trace_map(i, temp_spidx), '.r');
    hold on;
end

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(trace_map, 1) + 0.5]);

ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

% Save figure
saveas(gcf, [savefig_path 'Exemplary' f neuron '_tr_traces_average.png']);
saveas(gcf, [savefig_path 'Exemplary' f neuron '_tr_traces_average.pdf']);



%% Code to plot a raster plot
% Plot the raster of the current neuron

raster_map = NaN(size(spidx_map));
raster_map(find(~isnan(spidx_map))) = 1;

raster_map = raster_map.*(1:size(spidx_map, 1))';

% Clean up the nans
spidx_map(isnan(spidx_map)) = [];
raster_map(isnan(raster_map)) = [];

% Plot all of the raster points
plot(timeline(spidx_map), raster_map, '.k');

xlim([min(timeline) max(timeline)]);
ylim([0.5 size(trace_map, 1)]);
ax = gca;
set(ax, 'Color', 'none', 'Box', 'on', 'TickDir', 'out', 'linewidth', 0.2);

%% Get exemplary trace at 140 for V1
example_matfile = [pv_data_path '611284_V1_rec20210827_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 5;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['611284_V1_rec20210827_FOV1_140_60_.mat tr 5' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 15;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);



title('Exemplary V1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 140Hz Power Spectra');


%% Get 2nd exemplary trace at 140 for V1
example_matfile = [pv_data_path '611284_V1_rec20221108_FOV2_140_50_.mat'];
data = load(example_matfile);
trial_idx = 2;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['611284_V1_rec20221108_FOV2_140_50_.mat tr 2' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 24;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary V1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_140Hz_Trace2.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 140Hz Power Spectra');

%% Get exemplary M1 trace at 140

% Parameters for plot zoom in
base_zoom = [-600 -400];
stim_zoom = [10 210];
offset_zoom = [1400 1600];

example_matfile = [pv_data_path '617100_M1_rec20211111_FOV1_140_60_.mat'];
data = load(example_matfile);
trial_idx = 8;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure

figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
% Show file name
sgtitle(['617100_M1_rec20211111_FOV1_140_60_.mat tr 8'], 'Interpreter', 'none');

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 20;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 30, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -9;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-2, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);

% Plot boxes for each zoom in part
%baseline
dim = [];
dim(1) = base_idx(1);
dim(2) = -7;
dim(3) = range(base_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.base_color, 'LineStyle', '--');
hold on;

%stimulation
dim = [];
dim(1) = stim_ped_idx(1);
dim(2) = -7;
dim(3) = range(stim_ped_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.stim_color, 'LineStyle', '--');
hold on;

%offset
dim = [];
dim(1) = offset_idx(1);
dim(2) = -7;
dim(3) = range(offset_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.post_color, 'LineStyle', '--');

set(gca, 'Units', 'centimeters', 'Position', [5 20 10 3.00], 'PositionConstraint', 'innerposition');

title('Exemplary M1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.pdf']);


% Baseline zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
base_spike_idx = intersect(base_idx, spike_idx);

plot(detrend_traces(base_idx)./trace_noise, 'Color', Multi_func.base_color);
hold on;
plot(base_spike_idx - base_idx(1) + 1, detrend_traces(base_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy, ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 18]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Baseline Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_BaseZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_BaseZoomIn.pdf']);

% Stimulation zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
stim_spike_idx = intersect(stim_ped_idx, spike_idx);

plot(detrend_traces(stim_ped_idx)./trace_noise, 'Color', Multi_func.stim_color);
hold on;
plot(stim_spike_idx - stim_ped_idx(1) + 1, detrend_traces(stim_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Stimulation Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_StimZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_StimZoomIn.pdf']);


% Ofset zoom in of M1 140Hz trace
figure('renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
offset_spike_idx = intersect(offset_idx, spike_idx);

plot(detrend_traces(offset_idx)./trace_noise, 'Color', Multi_func.post_color);
hold on;
plot(offset_spike_idx - offset_idx(1) + 1, detrend_traces(offset_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 140 Hz Offset Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_OffsetZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_OffsetZoomIn.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 140Hz Power Spectra');

%------

%% Get 2nd exemplary M1 trace at 140
%example_matfile = [pv_data_path '617100_M1_rec20211110_FOV5_140_60_.mat'];
%trial_idx = 6;
example_matfile = [pv_data_path '617100_M1_rec20211110_FOV4_140_60_.mat'];
trial_idx = 3;
data = load(example_matfile);

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;


% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['617100_M1_rec20211110_FOV5_140_60_.mat tr 6'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 20;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '140 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -9;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary M1 140 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_140Hz_Trace2.pdf']);

%% Get exemplary trace at 40 Hz for V1
example_matfile = [pv_data_path '23072_V1_rec20220217_FOV3_40_220_.mat'];
data = load(example_matfile);
trial_idx = 3;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end

detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
spike_idx(find(spike_idx < front_frame_drop | spike_idx > back_frame_drop)) = [];
spike_idx = spike_idx - front_frame_drop + 1;
if length(spike_idx) == 0
    spike_idx = [NaN];
end

sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_time = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_time*sam_freq);
stim_idx = stim_idx - front_frame_drop + 1;

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['23072_V1_rec20220217_FOV3_40_220_.mat tr 3'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 15;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -10;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('Exemplary V1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace.pdf']);


%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 40Hz Power Spectra');


%% Get 2nd exemplary trace at 40 for V1
example_matfile = [pv_data_path '23072_V1_rec20220223_FOV2_40_250_.mat'];
data = load(example_matfile);
trial_idx = 8;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1

end 

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['23072_V1_rec20220223_FOV2_40_250_.mat tr 8'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 18;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -10;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('Exemplary V1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'V1_40Hz_Trace2.pdf']);


%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary V1 40Hz Power Spectra');

%% Get exemplary M1 trace at 40

% Parameters for plot zoom in
base_zoom = [-400 -200];
stim_zoom = [10 210];
offset_zoom = [1100 1300];

example_matfile = [pv_data_path '617100_M1_rec20211110_FOV3_40_60_.mat'];
data = load(example_matfile);
trial_idx = 9;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop),...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop));
end
detrend_traces = data.raw.trial{trial_idx}.raw_traces(front_frame_drop:back_frame_drop) - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_idx = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_idx*sam_freq);
stim_start = data.raw.trial{trial_idx}.raw_stimulation_time(1);

% Calculate the idx time for base, stim, and offset
base_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > base_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < base_zoom(2)./1000);
stim_ped_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > stim_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < stim_zoom(2)./1000);
offset_idx = find(data.align.trial{trial_idx}.camera_frame_time - stim_start > offset_zoom(1)./1000 & data.align.trial{trial_idx}.camera_frame_time - stim_start < offset_zoom(2)./1000);

% Generate figure
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1000]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['617100_M1_rec20211110_FOV3_40_60_.mat tr 9' ], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 22;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);

% Plot boxes for each zoom in part
%baseline
dim = [];
dim(1) = base_idx(1);
dim(2) = -7;
dim(3) = range(base_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.base_color, 'LineStyle', '--');
hold on;

%stimulation
dim = [];
dim(1) = stim_ped_idx(1);
dim(2) = -7;
dim(3) = range(stim_ped_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.stim_color, 'LineStyle', '--');
hold on;

%offset
dim = [];
dim(1) = offset_idx(1);
dim(2) = -7;
dim(3) = range(offset_idx);
dim(4) = 25;
rectangle('Position', dim, 'EdgeColor', Multi_func.post_color, 'LineStyle', '--');

set(gca, 'Units', 'centimeters', 'Position', [5 20 9 3.00], 'PositionConstraint', 'innerposition');

title('Exemplary M1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.pdf']);


% Baseline zoom in of M1 40Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
base_spike_idx = intersect(base_idx, spike_idx);

plot(detrend_traces(base_idx)./trace_noise, 'Color', Multi_func.base_color);
hold on;
plot(base_spike_idx - base_idx(1) + 1, detrend_traces(base_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy, ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 18]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Baseline Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_BaseZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_BaseZoomIn.pdf']);

% Stimulation zoom in of M1 40Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
stim_spike_idx = intersect(stim_ped_idx, spike_idx);

plot(detrend_traces(stim_ped_idx)./trace_noise, 'Color', Multi_func.stim_color);
hold on;
plot(stim_spike_idx - stim_ped_idx(1) + 1, detrend_traces(stim_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Stimulation Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_StimZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_StimZoomIn.pdf']);

% Ofset zoom in of M1 140Hz trace
figure('visible', 'on', 'renderer', 'painters', 'Position', [0 0 1200 1200]);
set(gca, 'Units', 'centimeters', 'Position', [15 20 2.86 3.00]);
offset_spike_idx = intersect(offset_idx, spike_idx);

plot(detrend_traces(offset_idx)./trace_noise, 'Color', Multi_func.post_color);
hold on;
plot(offset_spike_idx - offset_idx(1) + 1, detrend_traces(offset_spike_idx)./trace_noise, '.r', 'MarkerSize', 6);
hold on;

% Plot the SNR line reference
posx = -1;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 10, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 20; % Plotting 20 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 15]);

Multi_func.set_default_axis(gca);
axis off;

title('Exemplary M1 40 Hz Offset Zoom in');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_OffsetZoomIn.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_OffsetZoomIn.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 40Hz Power Spectra');


%% 2nd example M1 trace at 40
example_matfile = [pv_data_path '31556eartag_M1_rec20221214_FOV1_40_80_.mat'];
data = load(example_matfile);
trial_idx = 4;

if sophis_bleachdetrend == 2
    [baseline, coeff] = Multi_func.exp_fit_Fx_Base(data.raw.trial{trial_idx}.raw_traces,...
        round(data.align.trial{trial_idx}.camera_framerate));

elseif sophis_bleachdetrend == 1
    [baseline, coeff] = Multi_func.exp_fit_Fx(data.raw.trial{trial_idx}.raw_traces,...
        round(data.align.trial{trial_idx}.camera_framerate));
else
    [baseline, coeff] = Multi_func.exp_fit(data.raw.trial{trial_idx}.raw_traces);
end
detrend_traces = data.raw.trial{trial_idx}.raw_traces - baseline';
trace_noise = data.align.trial{trial_idx}.spike_info375.trace_noise;
spike_idx = data.align.trial{trial_idx}.spike_info375.spike_idx{1};
sam_freq = data.align.trial{trial_idx}.camera_framerate;
stim_idx = data.raw.trial{trial_idx}.raw_stimulation_time - data.raw.trial{trial_idx}.raw_camera_start_time;
stim_idx = round(stim_idx*sam_freq);

% Generate figure
figure('renderer', 'painters', 'Position', [0 0 1200 500]);

% Show file name
posx = 900;
posy = 20;
text(posx, posy, ['31556eartag_M1_rec20221214_FOV1_40_80_.mat tr 4'], 'Interpreter', 'none');
hold on;

%Plot the trace
plot(detrend_traces./trace_noise, 'k');
%set(gcf, 'color', 'none');
%set(gca, 'color', 'none');
axis off;
hold on;

% Plot spikes detected
plot(spike_idx, detrend_traces(spike_idx)./trace_noise, '.r', 'MarkerSize', 12);
hold on;

% Plot the stimulation time pulses
posx = 3;
posy = 18;
plot(stim_idx, repmat(posy, length(stim_idx), 1), '|k');
hold on;
plot([1, length(detrend_traces)], [posy posy], '-k');
hold on;
text(posx, posy + 1, '40 Hz Stimulation')

% Plot the SNR line reference
posx = -20;
posy = 0;
snr_scale = 5;
plot([posx posx], [posy posy+snr_scale], 'k', 'LineWidth', 2);
hold on;
ht = text(posx - 25, posy , ['Spike SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
hold on;

% Plot the timescale
posx = 0;
posy = -7;
time_scale = 500; % Plotting 500 ms
plot([posx posy+time_scale]*sam_freq/1000, [posy posy], 'k', 'LineWidth', 2);
hold on;
text(posx, posy-1, [num2str(time_scale) 'ms']);
ylim([posy-5 30]);
title('2nd Exemplary M1 40 Hz trace');
%saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.eps'], 'epsc');
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.png']);
saveas(gcf, [savefig_path 'Exemplary' f 'M1_40Hz_Trace2.pdf']);

%From same exemplary trace as above, show the individual spectrum
%signal = data.align.trial{trial_idx}.spike_info375.trace_ws;
%[wt, f] = Multi_func.get_power_spec(signal, sam_freq);
%
% Generate figure
%figure('renderer', 'painters', 'Position', [0 0 1200 500]);
%surface(data.align.trial{trial_idx}.camera_frame_time, f, abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%set(gca, 'color', 'none');
%title('Exemplary M1 40Hz Power Spectra');
