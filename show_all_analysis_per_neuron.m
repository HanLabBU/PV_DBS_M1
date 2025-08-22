clc;
clear all;
close all;

%%
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


%%
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
extra_trace = 3;

%% Loop through and test the power spectra filtering
%set(0, 'DefaultFigureWindowStyle', 'docked');
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through each stimulation parameter
    for f_stim= stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        pop_spec = cellfun(@(spec_data) mean(abs(spec_data), 3, 'omitnan'),...
            popul_data.neuron_hilbfilt, 'UniformOutput', false);
        
        avg_spec = double(mean(cat(3, pop_spec{:}), 3, 'omitnan'));
        
        time = mean(popul_data.trace_timestamps, 2, 'omitnan');

        %figure;

        %imagesc(avg_spec);
        %imagesc('XData', time', 'YData', Multi_func.entr_freqs, 'ZData', avg_spec);
        %imagesc('XData', 1:2316, 'YData', 1:200, 'ZData', avg_spec);
        %surface(time, Multi_func.entr_freqs, avg_spec, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');

        % tiled layout
        % 1. heatmap trials
        % 2. Spike Raster
        % 3. Sub Vm average
        % 4. Pulse Triggered average, for each period of stimulation
        % 5. POwer Spectra
        % 6. stim-Vm PLV
        % 7. 

        % Plot all single cell data in each tiledlayout
        for nr=1:length(popul_data.neuron_name)
            time = popul_data.trace_timestamps(:, nr);

            figure('Position', [0 0 1000 1500]);
            tiledlayout(6, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

            % Plotting all subVm heatmap
            nexttile;
            nr_vm = popul_data.all_trial_SubVm{nr};
            %imagesc(nr_vm');
            size(time)
            size([1:size(nr_vm, 2)])
            size(nr_vm')

            imagesc('XData', time, 'YData', 1:size(nr_vm, 2), 'CData', nr_vm');
            Multi_func.set_default_axis(gca);
            xlim([min(time) max(time)]);
            colorbar;
            title('All Trial Vm');

            % Plotting the spike raster
            nexttile;
            [row, col] = find(popul_data.all_trial_spike_rasters{nr});
            plot(time(row), col, '.');
            Multi_func.set_default_axis(gca);
            
            xlim([min(time) max(time)]);
            legend(['Num Spikes: ' num2str(sum(popul_data.all_trial_spike_rasters{nr} , 'all'))]);
            title('Spike Raster');

            % Plotting Sub Vm average
            nexttile;
            avg_vm = mean(nr_vm, 2, 'omitnan');
            sem_vm = std(nr_vm, 0, 2, 'omitnan')./sqrt(size(nr_vm, 2));
            fill_h = fill([time; flipud(time)], ...
                [avg_vm - sem_vm ; flipud(avg_vm + sem_vm)], ...
                [0.5 0.5 0.5]);
            Multi_func.set_fill_properties(fill_h);

            Multi_func.set_default_axis(gca);
            hold on;
            plot(time, avg_vm, '-k');

            xlim([min(time) max(time)]);
            title('Average Vm');
        
            % plotting the pulse triggered average for each period of stim
            nexttile;
            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width_time = mean(diff(popul_data.stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(popul_data.trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width_time./nr_avg_trace_time);
            
            % Grab all pulse Vm across all trials

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
            nr_pulse_vm = cell2mat(nr_norm_vm);
            
            avg_pulse_vm = mean(nr_pulse_vm, 2);
            sem_pulse_vm = std(nr_pulse_vm, 0, 2)./sqrt(size(nr_pulse_vm, 2));
            time_pul = ([0:size(nr_pulse_vm, 1) - 1] - extra_trace)*1000./avg_Fs;
            
            fill_h = fill([time_pul'; flipud(time_pul')], ...
                [avg_pulse_vm - sem_pulse_vm ; flipud(avg_pulse_vm + sem_pulse_vm)], ...
                [0.5 0.5 0.5]);
            Multi_func.set_fill_properties(fill_h);
            Multi_func.set_default_axis(gca);

            hold on;
            plot(time_pul, avg_pulse_vm);
            xline([0 nr_avg_pulse_width_time*1000]);
            title('Pulse Triggered Average');
            
            % Plotting the power spectra
            nexttile;
            
            % Get the indices define the periods within in the trace
            baseline_idx = find(popul_data.trace_timestamps(:, nr) < popul_data.stim_timestamps(1, nr));
            stim_idx = find(popul_data.trace_timestamps(:, nr) >= popul_data.stim_timestamps(1, nr) & ...
                            popul_data.trace_timestamps(:, nr) <= popul_data.stim_timestamps(end, nr));
            
            nr_power_spec = popul_data.all_trial_power_spec{nr};
            base_power = mean(nr_power_spec(:, baseline_idx, :), 2, 'omitnan');
            stim_power = mean(nr_power_spec(:, stim_idx, :), 2, 'omitnan');
            
            norm_power_spec = (nr_power_spec - base_power)./(base_power + stim_power);
            avg_power_spec = mean(norm_power_spec, 3, 'omitnan');
            
            freq = mean(popul_data.all_trial_spec_freq{nr}, 3, 'omitnan');

            surface(time, freq, avg_power_spec, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            xlim([min(time) max(time)]);
            ylim([min(freq) max(freq)]);
            colorbar;
            Multi_func.set_default_axis(gca);

            title('Sub Vm Power Spectra');

            % Plotting the stim-Vm PLV (may need to run 'single_cell_mod.m')
            nexttile;
            plot(Multi_func.entr_freqs, popul_data.stim_dbsvm_plvs_adj(:, nr));
            Multi_func.set_default_axis(gca);
            xlabel('Frequency (Hz)');
            title('stim-Vm PLV');
            
            sgtitle([f_region ' ' f_stim ' ' popul_data.neuron_name{nr}], 'Interpreter', 'none');
            saveas(gcf, [figure_path  'Neuronwise' f 'single' f popul_data.neuron_name{nr} '.png']);
            savefig(gcf, [figure_path  'Neuronwise' f 'single' f popul_data.neuron_name{nr} '.fig']);
        end
    end
end
