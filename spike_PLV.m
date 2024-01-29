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

%% Loop through and calculate spike-Vm PLV values for all region and conditions
freqs = [1:200];
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        base_plvs = [];
        stim_plvs = [];
        base_plvs_adjusted = [];
        stim_plvs_adjusted = [];
        
        tic;
        parfor nr = 1:length(stim_data.all_trial_SubVm)
            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            partial_apply = applyFunToColi(filt_trial, stim_data.all_trial_SubVm{nr});
            
            vm_phases{nr} = arrayfun(partial_apply, [1:size(stim_data.all_trial_SubVm{nr}, 2)]' , 'UniformOutput', false); 
            % 
            vm_phases{nr} = cat(3, vm_phases{nr}{:});
        
            time = stim_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(vm_phases{nr}, 3));

            base_idx = find(time < Multi_func.base_ped(2)/1000);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            
            spike_rasters = stim_data.all_trial_spike_rasters{nr};
            base_rasters = zeros(size(spike_rasters));
            stim_rasters = zeros(size(spike_rasters));
            
            base_rasters(base_idx) = spike_rasters(base_idx);
            stim_rasters(stim_idx) = spike_rasters(stim_idx);
            
            [PLV, PLV2] = Multi_func.spike_field_PLV(vm_phases{nr}, base_rasters, 7);             
            base_plvs(nr, :) = PLV;
            base_plvs_adjusted(nr, :) = PLV2;


            [PLV, PLV2] = Multi_func.spike_field_PLV(vm_phases{nr}, stim_rasters, 7);             
            stim_plvs(nr, :) = PLV;
            stim_plvs_adjusted(nr, :) = PLV2;
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).base_spikevm_plvs_adj = base_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_spikevm_plvs_adj = stim_plvs_adjusted;
        
        region_data.(f_region).(f_stim).base_spikevm_plvs = base_plvs;
        region_data.(f_region).(f_stim).stim_spikevm_plvs = stim_plvs;
        toc
    end
end

%% Loop and plot all of the spike-vm stuff
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        
        figure;
        plot(freqs, nanmean(stim_data.base_spikevm_plvs_adj, 1));
        hold on;
        plot(freqs, nanmean(stim_data.stim_spikevm_plvs_adj, 1));

        legend({'base', 'stim'});
        ylim([-0.07 0.5]);
        %xlim([2, 10]);
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');

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

%% Loop through and calculate spike-DBS PLV values for all region and conditions
freqs = [1:200];
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        stim_plvs = [];
        stim_plvs_adjusted = [];
        
        tic;
        for nr = 1%:length(stim_data.all_trial_SubVm)
            % Creates a dbs signal just using a one at pulse times
            %dbs_signal = zeros(size(stim_data.trace_timestamps(:, nr)));
            %pulse_time = stim_data.stim_timestamps(:, nr);
            %dbs_start = stim_data.stim_timestamps(1, nr) - stim_data.trace_timestamps(1, nr);
            %stim_idx =  ceil(((dbs_start + pulse_time)*avg_Fs) + 1); 
            %dbs_signal(stim_idx) = 1;
            %dbs_signal = repmat(dbs_signal, 1, size(stim_data.all_trial_SubVm{nr}, 2));
            pulse_time = stim_data.stim_timestamps(:, nr);
            period = 2*pi/(mean(diff(pulse_time)));
            dbs_signal = zeros(size(stim_data.trace_timestamps(:, nr)));
            time = stim_data.trace_timestamps(:, nr);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            % Still of a bit, but the 
            dbs_signal(stim_idx) = cos(period*time(stim_idx));
            
            %DEBUG
            pause;

            filt_trial = @(trial) (angle(Multi_func.filt_data(trial, freqs, avg_Fs)));
            
            applyFunToColi = @(func, mat) (@(col) func(mat(:, col)));
            partial_apply = applyFunToColi(filt_trial, dbs_signal);
            

            dbs_phases{nr} = arrayfun(partial_apply, [1:size(dbs_signal, 2)]' , 'UniformOutput', false); 
            % 
            dbs_phases{nr} = cat(3, dbs_phases{nr}{:});
        
            time = stim_data.trace_timestamps(:, nr);
            time = repmat(time, 1, size(dbs_phases{nr}, 3));

            base_idx = find(time < Multi_func.base_ped(2)/1000);
            stim_idx = find(time >= Multi_func.stim_ped(1)/1000 & time <= Multi_func.stim_ped(2)/1000);
            
            spike_rasters = stim_data.all_trial_spike_rasters{nr};
            stim_rasters = zeros(size(spike_rasters));
            stim_rasters(stim_idx) = spike_rasters(stim_idx);

            [PLV, PLV2] = Multi_func.spike_field_PLV(dbs_phases{nr}, stim_rasters, 7);             
            stim_plvs(nr, :) = PLV;
            stim_plvs_adjusted(nr, :) = PLV2;
        end

        % Save the variables to structs
        region_data.(f_region).(f_stim).stim_spikedbs_plvs_adj = stim_plvs_adjusted;
        region_data.(f_region).(f_stim).stim_spikedbs_plvs = stim_plvs;
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
        stim_data = data_bystim.(f_stim);
        
        figure;
        plot(freqs, nanmean(stim_data.stim_spikedbs_plvs_adj, 1));
        hold on;
        yline(0);

        legend({'base', 'stim'});
        %ylim([-0.07 0.5]);
        %xlim([2, 10]);
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');

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
%cellfun(neuron_plv, stim_data.all_trial_SubVm);
