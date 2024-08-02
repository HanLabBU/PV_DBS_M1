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
freqs = Multi_func.entr_freqs;
for f_region = fieldnames(region_data)'
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
        base_phase_vectors = [];
        stim_phase_vectors = [];
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
            partial_apply = applyFunToColi(filt_trial, stim_data.all_trial_SubVm{nr});
            % Iterate through each column trial and concatenate all o fthe results
            vm_phases{nr} = arrayfun(partial_apply, [1:size(stim_data.all_trial_SubVm{nr}, 2)]' , 'UniformOutput', false); 
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
            base_phase_vectors = [base_phase_vectors; norm_vecs];

            [PLV, PLV2, norm_vecs] = Multi_func.spike_field_PLV(vm_phases{nr}, stim_rasters, 0, 10);           
            stim_plvs(:, nr) = PLV(:);
            stim_plvs_adjusted(:, nr) = PLV2(:); %TODO I changed this to transpose the dimensions
            stim_phase_vectors = [stim_phase_vectors; norm_vecs];
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

%% Loop and plot all of the dbs-vm stuff
freqs = [1:200];
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Loop through stim conditions
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_data = data_bystim.(f_stim);
        
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

        % Plot stim data with error bars
        stim_plvs_mean = nanmean(stim_data.stim_dbsvm_plvs_adj, 1);
        stim_plvs_std = nanstd(stim_data.stim_dbsvm_plvs_adj, 1);
        num_stim_plvs = size(stim_data.stim_dbsvm_plvs_adj, 1);
        stim_plvs_sem = stim_plvs_std./sqrt(num_stim_plvs);
        
        fill_h = fill([freqs, flip(freqs)], [[stim_plvs_mean + stim_plvs_sem], flip(stim_plvs_mean - stim_plvs_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, stim_plvs_mean, 'color', Multi_func.stim_color);
        hold on;
        
        ax = gca;
        set(ax,'Xscale','log');
        legend({'base', 'stim'}, 'Location', 'west');
        ax.Units = 'centimeters';
        ax.InnerPosition = [2 2 3.91 3.24];

        %ylim([-0.07 0.5]);
        title([f_region(3:end) ' ' f_stim(3:end)], 'Interpreter', 'none');
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
