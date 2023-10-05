% This uses the modulation factor (from Aryn Gitis's paper):
%(FR_stim - FR_base)/(FR_stim + FR_base)


clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';

pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

exclude_200ms = 1;

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end
%Load the data
load(save_all_data_file);

% Calculate the modulation factor between baseline and stim for both Vm and Firing rate per neuron
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Generate figure
    figure('Position', [100 100 1500 1000]);
    %tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 17.56, 3.17]);
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
    
    for f_stim=stims'
        f_stim = f_stim{1};
        
        Vm_mod = [];
        FR_mod = [];
        % Loop through each neuron to calculate the baseline vs stim modulation factor
        for nr=1:size(data_bystim.(f_stim).neuron_Vm, 2)
            trace_time = data_bystim.(f_stim).trace_timestamps(:, nr);
            stim_time = data_bystim.(f_stim).stim_timestamps(:, nr);
            base_idx = find(trace_time < stim_time(1));
            stim_idx = find(trace_time >= stim_time(1) & trace_time <= stim_time(end));

            % calculate the Vm modulation factor
            base_vm = mean(data_bystim.(f_stim).neuron_Vm(base_idx, nr), 'omitnan');
            stim_vm = mean(data_bystim.(f_stim).neuron_Vm(stim_idx, nr), 'omitnan');
            cur_vm_mod = (stim_vm - base_vm)/(stim_vm + base_vm);
            Vm_mod(end + 1) = cur_vm_mod;

            % calculate the FR modulation factor
            Fs = data_bystim.(f_stim).framerate(nr);
            base_FR = mean(Fs*data_bystim.(f_stim).neuron_spikecounts_raster(base_idx, nr), 'omitnan');
            stim_FR = mean(Fs*data_bystim.(f_stim).neuron_spikecounts_raster(stim_idx, nr), 'omitnan');
            cur_fr_mod = (stim_FR - base_FR)/(stim_FR + base_FR);
            FR_mod(end + 1) = cur_fr_mod;
        end

        % Calculate the correlation regression for both variables
        p = polyfit(Vm_mod, FR_mod, 1);
        mdl = fitlm(Vm_mod, FR_mod)
        nexttile;
        plot(Vm_mod, FR_mod, '.');
        hold on;
        x1 = linspace(max(Vm_mod), min(Vm_mod), 100);
        plot(x1, polyval(p, x1));
        xlabel('Vm MF');
        ylabel('Firing Rate MF');
        title([f_stim(3:end) ' MF Correlation']);
        disp([f_region(3:end) ' ' f_stim(3:end) ' correlation coefficient']);
        
        %disp(xcorr(Vm_mod, FR_mod, 'coeff' ));
        %r = xcorr(Vm_mod, FR_mod, 'coeff' );
        %disp(xcorr(Vm_mod, FR_mod, 'coeff' ).^2);

    end

    sgtitle([f_region(3:end) ' MF Correlation']);
end
