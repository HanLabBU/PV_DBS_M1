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

% Parameters for frames to chop off
front_frame_drop = 15;
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

% Smoothing parameter for spike rate
% TODO could use a smaller window size
srate_win = 20;

% Do all region
all_region = 1;

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Load all of the data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
load(save_all_data_file);

% Combine all regions into 1
if all_region == 1
    region_data = Multi_func.combine_regions(region_data);
end

avg_Fs = mean(region_data.r_combine.f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%set figures off
%set(0,'DefaultFigureVisible','off');

spacing = 1.4;

% Show all of the neuron's trial-averaged Vm
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
        
    figure('Renderer', 'Painters', 'Position', [200 200 1200 740]);    
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');    
    
    for f_stim=stims'
        f_stim = f_stim{1};
     
        nexttile;
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)*1000; % Convert to ms
        % Plot each neuron's average Vm
        neuron_vm = data_bystim.(f_stim).neuron_Vm;
        neuron_vm = (neuron_vm - min(neuron_vm, [], 1))./ (max(neuron_vm, [], 1) - min(neuron_vm, [], 1));
        plot(timeline, neuron_vm + [spacing:spacing:size(neuron_vm, 2)*spacing], '-k');
        hold on;

        %-- Find points with higher than 2 std baseline
        base_vm = [];    
        % Grab the baseline sub Vm for all neurons    
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)    
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            base_vm = horzcat_pad(base_vm, neuron_vm(baseline_idx, i) );
        end 
        
        % Calculate the baseline std for each neuron
        avg_base_vm = mean(base_vm, 1, 'omitnan');
        base_std = std(base_vm, 0, 1, 'omitnan');

        [row, col] = find(neuron_vm > 2*base_std + avg_base_vm);
        
        plot(timeline(row), spacing*.75 + col*spacing, '.r', 'MarkerSize', 6);

        title(f_stim, 'Interpreter', 'none');
    end
end
