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
% Using a smaller window size
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
        
    figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);    
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');    
    
    for f_stim=stims'
        f_stim = f_stim{1};
     
        nexttile;
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2); 
        % Plot each neuron's average Vm
        neuron_vm = data_bystim.(f_stim).neuron_Vm;
        neuron_vm = (neuron_vm - min(neuron_vm, [], 1))./ (max(neuron_vm, [], 1) - min(neuron_vm, [], 1));
        plot(timeline, neuron_vm + [spacing:spacing:size(neuron_vm, 2)*spacing], '-k');
        hold on;

        %-- Find points with higher than 3 std baseline
        base_vm = [];    
        % Grab the baseline sub Vm for all neurons    
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)    
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            base_vm = horzcat_pad(base_vm, neuron_vm(baseline_idx, i) );
        end 
        
        % Calculate the baseline std for each neuron
        avg_base_vm = mean(base_vm, 1, 'omitnan');
        base_std = std(base_vm, 0, 1, 'omitnan');

        [row, col] = find(neuron_vm > 3*base_std + avg_base_vm);
        
        plot(timeline(row), spacing*.75 + col*spacing, '.b', 'MarkerSize', 6);
        xlim([-1 2.05]);
        Multi_func.set_default_axis(gca);
        
        title(f_stim, 'Interpreter', 'none');
    end

    % Save the figure
    saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_vm_base_comp.png']);
    %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_vm_base_comp.pdf']);
    %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_vm_base_comp.eps'], 'epsc');
end


% Show all of the neuron's trial-averaged Firing Rate
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
        
    figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);    
    tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');    
    
    for f_stim=stims'
        f_stim = f_stim{1};
     
        nexttile;
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2); 
        % Plot each neuron's average Firing Rate
        neuron_fr = data_bystim.(f_stim).neuron_srate_large; %TODO change this
        neuron_fr = (neuron_fr - min(neuron_fr, [], 1))./ (max(neuron_fr, [], 1) - min(neuron_fr, [], 1));
        plot(timeline, neuron_fr + [spacing:spacing:size(neuron_fr, 2)*spacing], '-k');
        hold on;

        %-- Find points with higher than 3 std baseline
        base_fr = [];    
        % Grab the baseline Firing rate for all neurons    
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)    
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            base_fr = horzcat_pad(base_fr, neuron_fr(baseline_idx, i) );
        end 
        
        % Calculate the baseline std for each neuron
        avg_base_fr = mean(base_fr, 1, 'omitnan');
        base_std = std(base_fr, 0, 1, 'omitnan');

        [row, col] = find(neuron_fr > 3*base_std + avg_base_fr);
        
        plot(timeline(row), spacing*.75 + col*spacing, '.b', 'MarkerSize', 6);
        xlim([-1 2.05]);
        Multi_func.set_default_axis(gca);
        
        title(f_stim, 'Interpreter', 'none');
    end

    % Save the figure
    saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_fr_base_comp.png']);
    %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_fr_base_comp.pdf']);
    %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_fr_base_comp.eps'], 'epsc');
end

% Plot each neuron's firing rate between baseline and stimulation
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
    
    for f_stim=stims'
        f_stim = f_stim{1};
        base_fr = [];
        stim_fr = [];

        for neuron = 1:length(data_bystim.(f_stim).all_trial_spikeidx)
            timeline = data_bystim.(f_stim).trace_timestamps(:, neuron); 
            
            % Find the region idxs
            base_idx = find(timeline < 0);
            stim_idx = find(timeline >= 0 & timeline <= 1);

            cur_nr = data_bystim.(f_stim).all_trial_spikeidx{neuron};
            
            base_raster = ismember(cur_nr, base_idx);
            stim_raster = ismember(cur_nr, stim_idx);

            % Calculate the spiking rate for 1 sec
            base_srate = sum(base_raster, 1)/1;
            stim_srate = sum(stim_raster, 1)/1;
           
            base_fr(end + 1) = mean(base_srate, 'omitnan');
            stim_fr(end + 1) = mean(stim_srate, 'omitnan');
        end

        % PLot the violins of the data
        figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
        data = [base_fr, stim_fr];
        %data = log2(data);
        labels = [repmat({'Base'}, 1, length(base_fr)), repmat({'Stim'}, 1, length(stim_fr))];
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);
        
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        hold on;
        
        Multi_func.set_default_axis(gca);
        set(gca, 'Units', 'centimeters', 'InnerPosition', [1, 1, 16.2/2, 20]);
        % Conditionally color lines between violinplots
        for i=1:length(base_fr)
            color = [0 0 0 0.4];
            if stim_fr(i) > base_fr(i)*1.10
                color = [[30, 2, 237]/255, 0.4];
            elseif stim_fr(i) < base_fr(i)*0.9
                color = [[235, 5, 28]/255, 0.4];
            end
        
            plot([1 2], ([base_fr(i), stim_fr(i)]), '-', 'Color', color);
            hold on;
        end
        ylabel('Firing Rate (Hz)');

        %Multi_func.set_default_axis(gca);
        title(['Firing Rate for ' f_stim(3:end)]);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.eps'], 'epsc');
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.pdf']);
    end
end
