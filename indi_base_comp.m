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

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate
% Using a smaller window size
srate_win = 3;

% Do all region
all_region = 0;

% Percentage to increase Vm and firing rate
thres_percent = 0.10;

% Set threshold higher

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Load all of the data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
load(save_all_data_file);

% Combine all regions into 1
if all_region == 1
    region_data = Multi_func.combine_regions(region_data);
end

field1 = fieldnames(region_data);
field1 = field1{1};
avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%set figures off
set(0,'DefaultFigureVisible','off');

spacing = 1.4;

% Show each neuron's trial-averaged Vm
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
    
    for f_stim=stims'
        
        f_stim = f_stim{1};
        figure('visible', 'on', 'Renderer', 'Painters', 'Position', [200 200 2000 1000]);    
        tiledlayout(size(data_bystim.(f_stim).trace_timestamps, 2), 2, 'TileSpacing', 'none', 'Padding', 'compact');    

        
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2); 

        % Loop through each neuron
        for nr=1:size(data_bystim.(f_stim).trace_timestamps, 2)
            nexttile;
            % Plot each neuron's average Vm
            neuron_vm = data_bystim.(f_stim).neuron_Vm(:, nr);
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) < data_bystim.(f_stim).stim_timestamps(1, nr));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) >= data_bystim.(f_stim).stim_timestamps(1, nr) & data_bystim.(f_stim).trace_timestamps(:, nr) <= data_bystim.(f_stim).stim_timestamps(end, nr));

            % Plotting all neuron's Vm normalized
            %neuron_vm = (neuron_vm - min(neuron_vm, [], 1))./ (max(neuron_vm, [], 1) - min(neuron_vm, [], 1));
            %plot(timeline, neuron_vm + [spacing:spacing:size(neuron_vm, 2)*spacing], '-k');
            base_vm = mean(neuron_vm(baseline_idx), 'omitnan');
            stim_vm = mean(neuron_vm(stim_idx), 'omitnan');
            
            % Need to perform linear interpolation here
            xq = 1:length(stim_idx)/length(baseline_idx):length(stim_idx);
            downsample = interp1(1:length(stim_idx), neuron_vm(stim_idx), xq);
            
            % Check for distribution difference
            [p,h] = signrank(neuron_vm(baseline_idx), downsample);

            % Plot the 2 distributions
            histogram(neuron_vm(baseline_idx), 100, 'FaceColor', [0.8500 0.3250 0.0980]);
            hold on;
            histogram(downsample, 100, 'FaceColor', [0.4660 0.6740 0.1880]);
            legend('Baseline', ['Stimulation ' num2str(p)]);
            nexttile;

            % Check for Vm increase
            color = [0 0 0 0.4];
            if h == 1
                if stim_vm > base_vm
                    color = [[30, 2, 237]/255, 0.4];
                elseif stim_vm < base_vm
                    color = [[235, 5, 28]/255, 0.4];
                end
            end
            plot(timeline, neuron_vm, 'Color', color);
            hold on;
            
            %-- Find points with higher than 2 std baseline
            
            % Calculate the baseline std for each neuron
            %avg_base_vm = mean(base_vm, 1, 'omitnan');
            %base_std = std(base_vm, 0, 1, 'omitnan');
            %[row, col] = find(neuron_vm > 2*base_std + avg_base_vm);
            %plot(timeline(row), max(neuron_vm) + 1, '.b', 'MarkerSize', 6);
            %hold on;
            
            xlim([-1 2.05]);
            ylim([min(neuron_vm) - 1, max(neuron_vm) + 1]);
            set(gca,'XTick',[]);
            Multi_func.set_default_axis(gca);
            
            % Plot the DBS onset line
            xline(0, '--');
            hold on;
            xline(1, '--');

        end

        sgtitle([f_region(3:end) '_' f_stim(3:end) ' Vm' ], 'Interpreter', 'non');

        % Save the figure
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_vm_base_comp.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_vm_base_comp.pdf']);
        %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_vm_base_comp.eps'], 'epsc');
    end
end

% Show all of the neuron's trial-averaged Firing Rate
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
        
    for f_stim=stims'
        f_stim = f_stim{1};
            
        figure('visible', 'on', 'Renderer', 'Painters', 'Position', [200 200 2000 1000]);    
        tiledlayout(size(data_bystim.(f_stim).trace_timestamps, 2), 2, 'TileSpacing', 'none', 'Padding', 'compact');    
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2); 
        
        for nr=1:size(data_bystim.(f_stim).trace_timestamps, 2)
            nexttile;
            % Plot each neuron's average Firing Rate
            neuron_fr = data_bystim.(f_stim).neuron_srate_50(:, nr); 
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) < data_bystim.(f_stim).stim_timestamps(1, nr));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) >= data_bystim.(f_stim).stim_timestamps(1, nr) & data_bystim.(f_stim).trace_timestamps(:, nr) <= data_bystim.(f_stim).stim_timestamps(end, nr));

            % Plotting normalized firing rate
            %neuron_fr = (neuron_fr - min(neuron_fr, [], 1))./ (max(neuron_fr, [], 1) - min(neuron_fr, [], 1));
            %plot(timeline, neuron_fr + [spacing:spacing:size(neuron_fr, 2)*spacing], '-k');
            base_fr = mean(neuron_fr(baseline_idx), 'omitnan');
            stim_fr = mean(neuron_fr(stim_idx), 'omitnan');
            
            % Downsample the stimulation data with linear interpolation
            xq = 1:length(stim_idx)/length(baseline_idx):length(stim_idx);
            downsample = interp1(1:length(stim_idx), neuron_fr(stim_idx), xq);

            [p, h] = signrank(neuron_fr(baseline_idx), downsample);
            
            % Plot the 2 distributions
            histogram(neuron_fr(baseline_idx), 100, 'FaceColor', [0.8500 0.3250 0.0980]);
            hold on;
            histogram(downsample, 100, 'FaceColor', [0.4660 0.6740 0.1880]);
            legend('Baseline', ['Stimulation ' num2str(p)]);
            nexttile;

            % Plot the firing rate with color indicating increase or decrease
            if h == 1
                color = [0 0 0 0.4];
                if stim_fr > base_fr
                    color = [[30, 2, 237]/255, 0.4];
                elseif stim_fr < base_fr
                    color = [[235, 5, 28]/255, 0.4];
                end
            end
            plot(timeline, neuron_fr, 'Color', color);
            hold on;

            %-- Find points with higher than 2 std baseline
            % Calculate the baseline std for each neuron
            %avg_base_fr = mean(base_fr, 1, 'omitnan');
            %base_std = std(base_fr, 0, 1, 'omitnan');
            %[row, col] = find(neuron_fr > 2*base_std + avg_base_fr);
            %plot(timeline(row), max(neuron_fr) + 1, '.b', 'MarkerSize', 6);
            %hold on;
            
            xlim([-1 2.05]);
            ylim([min(neuron_fr) - 1, max(neuron_fr) + 1]);
            set(gca,'XTick',[]);
            Multi_func.set_default_axis(gca);

            % Plot the DBS onset line
            xline(0, '--');
            hold on;
            xline(1, '--');
        end

        sgtitle([f_region(3:end) ' ' f_stim(3:end) ' Firing Rate' ], 'Interpreter', 'none');

        % Save the figure
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp.pdf']);
        %saveas(gcf, [figure_path 'Neuronwise/' f_region '_indi_fr_base_comp.eps'], 'epsc');
    end
end

% Violins of each neuron's Vm between baseline and stimulation
for f_region = fieldnames(region_data)'    
    f_region = f_region{1};    
    data_bystim = region_data.(f_region);    
    stims = fieldnames(data_bystim);    
    
    for f_stim=stims'
        f_stim = f_stim{1};
        base_vm = [];
        stim_vm = [];

        for neuron = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            timeline = data_bystim.(f_stim).trace_timestamps(:, neuron); 
            
            % Find the region idxs
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, neuron) < data_bystim.(f_stim).stim_timestamps(1, neuron));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, neuron) >= data_bystim.(f_stim).stim_timestamps(1, neuron) & data_bystim.(f_stim).trace_timestamps(:, neuron) <= data_bystim.(f_stim).stim_timestamps(end, neuron));

            neuron_vm = data_bystim.(f_stim).neuron_Vm(:, neuron);
           
            base_vm(end + 1) = mean(neuron_vm(baseline_idx), 'omitnan');
            stim_vm(end + 1) = mean(neuron_vm(stim_idx), 'omitnan');
        end

        % PLot the violins of the data
        figure(  'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
        data = [base_vm, stim_vm];
        %data = log2(data);
        labels = [repmat({'Base'}, 1, length(base_vm)), repmat({'Stim'}, 1, length(stim_vm))];
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);
        
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        hold on;
        
        Multi_func.set_default_axis(gca);
        set(gca, 'Units', 'centimeters', 'InnerPosition', [1, 1, 16.2/2, 20]);
        % Conditionally color lines between violinplots
        for i=1:length(base_vm)
            color = [0 0 0 0.4];
            if stim_vm(i) > base_vm(i)*(1 + thres_percent)
                color = [[30, 2, 237]/255, 0.4];
            elseif stim_vm(i) < base_vm(i)*(1 - thres_percent)
                color = [[235, 5, 28]/255, 0.4];
            end
        
            plot([1 2], ([base_vm(i), stim_vm(i)]), '-', 'Color', color);
            hold on;
        end
        ylabel('Vm (A.U.)');

        %Multi_func.set_default_axis(gca);
        title([f_region(3:end) ' Vm for ' f_stim(3:end)]);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_vm_base_comp_violin.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_vm_base_comp_violin.eps'], 'epsc');
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_vm_base_comp_violin.pdf']);
    end
end


% Violis of each neuron's firing rate between baseline and stimulation
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
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, neuron) < data_bystim.(f_stim).stim_timestamps(1, neuron));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, neuron) >= data_bystim.(f_stim).stim_timestamps(1, neuron) & data_bystim.(f_stim).trace_timestamps(:, neuron) <= data_bystim.(f_stim).stim_timestamps(end, neuron));

            cur_nr = data_bystim.(f_stim).all_trial_spikeidx{neuron};
            
            base_raster = ismember(cur_nr, baseline_idx);
            stim_raster = ismember(cur_nr, stim_idx);

            % Calculate the spiking rate for 1 sec
            base_srate = sum(base_raster, 1)/1;
            stim_srate = sum(stim_raster, 1)/1;
           
            base_fr(end + 1) = mean(base_srate, 'omitnan');
            stim_fr(end + 1) = mean(stim_srate, 'omitnan');
        end

        % PLot the violins of the data
        figure(  'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
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
            if stim_fr(i) > base_fr(i)*(1 + thres_percent)
                color = [[30, 2, 237]/255, 0.4];
            elseif stim_fr(i) < base_fr(i)*(1 - thres_percent)
                color = [[235, 5, 28]/255, 0.4];
            end
        
            plot([1 2], ([base_fr(i), stim_fr(i)]), '-', 'Color', color);
            hold on;
        end
        ylabel('Firing Rate (Hz)');

        %Multi_func.set_default_axis(gca);
        title([f_region(3:end) ' Firing Rate for ' f_stim(3:end)]);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.eps'], 'epsc');
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_indi_fr_base_comp_violin.pdf']);
    end
end
