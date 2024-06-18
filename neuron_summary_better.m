clear all;
close all;
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

% Get the first region field
field1 = fieldnames(region_data);
field1 = 'r_M1';
avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
display_names = 0;

%% Create the heatmaps for subthreshold Vm and spike raster
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Reset color_limits variable
    color_limits = [];    

    for f_stim = stims'
        f_stim = f_stim{1};
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';

        % Create a figure that includes: raw, spikes, and SubVm
        figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 40 40]);
        %tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'loose', 'Units', 'centimeters', 'InnerPosition', [4 5 8 8]);
        tiledlayout(1, 2, 'TileSpacing', 'none', 'Padding', 'loose', 'Units', 'centimeters', 'InnerPosition', [3 20 8.13 8.482]);
        
        % Plot the raw traces
        nexttile;
        
        % Calculate the sort order for Vm based on depolarization during stimulation phase
        nr_trans_dep = nanmean(data_bystim.(f_stim).neuron_Vm(avg_Fs:2*avg_Fs, :), 1);
        [~, I] = sort(nr_trans_dep);
        I = flip(I);

        neuron_bound = [0];
        data_map = [];
        for i = I %1:length(data_bystim.(f_stim).neuron_SubVm)
            % Non-linear normalized
            %cur_vm = data_bystim.(f_stim).all_trial_rawVm{i};
            %norm_vm = (cur_vm - min(cur_vm, [], 1))./(max(cur_vm, [], 1) - min(cur_vm, [], 1));
            % Normalized by spike amplitude
            cur_vm = data_bystim.(f_stim).all_trial_rawVm{i};
            norm_vm = cur_vm./data_bystim.(f_stim).neuron_spike_amp(i);

            data_map = [data_map; norm_vm'];
            num_neurons = size(data_bystim.(f_stim).all_trial_rawVm{i}, 2);
            neuron_bound(end + 1) = neuron_bound(end) + num_neurons;
        end
        imagesc('XData', timeline, 'YData', 1:size(data_map, 1), 'CData', data_map);

        a = colorbar;
        a.Label.String = 'Vm';

        % Set color limits if there was already
        if isempty(color_limits)
            color_limits = a.Limits;
        end
        caxis(color_limits);
        
        set(a, 'Location', 'westoutside')
        Multi_func.set_default_axis(gca);
        
        % Set neuron separation
        yline(neuron_bound + 0.5);
        neuron_mid = (neuron_bound(1:end - 1) + neuron_bound(2:end))/2;
        yticks(neuron_mid);
        yticklabels(1:length(neuron_mid))

        xlabel('Time from Stim onset(sec)');
        xlim([-.7 2.05]);
        ylim([0 size(data_map, 1)]);
        ylabel('Neuron (Black outlines trials of each neuron)');
        title('Normalized Traces', 'Interpreter', 'none');

        % Plot the raster plot
        nexttile;  
        neuron_bound = [0];
        index = 1;
        for fov = I %1:length(data_bystim.(f_stim).neuron_spikeidx)
            cur_color = [rand, rand, rand]*0.7;
            cur_fov = data_bystim.(f_stim).all_trial_spikeidx{fov};
            for tr = 1:size(cur_fov, 2)
                cur_spikeidx = cur_fov(:, tr);
                cur_spikeidx(isnan(cur_spikeidx)) = [];
                %plot(timeline(cur_spikeidx), repmat(index, length(cur_spikeidx), 1), '.k', 'MarkerSize', 3, 'color', cur_color);
                plot(timeline(cur_spikeidx), repmat(index, length(cur_spikeidx), 1), '.k', 'MarkerSize', 3);
                hold on;
                index = index + 1;
            end
            neuron_bound(end + 1) = index - 1;
            
            if display_names == 1
                % Plot the neuron label next to its raster
                text(2.2, index - (size(cur_fov, 2)/2), data_bystim.(f_stim).neuron_name{fov}, 'Interpreter', 'none');
                hold on;
            end
        end
        neuron_bound(end) = [];
        index = index - 1;
        yline(neuron_bound + 0.5);
        ax = gca;
        ax.YAxis.Visible = 'off';
        Multi_func.set_default_axis(gca);
        xlabel('Time from Stim onset(sec)');
        
        ylim([0 index]);
        set(gca, 'YTick', []);
        ylabel('Neuron Trials');

        % Adjust x axis if displaying names or not
        if display_names == 1
            xlim([-.7 4]);
        else
            xlim([-.7 2.05]);
        end

        title('Raster Plots', 'Interpreter', 'none');

        sgtitle([ f_region ' ' f_stim ' Neuronwise'], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.pdf']);
        savefig(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.fig']);
        
        %saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.eps'], 'epsc');
    end
end
