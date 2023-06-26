clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine how much before and after a stimulation pulse to take the average
trace_sur = 10; % This is ~6ms before and after stim pulse

all_regions = 1;

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
%Load the data
load(save_all_data_file);

if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

avg_Fs = mean(region_data.r_combine.f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;


for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 2000 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

    for f_stim=stims'
        f_stim = f_stim{1};
        
        % Store neuron's Vm for heatmap
        Vm_map = [];
        time_to_peak = [];

        % Find the peaks in the Vm
        avg_Vm = mean(data_bystim.(f_stim).neuron_Vm, 2, 'omitnan');
        nexttile;
        
        [pks, locs, w, p] = findpeaks(avg_Vm);
    
        % Plot the average Vm and each of the peak locations
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        plot(timeline, avg_Vm);
        hold on;
        plot(timeline(locs), avg_Vm(locs), 'or');
        hold on;
        
        % Plot the prominences from the peak value to the bottom
        for i=1:length(locs)
            plot([timeline(locs(i)), timeline(locs(i)) ], [avg_Vm(locs(i)) avg_Vm(locs(i)) - p(i)], '-g');
            hold on;
        end
        title(['Stim ' f_stim], 'Interpreter', 'none');

        % find the max peak value and the distance from 0

        % Get points within the stimulation period
        stim_idxs = find(timeline >= 0 & timeline <= 1);
        [c, ia] = ismember(locs, stim_idxs);
        locs = locs(c);
        pks = pks(c);
        p = p(c);
        prom_peak_loc = find(max(p) == p);
        plot([0 timeline(locs(prom_peak_loc)) ], [0 0 ]);

        disp(['Population Average Time from ' f_stim]);
        disp(num2str(timeline(locs(prom_peak_loc))*1000 ));

    end

    sgtitle('Population average peak');
end

% Calculate the peak by averaging the peak to each neuron's Vm, 
% This was not really working well because the peak for some neurons was way later
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    for f_stim=stims'
        f_stim = f_stim{1};
        
        
        %figure('Renderer', 'Painters', 'Position', [200 200 2000 700]);
        %tiledlayout(length(data_bystim.(f_stim).neuron_Vm, 2), 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Store neuron's Vm for heatmap
        Vm_map = [];
        time_to_peak = [];

        for nr = 1:size(data_bystim.(f_stim).neuron_Vm, 2)
            
            Vm_map(end + 1, :) = data_bystim.(f_stim).neuron_Vm(:, nr)';

            stim_time = data_bystim.(f_stim).stim_timestamps(:, nr);
            timeline = data_bystim.(f_stim).trace_timestamps(:, nr);
            
            cur_vm = data_bystim.(f_stim).neuron_Vm(:, nr);

            % Find the max prominescences 
            [pks, locs, w, p] = findpeaks(cur_vm);
            
            % Look at first 100ms of the stimulation period
            stim_idxs = find(timeline >= 0 & timeline <= 0.1);
            [c, ia] = ismember(locs, stim_idxs);
            locs = locs(c);
            pks = pks(c);
            p = p(c);
            prom_peak_loc = find(max(p) == p);

            time_to_peak(end + 1) = timeline(locs(prom_peak_loc));
 
            
            % DEBUG find the max point for each trace
            % Show the max Vm after stimulation 
            %nexttile;
            %figure('Renderer', 'Painters', 'Position', [500 500 2000 700]);
            %plot(timeline, cur_vm, '-k');
            %hold on;
            %plot(timeline(locs), cur_vm(locs), 'or');
            %hold on;
            %plot([timeline(locs)'; timeline(locs)'], [cur_vm(locs)'; cur_vm(locs)' - p'], '-g');
            %hold on;
            %plot([0 timeline(locs(prom_peak_loc)) ], [0 0 ], '-m');

        end
        
        % Print the histogram of the timepoints
        figure('Renderer', 'Painters', 'Position', [500 500 2000 700]);
        histogram(time_to_peak*1000, 0:0.2:100);
        title(f_stim, 'Interpreter', 'none');

        % Print out the time results
        disp([f_stim ' time to peak Vm']);
        mean_peak = mean(time_to_peak*1000, 'omitnan');
        std_peak = std(time_to_peak*1000, 'omitnan');
        disp([num2str(mean_peak) '±' num2str(std_peak)]);
    end
end

return;

%% Plot the histograms for all pulse to spike times
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 500]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    histogram(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    histogram(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:)*1000, 0:0.2:27);
    hold on;
    avg_first_pts = 1000*data_bystim.(f_stim{1}).all_first_pulse_spike_times(:);
    xline(mean(avg_first_pts, 'omitnan'), 'o' );
    hold on;
    avg_all_pts = 1000*data_bystim.(f_stim{1}).all_all_pulse_spike_times(:);
    xline(mean(avg_all_pts, 'omitnan'), 'b' );
    xlabel('Time (ms)');
    legend('All Pulses', 'First Pulses', 'Avg first', 'Avg all');
    title([f_stim{1}(3:end)], 'Interpreter', 'none');
    
    % Print the descriptive statistics for the spike to pulse stuff
    disp(['Printing stats for ' f_stim{1}]);
    disp(['Mean±std: ' num2str(mean(avg_all_pts, 'omitnan')) '±' num2str(std(avg_all_pts, 'omitnan'))]);
    disp('');
end
sgtitle('Spike to Pulse Times');
saveas(gcf, [figure_path 'Inter_Spike' f 'Hist_Spike_to_pulse.png']);
saveas(gcf, [figure_path 'Inter_Spike' f 'Hist_Spike_to_pulse.eps'], 'epsc');

%% Plot the violing plots pulse to spike times
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 500]);
tiledlayout(1, length(stims), 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through each stimulation parameter
for f_stim=stims'
    nexttile;
    data = [data_bystim.(f_stim{1}).all_all_pulse_spike_times(:)*1000; data_bystim.(f_stim{1}).all_first_pulse_spike_times(:)*1000];
    num_all_pulses = length(data_bystim.(f_stim{1}).all_all_pulse_spike_times(:));
    num_first_pulses = length(data_bystim.(f_stim{1}).all_first_pulse_spike_times(:));
    labels = [repmat({'All pulses'}, num_all_pulses, 1); repmat({'First pulses'}, num_first_pulses, 1)];
    violinplot2(data, labels);
    ylabel('Time (ms)');
    title([f_stim{1}(3:end)], 'Interpreter', 'none');
    set(gca, 'color', 'none');

    % Print the descriptive statistics for the spike to pulse stuff
    disp(['Printing stats for ' f_stim{1}]);
    disp(['Mean±std: ' num2str(mean(avg_all_pts, 'omitnan')) '±' num2str(std(avg_all_pts, 'omitnan'))]);
    disp('');
end
sgtitle('Spike to Pulse Times');
saveas(gcf, [figure_path 'Inter_Spike' f 'Violin_Spike_to_pulse.png']);
saveas(gcf, [figure_path 'Inter_Spike' f 'Violin_Spike_to_pulse.eps'], 'epsc');

%% Specific functions for determining which FOVs to look at
% Return matfiles by stimulation condition
function [cond_struct] = stim_cond(matfile_names)
    cond_struct = struct();
    
    % Loop through each matfilename and group by stimulation conditions
    for i=1:length(matfile_names)
            file_parts = split(matfile_names{i}, '_');
            stim = file_parts{5};
            
            % Create stimulation field if it does not exist
            if ~isfield(cond_struct, ['f_' stim])
                cond_struct.(['f_' stim]).names = {};
            end

            cond_struct.(['f_' stim]).names{end+1} = matfile_names{i};
    end
end
