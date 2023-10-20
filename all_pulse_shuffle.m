close all;
f = filesep;


%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';


% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate

% Specify seed
rng(100);

% Do all region
all_region = 0;
set(0,'DefaultFigureVisible','on');

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

% Get the timeline and average framerate
field1 = fieldnames(region_data);
field1 = field1{1};
avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;


% All pulse for Vm
% Just look at M1 data
for f_region = {'r_M1'} %fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters' ); %,'InnerPosition', [4, 20, 3.62, 5.16]);
    for f_stim=stims'
        f_stim = f_stim{1};
        all_pulse_vm = [];
        extra_trace = 3;
        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim).neuron_Vm, 2)
            
            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim).stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(data_bystim.(f_stim).trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width./nr_avg_trace_time);
 
            % Use fixed length of points for all frequencies
            %trace_width = 21;

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim).stim_timestamps(:, nr)'
                follow_trace_idx = find(pulse_time <= data_bystim.(f_stim).trace_timestamps(:, nr));
                
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                
                % Calculate end trace idx using the average inter-pulse time
                %end_trace_idx = find(pulse_time + trace_width >= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                % Using fixed number of pulse times
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;

                Vm_pulse_width = data_bystim.(f_stim).neuron_Vm(start_trace_idx:end_trace_idx, nr);
                all_pulse_vm = horzcat_pad(all_pulse_vm, Vm_pulse_width);
            end
        end
        fig_all_pulse_avg = mean(all_pulse_vm, 2, 'omitnan');

        % Shuffle and create null distribution of vm between pulses from each pulse window
        
        % Shuffle x amount of times
        wind_size = size(all_pulse_vm, 1);
        shuf_val_dist = [];
        for i=1:1000
            shuf_all_pulse_vm = [];
            for j=1:size(all_pulse_vm, 2)
                shuf_idx = randperm(wind_size);
                shuf_all_pulse_vm = horzcat_pad(shuf_all_pulse_vm, all_pulse_vm(shuf_idx, j));
            end

            % Calculate the average value
            shuf_avg_all_pulse = mean(shuf_all_pulse_vm, 2, 'omitnan');
            
            % This was the original value
            %shuf_val_dist(end + 1) = mean(avg_all_pulse, 'omitnan');;
            shuf_val_dist = horzcat_pad(shuf_val_dist, shuf_avg_all_pulse);
        end
        nexttile;
        histogram(shuf_val_dist(:), 100);
        hold on;
        low_perc = prctile(shuf_val_dist(:), 2.5);
        high_perc = prctile(shuf_val_dist(:), 97.5);
        xline(low_perc, 'b', 'DisplayName', ['Low Perc: ' num2str(low_perc)]);
        hold on;
        xline(high_perc, 'b', 'DisplayName', ['High Perc: ' num2str(high_perc)]);
        hold on;
        xline(max(fig_all_pulse_avg), 'k', 'DisplayName', ['Observ ' num2str(max(fig_all_pulse_avg))]);
        
        disp([f_region(3:end) ' ' f_stim(3:end) ' Vm']);
        disp(['Lower precent ' num2str(low_perc)]);
        disp(['higher precent ' num2str(high_perc)]);
        disp(['Max val is ' num2str(max(fig_all_pulse_avg))]);
        fprintf('\n\n');
    
        legend();
        title([ f_stim(3:end) ' Vm All Pulse Shuffled'], 'Interpreter', 'none');
    
    end
    sgtitle([f_region(3:end) ' Vm All Pulse Avg Shuffled'], 'Interpreter', 'none');
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.png']);
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.pdf']);
end

% Shuffling the firing rate 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters' ); %, 'InnerPosition', [4, 20, 3.5, 5]);
    for f_stim=stims'
        f_stim = f_stim{1};
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
        
        all_pulse_fr = [];
        extra_trace = 3;
        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim).neuron_spikecounts_raster, 2)

            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim).stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(data_bystim.(f_stim).trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width./nr_avg_trace_time);
                
            %Neuron framerate
            nr_rate = data_bystim.(f_stim).framerate(nr);

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim).stim_timestamps(:, nr)'
                follow_trace_idx = find(pulse_time <= data_bystim.(f_stim).trace_timestamps(:, nr));
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;
                %end_trace_idx = find(pulse_time + nr_avg_pulse_width >= data_bystim.(f_stim).trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                fr_pulse_width = data_bystim.(f_stim).neuron_spikecounts_raster(start_trace_idx:end_trace_idx, nr)*nr_rate;
                
                all_pulse_fr = horzcat_pad(all_pulse_fr, fr_pulse_width);
            end
        end
        fig_all_pulse_avg = mean(all_pulse_fr, 2, 'omitnan');

        % Shuffle x amount of times
        wind_size = size(all_pulse_fr, 1);
        shuf_val_dist = [];
        for i=1:1000 
            shuf_all_pulse_fr = [];
            for j=1:size(all_pulse_fr, 2)
                shuf_idx = randperm(wind_size);
                shuf_all_pulse_fr = horzcat_pad(shuf_all_pulse_fr, all_pulse_fr(shuf_idx, j));
            end

            % Calculate the average value
            shuf_avg_all_pulse = mean(shuf_all_pulse_fr, 2, 'omitnan');
            
            % This was the original value
            %shuf_val_dist(end + 1) = mean(avg_all_pulse, 'omitnan');;
            shuf_val_dist = horzcat_pad(shuf_val_dist, shuf_avg_all_pulse);
        end
        nexttile;
        histogram(shuf_val_dist(:), 100);
        hold on;
        low_perc = prctile(shuf_val_dist(:), 2.5);
        high_perc = prctile(shuf_val_dist(:), 97.5);
        xline(low_perc, 'b', 'DisplayName', ['Low Perc: ' num2str(low_perc)]);
        hold on;
        xline(high_perc, 'b', 'DisplayName', ['High Perc: ' num2str(high_perc)]);
        hold on;
        xline(max(fig_all_pulse_avg), 'k', 'DisplayName', ['Observ ' num2str(max(fig_all_pulse_avg))]);
        
        disp([f_region(3:end) ' ' f_stim(3:end) ' firing rate' ]);
        disp(['Lower precent ' num2str(low_perc)]);
        disp(['higher precent ' num2str(high_perc)]);
        disp(['Max val is ' num2str(max(fig_all_pulse_avg))]);
        fprintf('\n\n');

        %cur_srate = mean(all_pulse_fr, 2, 'omitnan');
        legend();
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Firing Rate all pulse average'], 'Interpreter', 'none');
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.png']);
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.eps'], 'epsc');
end
