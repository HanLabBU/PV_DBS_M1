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

% Parameters for frames to chop ofccfront_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
%pv_data_path = Multi_func.save_plot;

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine whether to combine all regions as one data
all_regions = 0;
%%% END Modification

%TODO add Post stimulation distribution period as well

% Check that the server path exists
%if ~isfolder(local_root_path)
%    disp('Server rootpath does not exist!!!');
%    return;
%end

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
%Load the data
load(save_all_data_file);


% Check if combining all of the regions or not
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

% grab the average framerate
field1 = fieldnames(region_data);
field1 = field1(1);
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');
%timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

% Look at spike phase Vm at stimulation frequency

%% Plotting spike-Vm phase filtered at stim frequency
stimfreq_phase_stats = struct();
stats_log = [figure_path 'Phase' f 'Stim_freq_stats'];
if exist(stats_log), delete(stats_log), end;
diary(stats_log);
diary off
% Loop through all regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    % Look at spike phase at Vm stimulation frequency
    % Go through and average each spike phase
    stims = fieldnames(data_bystim);
    for f_stim=stims'
        f_stim = f_stim{1};
        stim_num = str2num(f_stim(3:end));
        
        figure('Renderer', 'Painters', 'Position', [200 200 2000 700])
        tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 9.38, 4.94]);
    
        base_phases = [];
        stim_phases = [];
    
        %Loop through each neuron 
        for nr = 1:length(data_bystim.(f_stim).neuron_hilbfilt)
            base_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) < data_bystim.(f_stim).stim_timestamps(1, nr));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) >= data_bystim.(f_stim).stim_timestamps(1, nr) & ...
                            data_bystim.(f_stim).trace_timestamps(:, nr) <= data_bystim.(f_stim).stim_timestamps(end, nr));
    
            hilb_nr = data_bystim.(f_stim).neuron_hilbfilt{nr};
    
            % Loop through each trial
            for tr = 1:size(hilb_nr, 3)
                % Grab the spike idx from current trial
                cur_tr_spikeidx = data_bystim.(f_stim).all_trial_spikeidx{nr};
                cur_tr_spikeidx = cur_tr_spikeidx(:, tr);
    
                % Get baseline spikes
                base_spikeidx = intersect(cur_tr_spikeidx, base_idx);
                % Get stimulation spikes
                stim_spikeidx = intersect(cur_tr_spikeidx, stim_idx);
                
                % Grab the phases of each trial period
                base_phases = [base_phases, angle(hilb_nr(stim_num, base_spikeidx, tr))];
                stim_phases = [stim_phases, angle(hilb_nr(stim_num, stim_spikeidx, tr))];
                
            end
        end
    
        % Save circular staistics
        diary on
        fprintf('\n\n');
        disp(['Stim Frequency Phase Circular Stats' f_region ' ' f_stim]);
        stimfreq_phase_stats.(f_region).(f_stim).p = circ_wwtest(base_phases(:), stim_phases(:));
        diary off

        % Plot the polar histgrams
        nexttile;
        edges = linspace(0, 2*pi, 24);
        %polarhistogram(base_phases, edges, 'Normalization', 'probability', 'FaceColor', Multi_func.base_color, 'FaceAlpha', 0.3);
        polarhistogram(base_phases, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.base_color);
        title('Base');
        set(gca, 'Color', 'none');

        nexttile;
        %polarhistogram(stim_phases, edges, 'Normalization', 'probability', 'FaceColor', Multi_func.stim_color, 'FaceAlpha', 0.3);
        polarhistogram(stim_phases, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.stim_color);
        title('Stim');
    
        set(gca, 'Color', 'none');
        sgtitle([f_region(3:end) ' Filtered at Stim Frequency ' f_stim(3:end) ' Hz'], 'Interpreter', 'none');
     
        saveas(gcf, [figure_path 'Phase/' f_region '_' f_stim '_Spike_Phase_StimFreq.png']);
        saveas(gcf, [figure_path 'Phase/' f_region '_' f_stim '_Spike_Phase_StimFreq.pdf']);
        %saveas(gcf, [figure_path 'Phase/' f_stim '_Spike_Phase_StimFreq.eps']);
    end
end


%% Plotting phase for 2-10Hz filtered
lowfreq_phase_stats = struct();
region_data.(f_region).(f_stim).LowFreq = struct();
stats_log = [figure_path 'Phase' f 'Low_freq_stats'];
if exist(stats_log), delete(stats_log), end;
diary(stats_log);
diary off
% Loop through all regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    % Look at spike phase at Vm stimulation frequency
    % Go through and average each spike phase
    stims = fieldnames(data_bystim);

    for f_stim=stims'
        f_stim = f_stim{1};

        % Specify low frequency range to filter
        low_range = Multi_func.theta_range; % Should be the same as in Multi Func[2, 10]; % Delta/Theta frequency
        
        figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
        tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 11.50, 6.20]);
    
        base_phases = [];
        stim_phases = [];
    
        %Loop through each neuron 
        for nr = 1:length(data_bystim.(f_stim).all_trial_spikeidx)
            base_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) < data_bystim.(f_stim).stim_timestamps(1, nr));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) >= data_bystim.(f_stim).stim_timestamps(1, nr) & ...
                            data_bystim.(f_stim).trace_timestamps(:, nr) <= data_bystim.(f_stim).stim_timestamps(end, nr));
    

            % Loop through each trial
            for tr = 1:size(data_bystim.(f_stim).all_trial_SubVm{nr}, 2)
                % Grab the spike idx from current trial
                cur_tr_spikeidx = data_bystim.(f_stim).all_trial_spikeidx{nr};
                cur_tr_spikeidx = cur_tr_spikeidx(:, tr);
    
                % Get baseline spikes
                base_spikeidx = intersect(cur_tr_spikeidx, base_idx);
                % Get stimulation spikes
                stim_spikeidx = intersect(cur_tr_spikeidx, stim_idx);
    
                % Get the hilbert of the filtered 2-10Hz subVm
                cur_hilb = Multi_func.filt_range(data_bystim.(f_stim).all_trial_SubVm{nr}(:, tr), low_range, avg_Fs)';

                % Grab the phases of each period
                base_phases = [base_phases, angle(cur_hilb(base_spikeidx))];
                stim_phases = [stim_phases, angle(cur_hilb(stim_spikeidx))];
                
            end
        end
        
        % Add theta phases to region_data struct
        region_data.(f_region).(f_stim).LowFreq.base_phases = base_phases;
        region_data.(f_region).(f_stim).LowFreq.base_phases = stim_phases;

        % Save circular statistics
        diary on;
        disp(['Low Frequency Phase Circular Stats' f_region ' ' f_stim]);
        lowfreq_phase_stats.(f_region).(f_stim).p = circ_wwtest(base_phases(:), stim_phases(:))
        diary off;

        % Plot the polar histgrams
        nexttile;
        edges = linspace(0, 2*pi, 24);
        %polarhistogram(base_phases, edges, 'Normalization', 'probability', 'FaceColor', Multi_func.base_color, 'FaceAlpha', 0.3);
        polarhistogram(base_phases, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.base_color);
        rlim([0 .15]);
        title('Base');
        set(gca, 'Color', 'none');
    
        nexttile;
        %polarhistogram(stim_phases, edges, 'Normalization', 'probability', 'FaceColor', Multi_func.stim_color, 'FaceAlpha', 0.3);
        polarhistogram(stim_phases, edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.stim_color);
        rlim([0 .15]);
        title('Stim');
        set(gca, 'Color', 'none');
    
        sgtitle([f_region(3:end) ' 2-10Hz Delta/Theta Frequency for stim ' f_stim(3:end)], 'Interpreter', 'none');
     
        saveas(gcf, [figure_path 'Phase/' f_region '_' f_stim '_Spike_Phase_deltaTheta.png']);
        saveas(gcf, [figure_path 'Phase/' f_region '_' f_stim '_Spike_Phase_deltaTheta.pdf']);
        %saveas(gcf, [figure_path 'Phase/' f_stim '_Spike_Phase_deltaTheta.eps']);
    end
end

%% Plot individual neurons low frequency power during baseline period
low_freqs = [1:50] %Multi_func.theta_frequencies;
% Loop through all regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region); %

    % Look at spike phase at Vm stimulation frequency
    % Go through and average each spike phase
    stims = fieldnames(data_bystim);

    for f_stim=stims'
        f_stim = f_stim{1};
        
        % Save the power and relative neuron names
        base_all_power_norm = [];
        power_norm_neuron_name = {};
        
        % This is a parallel array that will indicate whether is has low frequency stuff or not
        has_low_freq = [];

        % Loop through each neuron
        for nr=1:size(data_bystim.(f_stim).neuron_Vm, 2)
            base_idx = find(data_bystim.(f_stim).trace_timestamps(:, nr) < data_bystim.(f_stim).stim_timestamps(1, nr));
            nr_filt_data = data_bystim.(f_stim).neuron_hilbfilt{nr};
            neuron_data = data_bystim.(f_stim);

            base_power = [];
            vm_map = [];
            % Loop through each trial
            for tr=1:size(nr_filt_data, 3)
                trace_noise = neuron_data.all_trial_trace_noise{nr}(tr);

                % Add the trace with SBR as value
                vm_map = [vm_map; neuron_data.all_trial_rawVm{nr}(:, tr)'./trace_noise];
     
                base_power(:, end + 1) = nanmean(abs(nr_filt_data(low_freqs, base_idx, tr)), 2);
            end
            
            % Calculate the average power
            avg_power = nanmean(base_power, 2);
            std_power = std(base_power, 0, 2, 'omitnan');
            num_trials = size(base_power, 2);
            sem_power = std_power./sqrt(num_trials);
            
            timeline = neuron_data.trace_timestamps(:, nr);
            freqs = 1:size(base_power, 1);

            % z-score the power across frequency
            base_all_power_norm(:, end + 1) = zscore(avg_power(:, end), [], 1);
            power_norm_neuron_name{end + 1} = neuron_data.neuron_name{nr};

            % Determine if neuron has low frequency oscillations
            has_low_freq(end + 1) = nanmean(base_all_power_norm(Multi_func.theta_frequencies, end), 1) > 1;

            % Make a figure and plot all of the poewrs
            %figure('Position', [0 0 2000 1800]);
            %tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

            %nexttile;
            %imagesc('XData', timeline, 'YData', 1:size(vm_map, 1), 'CData', vm_map);
            %
            %nexttile;

            %fill_h = fill([freqs'; flip(freqs)'], [avg_power + sem_power; flipud(avg_power - sem_power)], [0.5, 0.5, 0.5]);
            %hold on;

            %plot(freqs, avg_power' );
            %xticks(0:10:150);
            %ylabel('Power');
            %xlabel('Frequency (Hz)');


            % %Plotting the zscored power across frequency
            %nexttile;
            %plot(freqs, base_all_power_norm(:, end));

            %sgtitle(neuron_data.neuron_name{nr}, 'Interpreter', 'none');
            %saveas(gcf, [figure_path 'Spectra' f 'Individual' f ...
            %    neuron_data.neuron_name{nr} '_power_spec_density.png']);
            
            % End Making figure

        end
        region_data.(f_region).(f_stim).has_low_freq = has_low_freq;
        has_low_freq
        
        % Theta frequency variable: Multi_func.theta_frequencies
        low_freq_power_norm = nanmean(base_all_power_norm(Multi_func.theta_frequencies, :), 1);
        [~, I] = sort(low_freq_power_norm);
        low_freq_power_norm = low_freq_power_norm(I);
        power_norm_neuron_name = power_norm_neuron_name(I);

        % Plot normalized base power
        figure('Position', [0 0 1800 1800]);
        plot(low_freq_power_norm, '.');
        hold on;
        yline(1);
        ylabel('2-10Hz zscored power');
        xlabel('Neuron Number');
        xticks(1:length(power_norm_neuron_name));
        xticklabels(power_norm_neuron_name);
        set(gca,'TickLabelInterpreter','none');
        title([f_region ' ' f_stim ' 2-10Hz zscored sorted'], 'Interpreter', 'none');
    end
end


%% Create bar plots of the broad-band spike phase PLV



%DEBUG
%nexttile;
%Fn = avg_Fs/2;
%FB = [140*0.8 140*1.2];
%[B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
%
%spike_idxs = [base_spikeidx; stim_spikeidx];

%filt_sig = filtfilt(B, A, data_bystim.(f_stim).all_trial_SubVm{nr}(:, tr));
%plot(timeline, filt_sig);
%hold on;
%plot(timeline(spike_idxs), filt_sig(spike_idxs), '|r');
%hold on;
%
%for i = 1:length(spike_idxs)
%    %TODO plot the phase right above the spike indicator
%    text(timeline(spike_idxs(i)), filt_sig(spike_idxs(i)) + 0.5, num2str(rad2deg(angle(hilb_nr(stim_num, spike_idxs(i), tr))) + 360 ));
%    hold on;
%end    
