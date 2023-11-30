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

% Do all region
all_region = 0;


set(0,'DefaultFigureVisible','off');

%%% END Modification

% Seed number for random number generation
rng(100);

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

%ses = dir([pv_data_path '*.mat']);
%
%all_matfiles = {ses.name};
%
%% Select matfiles by brain region
%[region_matfiles] = Multi_func.find_region(all_matfiles);
%region_data = struct();
%all_Fs = [];
%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%    % Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    % Select matfiles by stim condition for given region
%    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);
%
%    % Loop through each field of the struct and concatenate everything together
%    % Store trace aspect data by each individual stimulation condition
%    data_bystim = struct();
%    % Store all of the calculated sampling frequencies
%
%    % Loop through each stimulation condition
%    for f_stim = fieldnames(matfile_stim)'
%        f_stim = f_stim{1};
%        matfiles = matfile_stim.(f_stim).names;    
%    
%        % Initialize field subthreshold array
%        data_bystim.(f_stim) = struct();
%        data_bystim.(f_stim).neuron_Vm = [];
%        data_bystim.(f_stim).neuron_srate = [];
%        data_bystim.(f_stim).stim_timestamps = [];
%        data_bystim.(f_stim).trace_timestamps = [];
%
%        % Loop through each matfile of the current stimulation condition
%        for matfile = matfiles
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile{1}]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%            cur_fov_Fs = [];
%            cur_fov_subVm = [];
%            cur_fov_srate = [];
%            cur_fov_raster = [];
%            cur_fov_stim_time = [];
%            cur_fov_trace_time = [];
%
%            % Loop through each ROI
%            for roi_idx=1:size(trial_data.detrend_traces, 2)
%                %Determine whether this roi is to be ignored for this particular trial
%                ri = strsplit(matfile{1}, '_');
%                try
%                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
%                catch
%                    trial_ignr_list = [];
%                end
%
%                % Remove trials from trial idx list
%                trial_idxs = setdiff(trial_idxs, trial_ignr_list);
%
%                % Skip ROI if there are at most 2 trials
%                if length(trial_idxs) <= 2
%                    continue;
%                end
%
%                % Loop through each trial                
%                for tr_idx=trial_idxs        
%                    trial_data = data.align.trial{tr_idx};
%                    raw_trial_data = data.raw.trial{tr_idx};
%
%                    % Store the camera framerate
%                    all_Fs(end+1) = trial_data.camera_framerate;
%                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;
%
%                    % Grab the subthreshold Vm
%                    % Chop the respective frames
%                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
%                    [baseline, coeff] = Multi_func.exp_fit_F(cur_trace_ws', round(trial_data.camera_framerate));
%                    detrend_subVm = cur_trace_ws - baseline;
%                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');
%                    
%                    % Calculate the spike rate
%                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
%                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
%                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win, 1, 1);
%                    cur_fov_srate = horzcat_pad(cur_fov_srate, cur_spikerate');            
%                    
%                    % Store the raster plot
%                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');
%
%                    % Store all of the timestamp info
%                    stim_start = raw_trial_data.raw_stimulation_time(1);
%                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, ...
%                                                    raw_trial_data.raw_stimulation_time(1:str2num(ri{5}) ) - stim_start);
%                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, ...
%                                                    trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
%                    % DEBUG
%                    if max(cur_fov_trace_time(:, end)) < 2
%                        stim_trace_sz = size(cur_fov_stim_time)
%                        trace_time_sz = size(cur_fov_trace_time)
%                        pause;
%                        matfile{1}
%                    end
%                end % End looping through each neuron
%            end
%            
%            % Skip rest of the calculations if the subthreshold Vm is nan
%            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
%                continue;
%            end
%
%            % Average for each neuron and save the subthreshold Vm
%            temp = data_bystim.(f_stim).neuron_Vm;
%            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
%            % Store average spike rate for each neuron
%            temp = data_bystim.(f_stim).neuron_srate;
%            data_bystim.(f_stim).neuron_srate = horzcat_pad(temp, nanmean(cur_fov_srate, 2));
%            % Store the timestamp data
%            temp = data_bystim.(f_stim).stim_timestamps;
%            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
%            temp = data_bystim.(f_stim).trace_timestamps;
%            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
%
%        end % End looping through FOVs of a condition
%    end
%
%    % Save the VM to the specific region
%    region_data.(f_region).data_bystim = data_bystim;
%end

% Load all of the data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
load(save_all_data_file);

% Combine all regions into 1
if all_region == 1
    region_data = Multi_func.combine_regions(region_data);
end

% Get the first region field
field1 = fieldnames(region_data);
field1 = field1{1};
avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%set figures off
set(0,'DefaultFigureVisible','off');

%% Spike rate showing first few pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)*1000;
        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate_3, 2, 'omitnan');
        std_srate = std(data_bystim.(f_stim{1}).neuron_srate_3, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate_3, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate_3, 1);
        sem_srate = std_srate./sqrt(num_neurons);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;

        % Plot the DBS stimulation time pulses
        xline(nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the semi-significant points on the subthreshold Vm
        base_srate = [];
        % Grab the baseline sub srate for all neurons
        for i = 1:size(data_bystim.(f_stim{1}).trace_timestamps, 2)
            baseline_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) < data_bystim.(f_stim{1}).stim_timestamps(1, i));
            base_srate = horzcat_pad(base_srate, data_bystim.(f_stim{1}).neuron_srate_3(baseline_idx, i) );
        end
        base_srate = mean(base_srate, 2, 'omitnan');
        std_baseline = std(base_srate, 0, 'omitnan');
        sig_idx = find(cur_srate > (2*std_baseline + mean(base_srate, 'omitnan')));
            
        % Calculate a good height for the points
        thres = 0.2;
        height = range(cur_srate);
        height = (1+thres)*height;
        %plot(timeline(sig_idx), repmat(height, 1, length(sig_idx)), '.b', 'MarkerSize', 10);
        
        % Plotting the 2 std line 
        yline(2*std_baseline, '--', 'Color', [0 0 0 0.5]);

        % Increase timescale resolution
        xlim([0 - 50, 0 + 100]);
        x = gca; x = x.XAxis;
        Multi_func.set_spacing_axis(x, 50, 1);
        Multi_func.set_default_axis(gca);
        ylabel('Firing Rate(Hz)');
        xlabel('Time from onset(ms)');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Population Spike rate with first pulse'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_FR.pdf']);
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_FR.eps'], 'epsc');
end

%% Sub Vm averaged across all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 2, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.62*2, 5.16]);
    for f_stim=stims'
        
        Vm_avg = [];
        all_pulse_vm = [];
        all_trans_pulse_vm = [];
        trans_vm_avg = [];
        all_sus_pulse_vm = [];
        sus_vm_avg = [];
        extra_trace = 3;
        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim{1}).neuron_Vm, 2)
            
            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim{1}).stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(data_bystim.(f_stim{1}).trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width./nr_avg_trace_time);
 
            % Use fixed length of points for all frequencies
            %trace_width = 21;

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim{1}).stim_timestamps(:, nr)'
                
                follow_trace_idx = find(pulse_time <= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                
                % Calculate end trace idx using the average inter-pulse time
                %end_trace_idx = find(pulse_time + trace_width >= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                % Using fixed number of pulse times
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;

                Vm_pulse_width = data_bystim.(f_stim{1}).neuron_Vm(start_trace_idx:end_trace_idx, nr);
                Vm_avg = horzcat_pad(Vm_avg, Vm_pulse_width);

                % Store all of the Vm pulses
                all_pulse_vm = horzcat_pad(all_pulse_vm, Vm_pulse_width);

                % Store all of the transient Vm pulses
                if pulse_time < Multi_func.trans_ped(2)/1000
                    all_trans_pulse_vm = horzcat_pad(all_trans_pulse_vm, Vm_pulse_width);
                    trans_vm_avg = horzcat_pad(trans_vm_avg, Vm_pulse_width);
                end

                % Store all of the sustained Vm pulses
                if pulse_time >= Multi_func.sus_ped(1)/1000 && pulse_time <= Multi_func.sus_ped(2)/1000
                    all_sus_pulse_vm = horzcat_pad(all_sus_pulse_vm, Vm_pulse_width);
                    sus_vm_avg = horzcat_pad(sus_vm_avg, Vm_pulse_width);
                end
            end
        end


        % Calculate the whole stimulation period average
        cur_subVm = mean(Vm_avg, 2, 'omitnan');
        std_subVm = std(Vm_avg, 0, 2, 'omitnan');
        num_pulses = size(Vm_avg, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_subVm, 1);
        sem_subVm = std_subVm./sqrt(num_pulses);
        nexttile;
        
        % Use the average trace framerate to calculate the timelinefigure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
        timeline = [ [0:size(Vm_avg, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        f = fill([timeline, flip(timeline)], [cur_subVm + sem_subVm; flipud(cur_subVm - sem_subVm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_subVm, 'k', 'LineWidth', 1);
        hold on;

        % Shuffle all pulse average and plot average lines
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

        % Calculate percentile values
        low_perc = prctile(shuf_val_dist(:), 2.5);
        high_perc = prctile(shuf_val_dist(:), 97.5);
        shuf_mean = mean(shuf_val_dist(:), 'omitnan');
        
        % Plot the shuffled values as dashed horizontal lines
        %yline([low_perc, high_perc], '--');
        %hold on;
        %yline(shuf_mean);

        % Plot the percentiles as a different colored shading
        shade_yvals = [repmat(high_perc, 1, length(timeline)), repmat(low_perc, 1, length(timeline))];
        f = fill([timeline, flip(timeline)], shade_yvals, [0.62 0.71 1]);
        Multi_func.set_fill_properties(f);
        hold on;
        yline(shuf_mean, '--');

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width*1000:nr_avg_pulse_width*1000], 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = 0;
        posy = 0;
        timelength = 0.5;
        %plot([posx, posx + timelength], [posy posy], 'k', 'LineWidth', 2);
        %text(posx, posy - 0.2, [num2str(timelength*(1/avg_Fs)) 'ms']);

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Vm');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);

        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title(f_stim{1}(3:end), 'Interpreter', 'none');
        %-- End of plotting the pulse averages for whole stim time


        %-- Calculate the transient stimulation period average
        cur_subVm = mean(trans_vm_avg, 2, 'omitnan');
        std_subVm = std(trans_vm_avg, 0, 2, 'omitnan');
        num_pulses = size(trans_vm_avg, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_subVm, 1);
        sem_subVm = std_subVm./sqrt(num_pulses);
        nexttile;
        
        % Use the average trace framerate to calculate the timelinefigure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
        timeline = [ [0:size(trans_vm_avg, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        f = fill([timeline, flip(timeline)], [cur_subVm + sem_subVm; flipud(cur_subVm - sem_subVm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_subVm, 'k', 'LineWidth', 1);
        hold on;

        % Shuffle all pulse average and plot average lines
        wind_size = size(all_trans_pulse_vm, 1);
        shuf_val_dist = [];
        for i=1:1000
            shuf_all_trans_pulse_vm = [];
            for j=1:size(all_trans_pulse_vm, 2)
                shuf_idx = randperm(wind_size);
                shuf_all_trans_pulse_vm = horzcat_pad(shuf_all_trans_pulse_vm, all_trans_pulse_vm(shuf_idx, j));
            end

            % Calculate the average value
            shuf_avg_all_pulse = mean(shuf_all_trans_pulse_vm, 2, 'omitnan');
            
            % This was the original value
            %shuf_val_dist(end + 1) = mean(avg_all_pulse, 'omitnan');;
            shuf_val_dist = horzcat_pad(shuf_val_dist, shuf_avg_all_pulse);
        end

        % Calculate percentile values
        low_perc = prctile(shuf_val_dist(:), 2.5);
        high_perc = prctile(shuf_val_dist(:), 97.5);
        shuf_mean = mean(shuf_val_dist(:), 'omitnan');
        
        % Plot the shuffled values as dashed horizontal lines
        %yline([low_perc, high_perc], '--');
        %hold on;
        %yline(shuf_mean);

        % Plot the percentiles as a different colored shading
        shade_yvals = [repmat(high_perc, 1, length(timeline)), repmat(low_perc, 1, length(timeline))];
        f = fill([timeline, flip(timeline)], shade_yvals, [0.62 0.71 1]);
        Multi_func.set_fill_properties(f);
        hold on;
        yline(shuf_mean, '--');

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width*1000:nr_avg_pulse_width*1000], 'Color', [170, 176, 97]./255, 'LineWidth', 2);
        hold on;

        % Plot the timescale bar
        posx = 0;
        posy = 0;
        timelength = 0.5;
        %plot([posx, posx + timelength], [posy posy], 'k', 'LineWidth', 2);
        %text(posx, posy - 0.2, [num2str(timelength*(1/avg_Fs)) 'ms']);

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        Multi_func.set_default_axis(gca);
        ylabel('Vm');
        xlabel('Time from pulse(ms)');
        %ylim([0, 7]);

        % Remove x-axis and right y-axis
        %set(gca,'xtick',[]);
        title([f_stim{1}(3:end) ' Transient'], 'Interpreter', 'none');
        %-- End of transient plotting the pulse averages

    end
    sgtitle([f_region(3:end) ' Sub Vm all pulse average'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.pdf']);
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_Vm.eps'], 'epsc');
end

%% Spike rate averaged across all pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    for f_stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        
        FR_avg = [];
        extra_trace = 3;
        all_pulse_fr = [];

        % Looping through each neuron
        for nr = 1:size(data_bystim.(f_stim{1}).neuron_spikecounts_raster, 2)

            % Calculate the number of trace idxs between pulses
            nr_avg_pulse_width = mean(diff(data_bystim.(f_stim{1}).stim_timestamps(:, nr) ), 'omitnan');
            nr_avg_trace_time = mean(diff(data_bystim.(f_stim{1}).trace_timestamps(:, nr) ), 'omitnan');
            trace_width = ceil(nr_avg_pulse_width./nr_avg_trace_time);
                
            %Neuron framerate
            nr_rate = data_bystim.(f_stim{1}).framerate(nr);

            % Loop through each stimulation time pulses
            for pulse_time = data_bystim.(f_stim{1}).stim_timestamps(:, nr)'
                follow_trace_idx = find(pulse_time <= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                start_trace_idx = follow_trace_idx(1) - extra_trace;
                end_trace_idx = follow_trace_idx(1) + trace_width + extra_trace;
                %end_trace_idx = find(pulse_time + nr_avg_pulse_width >= data_bystim.(f_stim{1}).trace_timestamps(:, nr));
                %end_trace_idx = end_trace_idx(end);
                
                fr_pulse_width = data_bystim.(f_stim{1}).neuron_spikecounts_raster(start_trace_idx:end_trace_idx, nr)*nr_rate;
                
                FR_avg = horzcat_pad(FR_avg, fr_pulse_width);
                all_pulse_fr = horzcat_pad(all_pulse_fr, fr_pulse_width);
            end
        end

        cur_srate = mean(FR_avg, 2, 'omitnan');
        std_srate = std(FR_avg, 0, 2, 'omitnan');
        num_pulses = size(FR_avg, 2);
        %num_points = size(data_bystim.(f_stim{1}).neuron_srate_20, 1);
        sem_srate = std_srate./sqrt(num_pulses);
        
        timeline = [ [0:size(FR_avg, 1) - 1] - extra_trace]*1000./avg_Fs;
        
        nexttile;
        f = fill([timeline, flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 0.3);
        hold on;

        % Perform shuffling for all firing rate
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

        % Calculate percentile values
        low_perc = prctile(shuf_val_dist(:), 2.5);
        high_perc = prctile(shuf_val_dist(:), 97.5);
        shuf_mean = mean(shuf_val_dist(:), 'omitnan');

        % Plot the shuffled values
        % Plotting using dashed lines
        %yline([low_perc, high_perc], '--');
        %hold on;
        %yline(shuf_mean);
        shade_yvals = [repmat(high_perc, 1, length(timeline)), repmat(low_perc, 1, length(timeline))];
        f = fill([timeline, flip(timeline)], shade_yvals, [0.62 0.71 1]);
        Multi_func.set_fill_properties(f);
        hold on;
        yline(shuf_mean, '--');

        % Plot the DBS stimulation time pulses
        xline([0:nr_avg_pulse_width*1000:nr_avg_pulse_width*1000], 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the timescale bar
        %posx = 1;
        %posy = 0;
        %plot([posx, posx + 4], [posy posy], 'k', 'LineWidth', 2);
        %text(posx, posy - 0.1, [num2str(4) 'ms']);

        % Increase timescale resolution
        %xlim([0 - .100, 0 + .100]);
        x = gca; x = x.XAxis;
        %Multi_func.set_spacing_axis(x, 6, 1);
        Multi_func.set_default_axis(gca);
        ylabel('Firing Rate(Hz)');
        xlabel('Time from pulse(ms)');
        
        %ylim([-1 10]);

        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Firing Rate all pulse average'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_All_Pulse_Avg_FR.eps'], 'epsc');
end

%% Spike rate showing display all pulses
%for f_region = fieldnames(region_data)'
%    f_region = f_region{1};
%    data_bystim = region_data.(f_region);
%    stims = fieldnames(data_bystim);
%    
%    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%    for f_stim=stims'
%        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
%        cur_srate = mean(data_bystim.(f_stim{1}).neuron_srate_20, 2, 'omitnan');
%        std_srate = std(data_bystim.(f_stim{1}).neuron_srate_20, 0, 2, 'omitnan');
%        num_neurons = size(data_bystim.(f_stim{1}).neuron_srate_20, 2);
%        %num_points = size(data_bystim.(f_stim{1}).neuron_srate, 1);
%        sem_srate = std_srate./sqrt(num_neurons);
%        nexttile;
%        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
%        Multi_func.set_fill_properties(f);
%        hold on;
%        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
%        hold on;
%
%        % Plot the DBS stimulation time pulses
%        stim_time = nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2);
%        xline(stim_time, 'Color', [170, 176, 97]./255, 'LineWidth', 2);
%        hold on;
%
%        % Plot the timescale bar
%        posx = -.100;
%        posy = 0;
%        plot([posx, posx + 0.050], [posy posy], 'k', 'LineWidth', 2);
%        text(posx, posy - 0.2, '50ms');
%
%        % Increase timescale resolution
%        xlim([0 - .10, max(stim_time) + .10]);
%        a = gca;
%        a.XAxis.Visible = 'off';
%        set(gca, 'color', 'none');
%        ylabel('Firing Rate(Hz)');
%        title(f_stim{1}(3:end), 'Interpreter', 'none');
%    end
%    sgtitle([f_region '_' num2str(srate_win)  ' Population Spike rate with all pulse'], 'Interpreter', 'none');
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.png']);
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.pdf']);
%    saveas(gcf, [figure_path 'Average/' f_region '_Display_All_Pulse_Trig_FR.eps'], 'epsc');
%end

%% Subthreshold Vm first few pulses
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.40, 4.96]);
    for f_stim=stims'
        
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)*1000; % Convert to ms
        cur_Vm = mean(data_bystim.(f_stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(f_stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(f_stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        %num_points = size(data_bystim.(f_stim{1}).neuron_Vm, 1);
        nexttile;

        % Standard Error
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        
        % Plot Vm
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        hold on;
        
        % Plot the DBS stimulation time pulses
        xline(nanmean(data_bystim.(f_stim{1}).stim_timestamps, 2)*1000, 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);
        hold on;

        % Plot the semi-significant points on the subthreshold Vm

        base_Vm = [];
        % Grab the baseline sub Vm for all neurons
        for i = 1:size(data_bystim.(f_stim{1}).trace_timestamps, 2)
            baseline_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) < data_bystim.(f_stim{1}).stim_timestamps(1, i));
            base_Vm = horzcat_pad(base_Vm, data_bystim.(f_stim{1}).neuron_Vm(baseline_idx, i) );
        end
        
        % Calculate each Vm point's significance based on the population average of all neurons
        base_Vm = mean(base_Vm, 2, 'omitnan');
        std_baseline = std(base_Vm, 0, 'omitnan'); 
        sig_idx = find(cur_Vm > (2*std_baseline + mean(base_Vm, 'omitnan')));
            
        % Calculate a good height for marking significance
        thres = 0.2;
        height = range(cur_Vm);
        height = (1+thres)*height;
        %plot(timeline(sig_idx), repmat(height, 1, length(sig_idx)), '.b', 'MarkerSize', 8);
        
        % Plotting the 2 std line 
        yline(2*std_baseline, '--', 'Color', [0 0 0 0.5]);

        % Increase timescale resolution
        xlim([0 - 50, 0 + 100]);
        ylabel('Vm');
        xlabel('Time from onset(ms)');

        Multi_func.set_default_axis(gca);
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region(3:end) ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Vm.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Vm.pdf']);
    %saveas(gcf, [figure_path 'Average/' f_region '_First_Pulse_Vm.eps'], 'epsc');
end

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
