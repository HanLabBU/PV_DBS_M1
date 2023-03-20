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

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Figures' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate
srate_win = 50;

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);

all_matfiles = {ses.name};

% Select matfiles by brain region
[region_matfiles] = Multi_func.find_region(all_matfiles);
region_data = struct();
all_Fs = [];
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    %% Select matfiles by stim specific conditions for all regions
    %[matfile_stim] = stim_cond(all_matfiles); 
    %% Select matfiles by stim condition for given region
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);

    %% Loop through each field of the struct and concatenate everything together
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies

    % Loop through each stimulation condition
    for f_stim = fieldnames(matfile_stim)'
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    
    
        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();
        data_bystim.(f_stim).neuron_Vm = [];
        data_bystim.(f_stim).neuron_srate = [];
        data_bystim.(f_stim).neuron_spec_power = [];
        data_bystim.(f_stim).neuron_spec_freq = [];
        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_subVm = [];
            cur_fov_srate = [];
            cur_fov_raster = [];
            cur_fov_stim_time = [];
            cur_fov_trace_time = [];

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    %Determine whether this roi is to be ignored for this particular trial
                    ri = strsplit(matfile{1}, '_');
                    try
                        trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                    catch
                        trial_ignr_list = [];
                    end
    
                    % Check if current trial is in the ignore list
                    if ismember(tr_idx, trial_ignr_list)
                        continue;
                    end
                    
                    % If the trial data is empty, that means it was skipped
                    if isempty(trial_data)
                        continue;
                    end

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;

                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');
                    
                    % Calculate the spike rate
                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win, 1, 1);
                    cur_fov_srate = horzcat_pad(cur_fov_srate, cur_spikerate');            
                    
                    % Store the raster plot
                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
                    % DEBUG
                    if max(cur_fov_trace_time(:, end)) < 2
                        stim_trace_sz = size(cur_fov_stim_time)
                        trace_time_sz = size(cur_fov_trace_time)
                        pause;
                        matfile{1}
                    end
                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_Vm;
            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
            % Store average spike rate for each neuron
            temp = data_bystim.(f_stim).neuron_srate;
            data_bystim.(f_stim).neuron_srate = horzcat_pad(temp, nanmean(cur_fov_srate, 2));
            % Store the timestamp data
            temp = data_bystim.(f_stim).stim_timestamps;
            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
            temp = data_bystim.(f_stim).trace_timestamps;
            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
            % Calculate and save frequency data
            [wt, f] = get_power_spec(nanmean(cur_fov_subVm, 2)', nanmean(cur_fov_Fs));
            temp = data_bystim.(f_stim).neuron_spec_power;
            data_bystim.(f_stim).neuron_spec_power = cat(3, temp, wt);
            temp = data_bystim.(f_stim).neuron_spec_freq;
            data_bystim.(f_stim).neuron_spec_freq = cat(3, temp, f);

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

%% Region separated analysis
% Plot subthreshold Vm by frequency stimulation and brain stimulation
%for f_region = fieldnames(region_data)'
%    f_region = f_region{1};
%    data_bystim = region_data.(f_region).data_bystim;
%
%    %TODO find the stimulation onset with the timestamps from the stimulation and camera
%    stims = fieldnames(data_bystim);
%    avg_Fs = nanmean(all_Fs);
%    % Calculate the trial timeline, convert the idx's into timestamps
%    % Adding 4 here to account for the timestamps dropped by the alignment and motion correction
%    timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;
%    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%    for stim=stims'
%        nexttile;
%        plot(timeline, nanmean(data_bystim.(stim{1}).neuron_Vm, 2));
%        title(stim, 'Interpreter', 'none');
%    end
%    sgtitle(['Average Sub Vm for ' num2str(f_region)], 'Interpreter', 'none');
%end
%
%% Plot the subthreshold frequency spectrum for each region
%for f_region = fieldnames(region_data)'
%    f_region = f_region{1};
%    data_bystim = region_data.(f_region).data_bystim;
%
%    %TODO find the stimulation onset with the timestamps from the stimulation and camera
%    stims = fieldnames(data_bystim);
%    avg_Fs = nanmean(all_Fs);
%    % Calculate the trial timeline, convert the idx's into timestamps
%    % Adding 4 here to account for the timestamps dropped by the alignment and motion correction
%    timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;
%    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%    for stim=stims'
%        nexttile;
%        
%        % All frequency weights
%        all_wt = [];
%        freqLimits = [0 250];
%
%        % Loop through each neuron subthreshold trace
%        for i = 1:size(data_bystim.(stim{1}).neuron_Vm, 2)
%            cur_subVm = data_bystim.(stim{1}).neuron_Vm(:, i);
%            if isnan(cur_subVm)
%                continue;
%            end
%            fb = cwtfilterbank(SignalLength=length(cur_subVm),...
%                               SamplingFrequency=avg_Fs,...
%                               FrequencyLimits=freqLimits);
%            [wt, f] = cwt(cur_subVm, FilterBank=fb);
%            all_wt = cat(3, all_wt, wt);
%        end
%            
%        contourf(timeline, f, nanmean(abs(wt), 3), 'edgecolor', 'none');
%        colorbar;
%        title(stim, 'Interpreter', 'none');
%    end
%    sgtitle(['Average Sub Vm for ' num2str(f_region)], 'Interpreter', 'none');
%end

avg_Fs = nanmean(all_Fs);
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

%%% All PV neurons summary
%%Raw subthreshold spectra
%stims = fieldnames(data_bystim);
%figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%    surface(timeline, nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3), nanmean(abs(data_bystim.(f_stim{1}).neuron_spec_power), 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%    colorbar;
%    title(f_stim{1}(3:end), 'Interpreter', 'none');
%end
%sgtitle('Spectra from averaged Sub Vm Trace');
%
%% Subthreshold spectra with baseline period ratio-normalization for each neuron and 
%% Averaged afterwards
%stims = fieldnames(data_bystim);
%figure('visible', 'off','Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%   
%    % Normalize ratio
%    % First calculate by doing the ratios for each neuron first and then averaging
%    cur_spec_pow = data_bystim.(f_stim{1}).neuron_spec_power;
%
%    % Plot the starting time point for each neuron
%    sz = size(data_bystim.(f_stim{1}).trace_timestamps);
%
%    % Loop throug each neuron
%    for i = 1:size(data_bystim.(f_stim{1}).trace_timestamps, 2)
%        baseline_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) < data_bystim.(f_stim{1}).stim_timestamps(1, i));
%        %plot(baseline_idx, i*5, '.r');
%        %hold on;
%        base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
%        cur_spec_pow(:, :, i) = abs(cur_spec_pow(:, :, i))./base_power;
%    end
%    surface(nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)', ... 
%            nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3), ...
%            nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%    colorbar;
%    set(gca, 'color', 'none');
%    %avg_power = nanmean(data_bystim.(f_stim{1}).neuron_spec_power, 3);
%
%    title(f_stim{1}(3:end), 'Interpreter', 'none');
%end
%sgtitle('Spectra with baseline-ratio normalized individually Sub Vm Trace');
%
%% Plot the average spectra for baseline, stimulation, and offset periods
%stims = fieldnames(data_bystim);
%figure('visible', 'off','Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%   
%    % Normalize ratio
%    % First calculate by doing the ratios for each neuron first and then averaging
%    cur_spec_pow = data_bystim.(f_stim{1}).neuron_spec_power;
%
%    % Plot the starting time point for each neuron
%    sz = size(data_bystim.(f_stim{1}).trace_timestamps);
%    
%    % Store the baseline power spectra values 
%    all_base_power = [];
%    all_stim_power = [];
%    all_offset_power = [];
%
%    % Loop throug each neuron
%    for i = 1:size(data_bystim.(f_stim{1}).trace_timestamps, 2)
%        baseline_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) < data_bystim.(f_stim{1}).stim_timestamps(1, i));
%        stim_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) >= data_bystim.(f_stim{1}).stim_timestamps(1, i) & ...
%                        data_bystim.(f_stim{1}).trace_timestamps(:, i) <= data_bystim.(f_stim{1}).stim_timestamps(end, i));
%                    %TODO I think the offset idx here is actually wrong?
%        offset_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) > data_bystim.(f_stim{1}).stim_timestamps(end, i));
%        
%        base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
%        all_base_power = horzcat(all_base_power, base_power);
%        stim_power = nanmean(abs(cur_spec_pow(:, stim_idx, i)), 2);
%        all_stim_power = horzcat(all_stim_power, stim_power);
%        offset_power = nanmean(abs(cur_spec_pow(:, offset_idx, i)), 2);
%        all_offset_power = horzcat(all_offset_power, offset_power);
%    end
%
%    % Get the average frequency 
%    spec_freq = nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3);
% 
%    % TODO keep working here
%    %plot(nanmean());
%    set(gca, 'color', 'none');
%    %avg_power = nanmean(data_bystim.(f_stim{1}).neuron_spec_power, 3);
%
%    title(f_stim{1}(3:end), 'Interpreter', 'none');
%end
%sgtitle('');


% Subthreshold spectra with (x - A)/(A + B) normalization for each neuron and 
% Averaged afterwards
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        nexttile;
       
        % Normalize ratio
        % First calculate by doing the ratios for each neuron first and then averaging
        cur_spec_pow = data_bystim.(f_stim{1}).neuron_spec_power;
    
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim{1}).trace_timestamps);
    
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim{1}).trace_timestamps, 2)
            % Get the indices define the periods within in the trace
            baseline_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) < data_bystim.(f_stim{1}).stim_timestamps(1, i));
            stim_idx = find(data_bystim.(f_stim{1}).trace_timestamps(:, i) >= data_bystim.(f_stim{1}).stim_timestamps(1, i) & ...
                            data_bystim.(f_stim{1}).trace_timestamps(:, i) <= data_bystim.(f_stim{1}).stim_timestamps(end, i));
            %plot(baseline_idx, i*5, '.r');
            %hold on;
            base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
            stim_power = nanmean(abs(cur_spec_pow(:, stim_idx, i)), 2);
            cur_spec_pow(:, :, i) = (abs(cur_spec_pow(:, :, i)) - base_power)./(base_power + stim_power);
        end
        surface(nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)', ... 
                nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3), ...
                nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        a = colorbar;
        a.Label.String = 'Power (A.U.)';
    
        set(gca, 'color', 'none');
        %avg_power = nanmean(data_bystim.(f_stim{1}).neuron_spec_power, 3);
        xlabel('Time from Stim onset(sec)');
        ylabel('Freq (Hz)');
        title(f_stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([ f_region ' Spectra with (x - A)/(A + B) normalization individually Sub Vm Trace'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Spectra/' f_region 'A_B_Normalization_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region 'A_B_Normalization_Spectra.eps'], 'epsc');
end

% Peform normalization after average subthreshold
%figure('visible', 'off','Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%   
%    % Get the average spectra
%    cur_spec_pow = nanmean(abs(data_bystim.(f_stim{1}).neuron_spec_power), 3);
%    % Average the stimulation start time
%    cur_stim_start = nanmean(data_bystim.(f_stim{1}).stim_timestamps(1, :));
%    % Average trace timestamps
%    cur_trace_time = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
%
%    % Grab the baseline spec power
%    baseline_idx = find(cur_trace_time < cur_stim_start);
%    base_power = nanmean(abs(cur_spec_pow(:, baseline_idx)), 2);
%    cur_spec_pow = abs(cur_spec_pow)./base_power;
%
%    surface(cur_trace_time', ... 
%            nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3), ...
%            cur_spec_pow, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%    colorbar;
%    %avg_power = nanmean(data_bystim.(f_stim{1}).neuron_spec_power, 3);
%
%    title(f_stim{1}(3:end), 'Interpreter', 'none');
%end
%sgtitle('Spectra with baseline-ratio normalized from the averages Sub Vm Trace');
%
%% Raw subthreshold Vm spectra
%figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
%tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%% Loop through each stimulation parameter
%for f_stim=stims'
%    nexttile;
%    surface(timeline, nanmean(data_bystim.(f_stim{1}).neuron_spec_freq, 3), nanmean(abs(data_bystim.(f_stim{1}).neuron_spec_power), 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%    colorbar;
%    title(f_stim{1}(3:end), 'Interpreter', 'none');
%end
%sgtitle('Spectra from averaged Sub Vm Trace');


%% Full collective spike rate over time
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2);
        cur_srate = mean(data_bystim.(stim{1}).neuron_srate, 2, 'omitnan');
        std_srate = std(data_bystim.(stim{1}).neuron_srate, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(stim{1}).neuron_srate, 2);
        %num_points = size(data_bystim.(stim{1}).neuron_srate, 1);
        sem_srate = cur_srate./sqrt(num_neurons);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_srate + sem_srate; flipud(cur_srate - sem_srate)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_srate, 'k', 'LineWidth', 1);
        set(gca, 'color', 'none')
        title(stim{1}(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Average Spike rate'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Summary_Continuous_FiringRate.eps']);
end

%DEBUG
return;

%% Subthreshold Vm
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    for stim=stims'
        timeline = nanmean(data_bystim.(f_stim{1}).trace_timestamps, 2)';
        cur_Vm = mean(data_bystim.(stim{1}).neuron_Vm, 2, 'omitnan');
        std_Vm = std(data_bystim.(stim{1}).neuron_Vm, 0, 2, 'omitnan');
        num_neurons = size(data_bystim.(stim{1}).neuron_Vm, 2);
        sem_Vm = std_Vm./sqrt(num_neurons);
        num_points = size(data_bystim.(stim{1}).neuron_Vm, 1);
        nexttile;
        f = fill([timeline; flip(timeline)], [cur_Vm + sem_Vm; flipud(cur_Vm - sem_Vm)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(timeline, cur_Vm, 'k', 'LineWidth', 1);
        title(stim{1}(3:end), 'Interpreter', 'none');
        set(gca, 'color', 'none')
    end
    sgtitle([f_region ' Average subthreshold Vm'], 'Interpreter', 'none');
    saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.png']);
    saveas(gcf, [figure_path 'Average/' f_region '_Average_sub_thres.eps']);
end



% Calculate cwt for input signal and 
function [wt, f] = get_power_spec(signal, samp_freq)
    freqLimits = [0 150];
    fb = cwtfilterbank(SignalLength=length(signal),...
                       SamplingFrequency=samp_freq,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(signal, FilterBank=fb);
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
