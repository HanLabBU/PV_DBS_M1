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
srate_win = 100;

% Parameter to determine whether to combine all regions as one data
all_regions = 1;

% Define transient and sustained periods
trans_ped = [0, 150];
sus_ped = [150, 1000];
offset_trans_ped = [1000, 1150];

%%% END Modification

% TODO need to clean up this script file to use the save data matfile with the modified baseline period

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
%    %% Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    %% Select matfiles by stim condition for given region
%    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);
%
%    %% Loop through each field of the struct and concatenate everything together
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
%        data_bystim.(f_stim).neuron_spec_power = [];
%        data_bystim.(f_stim).neuron_spec_freq = [];
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
%                % Remove ignored trials from trial idx list
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
%                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
%                    detrend_subVm = cur_trace_ws - baseline;
%                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');
%
%                    % Store all of the timestamp info
%                    stim_start = raw_trial_data.raw_stimulation_time(1);
%                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
%                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
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
%            % Store the timestamp data
%            temp = data_bystim.(f_stim).stim_timestamps;
%            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
%            temp = data_bystim.(f_stim).trace_timestamps;
%            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
%            
%            % Calculate and save frequency data
%            [wt, f] = get_power_spec(nanmean(cur_fov_subVm, 2)', nanmean(cur_fov_Fs));
%            temp = data_bystim.(f_stim).neuron_spec_power;
%            data_bystim.(f_stim).neuron_spec_power = cat(3, temp, wt);
%            temp = data_bystim.(f_stim).neuron_spec_freq;   
%            data_bystim.(f_stim).neuron_spec_freq = cat(3, temp, f);
%
%        end % End looping through FOVs of a condition
%    end
%
%    % Save the VM to the specific region
%    region_data.(f_region).data_bystim = data_bystim;
%end

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
%Load the data
load(save_all_data_file);

% Check if combining all of the regions or not
%if all_regions == 1
%    region_data = Multi_func.combine_regions_old(region_data);
%end

if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

% Calculate the sampling frequency from all of the 
avg_Fs = mean(region_data.r_combine.f_40.framerate, 'omitnan');

% Subthreshold time series spectra with (x - A)/(A + B) normalization for each neuron and 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 900 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
        
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        
        % Normalize ratio
        % First calculate by doing the ratios for each neuron first and then averaging
        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
    
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
    
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            % Get the indices define the periods within in the trace
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= data_bystim.(f_stim).stim_timestamps(1, i) & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= data_bystim.(f_stim).stim_timestamps(end, i));
            %plot(baseline_idx, i*5, '.r');
            %hold on;
            base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
            stim_power = nanmean(abs(cur_spec_pow(:, stim_idx, i)), 2);
            cur_spec_pow(:, :, i) = (abs(cur_spec_pow(:, :, i)) - base_power)./(base_power + stim_power);
        end

        surface(timeline, ... 
                nanmean(data_bystim.(f_stim).neuron_spec_freq, 3), ...
                nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        % Add DBS bar
        hold on;
        Multi_func.plot_dbs_bar([0, 1], 155, f_stim);

        a = colorbar;
        a.Ticks = linspace(a.Limits(1), a.Limits(2), 4);
        a.TickLabels = num2cell(round(linspace(a.Limits(1), a.Limits(2), 4), 1));
        a.Label.String = 'Relative Power';
        Multi_func.set_default_axis(gca);
        xlabel('Time from Stim onset(sec)');
        xlim([-1 2.05]);
        ylabel('Freq (Hz)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([ f_region ' Time Series Spectra with (x - A)/(A + B) normalization individually'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_Time_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_Time_Spectra.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_Time_Spectra.eps'], 'epsc');
    savefig(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_Time_Spectra.fig']);
end

% Subthreshold time series spectra with (x - A)/(A + B) normalization for each neuron and 
% zoom in on the 50 Hz
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('Renderer', 'Painters', 'Position', [200 200 650 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
        
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        
        % Normalize ratio
        % First calculate by doing the ratios for each neuron first and then averaging
        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
    
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
    
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            % Get the indices define the periods within in the trace
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= data_bystim.(f_stim).stim_timestamps(1, i) & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= data_bystim.(f_stim).stim_timestamps(end, i));
            %plot(baseline_idx, i*5, '.r');
            %hold on;
            base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
            stim_power = nanmean(abs(cur_spec_pow(:, stim_idx, i)), 2);
            cur_spec_pow(:, :, i) = (abs(cur_spec_pow(:, :, i)) - base_power)./(base_power + stim_power);
        end

        surface(timeline, ... 
                nanmean(data_bystim.(f_stim).neuron_spec_freq, 3), ...
                nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        % Add DBS bar
        hold on;
        Multi_func.plot_dbs_bar([0, 1], 155, f_stim);

        a = colorbar;
        a.Ticks = linspace(a.Limits(1), a.Limits(2), 4);
        a.TickLabels = num2cell(round(linspace(a.Limits(1), a.Limits(2), 4), 1));
        a.Label.String = 'Relative Power';
        Multi_func.set_default_axis(gca);
        xlabel('Time from Stim onset(sec)');
        xlim([-1 2.05]);
        ylabel('Freq (Hz)');
        ylim([0, 50]);
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([ f_region ' Time Series Spectra with (x - A)/(A + B) normalization individually'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Norm_Zoom50_Time_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Norm_Zoom50_Time_Spectra.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Norm_Zoom50_Time_Spectra.eps'], 'epsc');
end

%% Raw sub Vm spectra
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 900 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        surface(timeline, nanmean(data_bystim.(f_stim).neuron_spec_freq, 3), nanmean(abs(data_bystim.(f_stim).neuron_spec_power), 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        a = colorbar;
        a.Label.String = 'Power (A.U.)';
        hold on;
        Multi_func.plot_dbs_bar([0, 1], 155, f_stim);
        Multi_func.set_default_axis(gca);
        xlabel('Time from Stim onset(sec)');
        xlim([-1 2.05]);
        ylabel('Freq (Hz)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Raw Sub Vm Spectra'], 'Interpreter', 'none');

    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_Spectra.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_Spectra.eps'], 'epsc');
end

%% Vm Spectra zscored across time
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
    
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';

        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
    
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            cur_spec_pow(:, :, i) = zscore(abs(cur_spec_pow(:, :, i)), [], 2);
        end

        surface(timeline, nanmean(data_bystim.(f_stim).neuron_spec_freq, 3), nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        a = colorbar;
        a.Label.String = 'Power (A.U.)';
    
        set(gca, 'color', 'none');
        %avg_power = nanmean(data_bystim.(f_stim).neuron_spec_power, 3);
        xlabel('Time from Stim onset(sec)');
        ylabel('Freq (Hz)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Sub Vm Spectra Z-scored time'], 'Interpreter', 'none');

    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_time_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_time_Spectra.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_time_Spectra.eps'], 'epsc');
end

%% Vm Spectra zscored across frequencies
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('visible', 'off', 'Renderer', 'Painters', 'Position', [200 200 1000 1000]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
    
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';

        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
    
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            cur_spec_pow(:, :, i) = zscore(abs(cur_spec_pow(:, :, i)), [], 1);
        end

        surface(timeline, nanmean(data_bystim.(f_stim).neuron_spec_freq, 3), nanmean(cur_spec_pow, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        a = colorbar;
        a.Label.String = 'Power (A.U.)';
    
        set(gca, 'color', 'none');
        %avg_power = nanmean(data_bystim.(f_stim).neuron_spec_power, 3);
        xlabel('Time from Stim onset(sec)');
        ylabel('Freq (Hz)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([f_region ' Sub Vm Spectra Z-scored frequency'], 'Interpreter', 'none');

    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_freq_Spectra.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_freq_Spectra.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_zscore_freq_Spectra.eps'], 'epsc');
end


% Period spectrum  with (x - A)/(A + B) normalization for each neuron and 
% Averaged afterwards
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('Renderer', 'Painters', 'Position', [200 200 900 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
      
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        
        % Normalize ratio
        % First calculate by doing the ratios for each neuron first and then averaging
        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
    
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
        
        base_spec = [];
        trans_spec = [];
        sus_spec = [];
	offset_trans_spec = [];
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            % Get the indices define the periods within in the trace
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
	    stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= data_bystim.(f_stim).stim_timestamps(1, i) & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= data_bystim.(f_stim).stim_timestamps(end, i));
            trans_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= trans_ped(1) & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= trans_ped(2)./1000);
            sus_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) > sus_ped(1)./1000 & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= sus_ped(2)./1000);
            offset_trans_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) > offset_trans_ped(1)./1000 & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= offset_trans_ped(2)./1000);
            %plot(baseline_idx, i*5, '.r');
            %hold on;
            base_power = nanmean(abs(cur_spec_pow(:, baseline_idx, i)), 2);
            stim_power = nanmean(abs(cur_spec_pow(:, stim_idx, i)), 2);
            cur_spec_pow(:, :, i) = (abs(cur_spec_pow(:, :, i)) - base_power)./(base_power + stim_power);

            % Store the baseline, stim, and offset spectrums
            base_spec = [base_spec, mean(cur_spec_pow(:, baseline_idx, i) , 2, 'omitnan')];
            trans_spec = [trans_spec, mean(cur_spec_pow(:, trans_idx, i) , 2, 'omitnan')];
            sus_spec = [sus_spec, mean(cur_spec_pow(:, sus_idx, i) , 2, 'omitnan')];
            offset_trans_spec = [offset_trans_spec, mean(cur_spec_pow(:, offset_trans_idx, i) , 2, 'omitnan')];
        end

        freqs = nanmean(data_bystim.(f_stim).neuron_spec_freq, 3);
        
        % Average base frequency spectra
        %avg_base = mean(base_spec, 2, 'omitnan');
        %std_base = std(base_spec, 0, 2, 'omitnan');
        %num_neurons = size(base_spec, 2);
        %sem_base = std_base./sqrt(num_neurons);
        %f = fill([freqs; flip(freqs)], [avg_base + sem_base; flipud(avg_base - sem_base)], [0.5 0.5 0.5]);
        %Multi_func.set_fill_properties(f);
        %hold on;
        %plot(freqs, avg_base, '-', 'DisplayName', 'Base');
        %hold on;
 
        % Average transient frequency spectra
        avg_trans = mean(trans_spec, 2, 'omitnan');
        std_trans = std(trans_spec, 0, 2, 'omitnan');
        num_neurons = size(trans_spec, 2);
        sem_trans = std_trans./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_trans + sem_trans; flipud(avg_trans - sem_trans)], [0.5 0.5 0.5], 'HandleVisibility','off');
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_trans, '-', 'DisplayName', 'Transient');
        hold on;

        % Average sustained frequency spectra
        avg_sus = mean(sus_spec, 2, 'omitnan');
        std_sus = std(sus_spec, 0, 2, 'omitnan');
        num_neurons = size(sus_spec, 2);
        sem_sus = std_sus./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_sus + sem_sus; flipud(avg_sus - sem_sus)], [0.5 0.5 0.5], 'HandleVisibility','off');
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_sus, '-', 'DisplayName', 'Sustained');
        hold on;
        
        % Average offset transient frequency spectra
        avg_offset_trans = mean(offset_trans_spec, 2, 'omitnan');
        std_offset_trans = std(offset_trans_spec, 0, 2, 'omitnan');
        num_neurons = size(offset_trans_spec, 2);
        sem_offset_trans = std_offset_trans./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_offset_trans + sem_offset_trans; flipud(avg_offset_trans - sem_offset_trans)], [0.5 0.5 0.5], 'HandleVisibility','off');
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_offset_trans, '-', 'DisplayName', 'Offset Transient');
        hold on;
        plot(1:150, repmat(0, 150, 1), '-k', 'HandleVisibility','off');
        legend('Location','eastoutside');
        
        % Format and label axes
        Multi_func.set_default_axis(gca);
        x = gca; x = x.YAxis;
        Multi_func.set_spacing_axis(x, 4, 2);
        %avg_power = nanmean(data_bystim.(f_stim).neuron_spec_power, 3);
        xlabel('Freq (Hz)');
        ylabel('Relative Power');
	    %xlim([0 20]);
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([ f_region ' Spectrum with (x - A)/(A + B) normalization'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_trans_sus_Spectrum.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_trans_sus_Spectrum.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_A_B_Normalization_trans_sus_Spectrum.eps'], 'epsc');
end

%% Period spectrum of Raw Sub Vm
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    figure('Renderer', 'Painters', 'Position', [200 200 900 700]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        nexttile;
      
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';
        
        % Normalize ratio
        % First calculate by doing the ratios for each neuron first and then averaging
        cur_spec_pow = data_bystim.(f_stim).neuron_spec_power;
    
        % Plot the starting time point for each neuron
        sz = size(data_bystim.(f_stim).trace_timestamps);
        
        base_spec = [];
        stim_spec = [];
        offset_spec = [];
        % Loop throug each neuron
        for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
            % Get the indices define the periods within in the trace
            baseline_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
            stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= data_bystim.(f_stim).stim_timestamps(1, i) & ...
                            data_bystim.(f_stim).trace_timestamps(:, i) <= data_bystim.(f_stim).stim_timestamps(end, i));
            offset_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) > data_bystim.(f_stim).stim_timestamps(end, i));
            cur_spec_pow = abs(cur_spec_pow);
            % Store the baseline, stim, and offset spectrums
            base_spec = [base_spec, mean(cur_spec_pow(:, baseline_idx, i) , 2, 'omitnan')];
            stim_spec = [stim_spec, mean(cur_spec_pow(:, stim_idx, i) , 2, 'omitnan')];
            offset_spec = [base_spec, mean(cur_spec_pow(:, offset_idx, i) , 2, 'omitnan')];
        end

        freqs = nanmean(data_bystim.(f_stim).neuron_spec_freq, 3);
        
        % Average base frequency spectra
        avg_base = mean(base_spec, 2, 'omitnan');
        std_base = std(base_spec, 0, 2, 'omitnan');
        num_neurons = size(base_spec, 2);
        sem_base = std_base./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_base + sem_base; flipud(avg_base - sem_base)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_base, '-', 'DisplayName', 'Base');
        hold on;
 
        % Average stimulation frequency spectra
        avg_stim = mean(stim_spec, 2, 'omitnan');
        std_stim = std(stim_spec, 0, 2, 'omitnan');
        num_neurons = size(stim_spec, 2);
        sem_stim = std_stim./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_stim + sem_stim; flipud(avg_stim - sem_stim)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_stim, '-', 'DisplayName', 'Stim');
        hold on;
        
        % Average offset frequency spectra
        avg_offset = mean(offset_spec, 2, 'omitnan');
        std_offset = std(offset_spec, 0, 2, 'omitnan');
        num_neurons = size(offset_spec, 2);
        sem_offset = std_offset./sqrt(num_neurons);
        f = fill([freqs; flip(freqs)], [avg_offset + sem_offset; flipud(avg_offset - sem_offset)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(f);
        hold on;
        plot(freqs, avg_offset, '-', 'DisplayName', 'Offset');
        
        legend('Location','eastoutside');
        Multi_func.set_default_axis(gca);
        x = gca; x = x.YAxis;
        Multi_func.set_spacing_axis(x, 4, 2);
        xlabel('Freq (Hz)');
        ylabel('Power (A.U.)');
        title(f_stim(3:end), 'Interpreter', 'none');
    end
    sgtitle([ f_region ' Spectrum with Raw Sub Vm'], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_SubVm_Spectrum.png']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_SubVm_Spectrum.pdf']);
    saveas(gcf, [figure_path 'Spectra/' f_region '_Raw_SubVm_Spectrum.eps'], 'epsc');
end

%% Functin to calculate the power spectra
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

