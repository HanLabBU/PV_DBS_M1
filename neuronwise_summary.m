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

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

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
        data_bystim.(f_stim).neuron_SubVm = {};
        data_bystim.(f_stim).neuron_spikeidx = {};
        data_bystim.(f_stim).neuron_rawVm = {};
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
            cur_fov_spikeidx = [];
            cur_fov_rawtraces = [];
            cur_fov_stim_time = [];
            cur_fov_trace_time = [];

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(matfile{1}, '_');
                try
                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                catch
                    trial_ignr_list = [];
                end
                
                % Remove trials from trial idx list
                trial_idxs = setdiff(trial_idxs, trial_ignr_list);

                % Skip if there is only 1 trial
                if length(trial_idxs) < 2
                    continue;
                end

                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};


                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;

                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');

                    % Grab the spike idxs
                    cur_spike_idx = trial_data.spike_info375.spike_idx{1};
                    cur_spike_idx(find(cur_spike_idx < front_frame_drop | cur_spike_idx > back_frame_drop)) = NaN;
                    cur_spike_idx = cur_spike_idx - front_frame_drop;
                    
                    % Add if spikes were detected
                    if length(cur_spike_idx)  == 0
                        cur_spike_idx = [NaN];
                    end

                    cur_fov_spikeidx = horzcat_pad(cur_fov_spikeidx, cur_spike_idx);
                    
                    % Grab the raw traces
                    cur_raw_trace = raw_trial_data.raw_traces(front_frame_drop:back_frame_drop, roi_idx);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_raw_trace, round(trial_data.camera_framerate));
                    detrend_Vm = cur_raw_trace - baseline';
                    cur_fov_rawtraces = horzcat_pad(cur_fov_rawtraces, detrend_Vm);

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_SubVm;
            data_bystim.(f_stim).neuron_SubVm{end + 1} = cur_fov_subVm;

            % Store the timestamp data
            temp = data_bystim.(f_stim).stim_timestamps;
            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
            temp = data_bystim.(f_stim).trace_timestamps;
            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
            
            % Save the spike idxs
            temp = data_bystim.(f_stim).neuron_spikeidx;
            data_bystim.(f_stim).neuron_spikeidx{end + 1} = cur_fov_spikeidx;

            % Save the raw traces
            temp = data_bystim.(f_stim).neuron_rawVm;
            data_bystim.(f_stim).neuron_rawVm{end + 1} = cur_fov_rawtraces;

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

% Check if combining all of the regions or not
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

% Calculate the sampling frequency from all of the 
avg_Fs = nanmean(all_Fs);

% Plot each Region and each frequency raster, raw Vm, and subthreshold Vm
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
      
        % Get the trace timestamps
        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2)';

        % Create a figure that includes: raw, spikes, and SubVm
        figure('Renderer', 'Painters', 'Position', [200 200 2000 1000]);
        tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        % Plot the raw traces
        nexttile;  
        data_map = [];
        for i = 1:length(data_bystim.(f_stim).neuron_SubVm)
            data_map = [data_map; data_bystim.(f_stim).neuron_rawVm{i}'];
            data_map(end + 1, :) = NaN(1, size(data_bystim.(f_stim).neuron_rawVm{i}, 1));
        end
        imagesc('XData', timeline, 'YData', 1:size(data_map, 1), 'CData', data_map);

        a = colorbar;
        a.Label.String = 'Vm';
    
        set(gca, 'color', 'none');
        xlabel('Time from Stim onset(sec)');
        xlim([-1 2.25]);
        ylim([0 size(data_map, 1)]);
        ylabel('Neuron Trials');
        title('Raw Traces', 'Interpreter', 'none');
        
        % Plot the Subthreshold Vm
        nexttile;  
        data_map = [];
        for i = 1:length(data_bystim.(f_stim).neuron_SubVm)
            data_map = [data_map; data_bystim.(f_stim).neuron_SubVm{i}'];
            data_map(end + 1, :) = NaN(1, size(data_bystim.(f_stim).neuron_SubVm{i}, 1));
        end
        imagesc('XData', timeline, 'YData', 1:size(data_map, 1), 'CData', data_map);

        a = colorbar;
        a.Label.String = 'Vm';
    
        set(gca, 'color', 'none');
        xlabel('Time from Stim onset(sec)');
        xlim([-1 2.25]);
        ylim([0 size(data_map, 1)]);
        ylabel('Neuron Trials');
        title('Subthreshold Vm', 'Interpreter', 'none');

        % Plot the raster plot
        nexttile;  
        index = 1;
        for fov = 1:length(data_bystim.(f_stim).neuron_spikeidx)
            cur_color = [rand, rand, rand]*0.7;
            cur_fov = data_bystim.(f_stim).neuron_spikeidx{fov};
            for tr = 1:size(cur_fov, 2)
                cur_spikeidx = cur_fov(:, tr);
                cur_spikeidx(isnan(cur_spikeidx)) = [];
                plot(timeline(cur_spikeidx), repmat(index, length(cur_spikeidx), 1), '|', 'color', cur_color);
                hold on;
                index = index + 1;
            end
            plot(timeline, repmat(index, length(timeline), 1), '-k');
            hold on;
            index = index + 1;
        end
        
        % DEBUG
        disp([f_region ' ' f_stim ' ' num2str(index)]);
        
        %a = colorbar;
        %a.Label.String = 'Vm';
    
        set(gca, 'color', 'none');
        xlabel('Time from Stim onset(sec)');
        %xlim([-1 2.25]);
        ylim([0 index]);
        ylabel('Neuron Trials');
        title('Raster Plots', 'Interpreter', 'none');

        sgtitle([ f_region ' ' f_stim ' Neuronwise'], 'Interpreter', 'none');

        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.png']);
        saveas(gcf, [figure_path 'Neuronwise/' f_region '_' f_stim '_Neuronwise_Plot.eps'], 'epsc');
    end
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

