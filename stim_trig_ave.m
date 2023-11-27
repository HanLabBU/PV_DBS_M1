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
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine how much before and after a stimulation pulse to take the average
trace_sur = 10; % This is ~6ms before and after stim pulse

%%% END Modification

% Check that the server path exists
if ~isfolder(server_root_path)
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
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names); %stim_cond(all_matfiles);
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
        data_bystim.(f_stim).neuron_stim_SubVm = [];
        data_bystim.(f_stim).neuron_stim_RawVm = [];
        data_bystim.(f_stim).neuron_trace_sur_time = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_base_inter = [];
            cur_fov_stim_SubVm = [];
            cur_fov_stim_RawVm = [];
            cur_fov_sur_traceTime = [];

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

                    % Grabbing raw traces
                    cur_raw_trace = raw_trial_data.raw_traces(front_frame_drop:back_frame_drop, roi_idx);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_raw_trace, round(trial_data.camera_framerate));
                    detrend_Vm = cur_raw_trace - baseline';
                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{5}));
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop);
                    
                    % Grab subthreshold Vm
                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;

                    % Get the trace idx when there are stimulation pulses
                    % Cannot think of a clever way to vectorize this unfortunately
                    for pulse=cur_stim_time'
                        stim_center = find(min(abs(cur_trace_time - pulse)) == abs(cur_trace_time - pulse));
                        cur_fov_stim_RawVm = horzcat_pad(cur_fov_stim_RawVm, detrend_Vm(stim_center - trace_sur: stim_center + trace_sur));
                        cur_fov_stim_SubVm = horzcat_pad(cur_fov_stim_SubVm, detrend_subVm(stim_center - trace_sur: stim_center + trace_sur)');
                        %TODO find a better way of getting the timestamp here
                        cur_fov_sur_traceTime = horzcat_pad(cur_fov_sur_traceTime, cur_trace_time(stim_center - trace_sur: stim_center + trace_sur) - pulse);
                    end
                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(~isnan(cur_fov_stim_SubVm(:))) < 1 || isempty(cur_fov_stim_SubVm)
                continue;
            end

            % Store both Subthreshold Vm and Raw VM 
            temp = data_bystim.(f_stim).neuron_stim_SubVm;
            data_bystim.(f_stim).neuron_stim_SubVm = horzcat_pad(temp, mean(cur_fov_stim_SubVm, 2 ,'omitnan'));
            temp = data_bystim.(f_stim).neuron_stim_RawVm;
            data_bystim.(f_stim).neuron_stim_RawVm = horzcat_pad(temp, mean(cur_fov_stim_RawVm, 2 ,'omitnan'));
            temp = data_bystim.(f_stim).neuron_trace_sur_time;
            data_bystim.(f_stim).neuron_trace_sur_time = horzcat_pad(temp, mean(cur_fov_sur_traceTime, 2 ,'omitnan'));

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

avg_Fs = nanmean(all_Fs);

% Plot the average stimulation pulse
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region).data_bystim;
    stims = fieldnames(data_bystim);
    figure('Renderer', 'Painters', 'Position', [200 200 2000 700]);
    tiledlayout(length(stims), 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each stimulation parameter
    for f_stim=stims'
        f_stim = f_stim{1};
        
        % Subthreshold Vm
        nexttile;
        timeline = mean(data_bystim.(f_stim).neuron_trace_sur_time, 2, 'omitnan')*1000;
        plot(timeline, data_bystim.(f_stim).neuron_stim_SubVm);
    
        set(gca, 'color', 'none');
        xlabel('Time from Stim pulse(ms)');
        ylabel('Vm');
        title([f_stim(3:end) ' Sub Vm'], 'Interpreter', 'none');
        
        % Raw Vm
        nexttile;
        timeline = mean(data_bystim.(f_stim).neuron_trace_sur_time, 2, 'omitnan')*1000;
        plot(timeline, data_bystim.(f_stim).neuron_stim_RawVm);
    
        set(gca, 'color', 'none');
        xlabel('Time from Stim pulse(ms)');
        ylabel('Vm');
        title([f_stim(3:end) ' Raw Vm'], 'Interpreter', 'none');
    end
    sgtitle([ f_region], 'Interpreter', 'none');
    
    saveas(gcf, [figure_path 'Stim_trig_avg/' f_region '_Pop_Stim_Trig.png']);
    saveas(gcf, [figure_path 'Stim_trig_avg/' f_region '_Pop_Stim_Trig.eps'], 'epsc');
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
