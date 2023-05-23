% Save all analysis data to a matfile that can be easily read for analysis
clear all;
close all;

f = filesep;

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

%TODO may need to change the variables to differentiate between the larger and finer resolution of spike rate
srate_win_large = 100;
srate_win_small = 20;

% Time periods for comparison of firing rate and sub Vm
trans_ped = [0, 150];
sus_ped = [150, 1000];
offset_trans_ped = [1000, 1150];

% Data path for all of the intermediate anaylsis data
save_all_data_file = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Interm_Data' f 'pv_data.mat'];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

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
    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);

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

        % Trace stuff
        data_bystim.(f_stim).framerate = [];
        data_bystim.(f_stim).all_trial_SubVm = {};
        data_bystim.(f_stim).neuron_Vm = [];
        data_bystim.(f_stim).all_trial_rawVm = {};
        data_bystim.(f_stim).all_trial_spikeidx = {};
        data_bystim.(f_stim).neuron_srate_large = [];
        data_bystim.(f_stim).neuron_srate_small = [];
         
        % Power spectra stuff
        data_bystim.(f_stim).neuron_spec_power = [];
        data_bystim.(f_stim).neuron_spec_freq = [];

        % Store transient information    
        data_bystim.(f_stim).neuron_trans_Vm = [];      
        data_bystim.(f_stim).neuron_sus_Vm = [];        
        data_bystim.(f_stim).neuron_trans_FR = [];     
        data_bystim.(f_stim).neuron_sus_FR = [];

        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];
        data_bystim.(f_stim).neuron_name = {};

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            matfile = matfile{1};
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            
            cur_fov_Fs = [];
            cur_fov_subVm = [];
            cur_fov_rawtraces = [];
            cur_fov_spikeidx = [];

            cur_fov_srate_large = [];
            cur_fov_raster = [];

            cur_fov_trans_Vm = [];      
            cur_fov_sus_Vm = [];        
            cur_fov_trans_FR = [];     
            cur_fov_sus_FR = [];
            cur_fov_base_Vm = [];
            cur_fov_base_FR = [];

            cur_fov_stim_time = [];
            cur_fov_trace_time = [];

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(matfile, '_');
                try
                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                catch
                    trial_ignr_list = [];
                end
                
                % Remove trials from trial idx list
                trial_idxs = setdiff(trial_idxs, trial_ignr_list);

                % Skip if there is at most 2 trials
                if length(trial_idxs) <= 2
                    continue;
                end
            
                %Add neuron name
                data_bystim.(f_stim).neuron_name{end + 1} = matfile;

                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;
                    
                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
                    cur_stim_time = raw_trial_data.raw_stimulation_time - stim_start;
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, cur_stim_time);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, cur_trace_time);

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

                    % Calculate spikerate
                    cur_raster = trial_data.spike_info375.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_large, 1, 1);
                    % Baseline subtract the mean
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_large = horzcat_pad(cur_fov_srate_large, cur_spikerate');            
                    
                    % Store the raster plot
                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');

                    % Grab transient and sustained period for spikes and Vm
                    base_idx = find(cur_trace_time < cur_stim_time(1));
                    trans_idx = find(cur_trace_time > 0 & cur_trace_time <= trans_ped(2)./1000);
                    sus_idx = find(cur_trace_time > trans_ped(2)./1000 & cur_trace_time <= cur_stim_time(end));
                    
                    cur_fov_base_Vm = horzcat_pad(cur_fov_base_Vm, mean(detrend_subVm(base_idx)', 'omitnan'));
                    cur_fov_trans_Vm = horzcat_pad(cur_fov_trans_Vm, mean(detrend_subVm(trans_idx)', 'omitnan'));
                    cur_fov_sus_Vm = horzcat_pad(cur_fov_sus_Vm, mean(detrend_subVm(sus_idx)', 'omitnan'));

                    cur_fov_base_FR = horzcat_pad(cur_fov_base_FR, sum(cur_raster(base_idx))');
                    cur_fov_trans_FR = horzcat_pad(cur_fov_trans_FR, sum(cur_raster(trans_idx))./(diff(trans_ped)./1000)');
                    cur_fov_sus_FR = horzcat_pad(cur_fov_sus_FR, sum(cur_raster(sus_idx))./(diff(sus_ped)./1000)');

                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Save all the subthreshold Vm
            data_bystim.(f_stim).all_trial_SubVm{end + 1} = cur_fov_subVm;

            % Save each neuron's average sub vm
            temp = data_bystim.(f_stim).neuron_Vm;
            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, mean(cur_fov_subVm, 2, 'omitnan'));
            
            % Store average spike rate for each neuron
            temp = data_bystim.(f_stim).neuron_srate_large;
            data_bystim.(f_stim).neuron_srate_large = horzcat_pad(temp, nanmean(cur_fov_srate_large, 2));

            % Save all trial raw traces
            temp = data_bystim.(f_stim).all_trial_rawVm;
            data_bystim.(f_stim).all_trial_rawVm{end + 1} = cur_fov_rawtraces;

            % Save all spike idxs
            temp = data_bystim.(f_stim).all_trial_spikeidx;
            data_bystim.(f_stim).all_trial_spikeidx{end + 1} = cur_fov_spikeidx;

            % Store the timestamp data
            temp = data_bystim.(f_stim).stim_timestamps;
            data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_fov_stim_time, 2));
            temp = data_bystim.(f_stim).trace_timestamps;
            data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_fov_trace_time, 2));
            
            % Save framerate
            data_bystim.(f_stim).framerate(end + 1) = mean(cur_fov_Fs, 'omitnan');
                
            % Save all of the transient and sustained information
            temp = data_bystim.(f_stim).neuron_trans_Vm;
            data_bystim.(f_stim).neuron_trans_Vm = horzcat_pad(temp, mean(cur_fov_trans_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
            temp = data_bystim.(f_stim).neuron_sus_Vm;
            data_bystim.(f_stim).neuron_sus_Vm = horzcat_pad(temp, mean(cur_fov_sus_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
            temp = data_bystim.(f_stim).neuron_trans_FR;
            data_bystim.(f_stim).neuron_trans_FR = horzcat_pad(temp, mean(cur_fov_trans_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));
            temp = data_bystim.(f_stim).neuron_sus_FR;
            data_bystim.(f_stim).neuron_sus_FR = horzcat_pad(temp, mean(cur_fov_sus_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));

            % Calculate and save frequency data
            [wt, f] = get_power_spec(nanmean(cur_fov_subVm, 2)', nanmean(cur_fov_Fs));
            temp = data_bystim.(f_stim).neuron_spec_power;
            data_bystim.(f_stim).neuron_spec_power = cat(3, temp, wt);
            temp = data_bystim.(f_stim).neuron_spec_freq;   
            data_bystim.(f_stim).neuron_spec_freq = cat(3, temp, f);

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region) = data_bystim;
end

% Check first if the matfile exists
save([save_all_data_file], 'region_data', '-v7.3');

% Calculate cwt for input signal and 
function [wt, f] = get_power_spec(signal, samp_freq)
    freqLimits = [0 150];
    fb = cwtfilterbank(SignalLength=length(signal),...
                       SamplingFrequency=samp_freq,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(signal, FilterBank=fb);
end
