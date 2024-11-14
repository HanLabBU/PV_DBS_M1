clear all;
close all;
clc;

f = filesep;

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/';
% Windows server
%local_root_path = 'Z:\';

% List path where all of the matfiles are stored
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
%data_path = [server_root_path 'EricLowet' f 'DBS' f 'github' f 'DBS_volt_data_github'];
data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'CA1_Data' f];

exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;

% Path for this data matfile
%save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'ca1_data.mat'];

% Append data to the interm data matfile
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];

%Variables to differentiate between the larger and finer resolution of spike rate
srate_win_100 = 100;
srate_win_50 = 50;
srate_win_20 = 20;
srate_win_10 = 10;
srate_win_3 = 3;

% Time periods for comparison of firing rate and sub Vm
base_ped = Multi_func.base_ped;
trans_ped =Multi_func.trans_ped;
sus_ped = Multi_func.sus_ped;
stim_ped = Multi_func.stim_ped;
offset_trans_ped = Multi_func.offset_trans_ped;

% Check that the server path exists
if ~isfolder(local_root_path) || ~isfolder(data_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([data_path '*mat']);
matfiles = {ses.name};

% Loop through different stimulation conditions
regions = {'r_CA1'};
stim_conditions = {'f_140', 'f_40'};

% Select matfiles by brain region
region_data = struct();
all_Fs = [];
%%
for f_region = regions'
    f_region = f_region{1};

    % Select matfiles by stim specific conditions for all regions
    [matfile_stim] = stim_cond(matfiles); 
    % Select matfiles by stim condition for given region

    % Loop through each field of the struct and concatenate everything together
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies

    % Loop through each stimulation condition
    for f_stim = stim_conditions %fieldnames(matfile_stim)' %
        f_stim = f_stim{1};

        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();

        % Trace stuff
        data_bystim.(f_stim).framerate = [];
        data_bystim.(f_stim).all_trial_SubVm = {};
        data_bystim.(f_stim).neuron_SubVm = [];
        data_bystim.(f_stim).neuron_RawVm = [];
        data_bystim.(f_stim).all_trial_rawVm = {};
        data_bystim.(f_stim).all_trial_trace_noise = {};

        
        % Spike stuff
        data_bystim.(f_stim).all_trial_spikeidx = {};
        data_bystim.(f_stim).all_trial_spike_amp = {};
        data_bystim.(f_stim).neuron_srate_100 = [];
        data_bystim.(f_stim).neuron_srate_50 = [];
        data_bystim.(f_stim).neuron_srate_20 = [];
        data_bystim.(f_stim).neuron_srate_10 = [];
        data_bystim.(f_stim).neuron_srate_3 = [];
        data_bystim.(f_stim).neuron_spikecounts_raster = [];
        data_bystim.(f_stim).all_trial_spike_rasters = {};
        data_bystim.(f_stim).neuron_spike_amp = [];

        % Power spectra stuff
        data_bystim.(f_stim).neuron_subvm_spec_power = [];
        data_bystim.(f_stim).neuron_subvm_spec_freq = [];
        data_bystim.(f_stim).neuron_rawvm_spec_power = [];
        data_bystim.(f_stim).neuron_rawvm_spec_freq = [];
        data_bystim.(f_stim).neuron_subvm_hilbfilt = {};
        data_bystim.(f_stim).neuron_rawvm_hilbfilt = {};
        data_bystim.(f_stim).all_trial_subvm_power_spec = {};        
        data_bystim.(f_stim).all_trial_subvm_spec_freq = {};        
        data_bystim.(f_stim).all_trial_rawvm_power_spec = {};        
        data_bystim.(f_stim).all_trial_rawvm_spec_freq = {};        

        % Store transient information    
        data_bystim.(f_stim).neuron_trans_Vm = [];      
        data_bystim.(f_stim).neuron_sus_Vm = [];        
        data_bystim.(f_stim).neuron_stim_Vm = [];        
        data_bystim.(f_stim).neuron_trans_FR = [];     
        data_bystim.(f_stim).neuron_sus_FR = [];
        data_bystim.(f_stim).neuron_stim_FR = [];

        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];

        % Neuron meta data
        data_bystim.(f_stim).neuron_name = {};
        data_bystim.(f_stim).trial_num = [];

        % Loop through each matfile of the current stimulation condition
        for matf= matfile_stim.(f_stim).names
            cur_mat = matf{1}
            % Read in the mat file of the current condition
            data = load([data_path cur_mat]);

            % double check that data has aligned stuff

            if ~isfield(data, 'align')
                continue;
            end

            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            
            cur_fov_Fs = [];
            cur_fov_subVm = [];
            cur_fov_rawtraces = [];
            cur_fov_tracenoises = [];
            cur_fov_spikeidx = [];
            cur_fov_spike_amp = [];

            cur_fov_srate_100 = [];
            cur_fov_srate_50 = [];
            cur_fov_srate_20 = [];
            cur_fov_srate_10 = [];
            cur_fov_srate_3 = [];
            cur_fov_srate_inst = [];
            cur_fov_raster = [];

            cur_fov_trans_Vm = [];      
            cur_fov_sus_Vm = []; 
            cur_fov_stim_Vm = [];
            cur_fov_trans_FR = [];     
            cur_fov_sus_FR = [];
            cur_fov_stim_FR = [];

            cur_fov_base_Vm = [];
            cur_fov_base_FR = [];

            cur_fov_stim_time = [];
            cur_fov_trace_time = [];

            cur_fov_subvm_wt = [];
            cur_fov_subvm_f = [];
            cur_fov_rawvm_wt = [];
            cur_fov_rawvm_f = [];

            cur_fov_subvm_hilbfilt = [];
            cur_fov_rawvm_hilbfilt = [];

            % Loop through each ROI
            for roi_idx=1 % :size(trial_data.detrend_traces, 2)
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(cur_mat, '_');
                
                try
                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{2}, 'rec')]).(ri{3}).(['f_' ri{4}]).(['ROI' num2str(roi_idx)]);
                catch
                    trial_ignr_list = [];
                end
                
                % Remove trials from trial idx list
                trial_idxs = setdiff(trial_idxs, trial_ignr_list);

                % Skip if there is at most 2 trials
                if length(trial_idxs) <= 2
                    continue;
                end
            

                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    % Skip if there are not any spikes
                    if isempty(trial_data.spike_info)
                        continue;
                    end

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_fov_Fs(end + 1) = trial_data.camera_framerate;
                    
                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{4})) - stim_start;
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, cur_stim_time);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, cur_trace_time);

                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, detrend_subVm');

                    % Grab hilbert transform coefficients of subthreshold Vm
                    [filt_sig] = filt_data(detrend_subVm', [1:1:200], mean(cur_fov_Fs));
                    cur_fov_subvm_hilbfilt(:, :, end + 1) = filt_sig;


                    % Grab the spike idxs
                    cur_spike_idx = trial_data.spike_info.spike_idx{1};
                    cur_spike_idx(find(cur_spike_idx < front_frame_drop | cur_spike_idx > back_frame_drop)) = NaN;
                    cur_spike_idx = cur_spike_idx - front_frame_drop + 1;
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

                    % Grab hilbert transform coefficients of all Vm
                    [filt_sig] = filt_data(detrend_Vm', [1:1:200], mean(cur_fov_Fs));
                    cur_fov_rawvm_hilbfilt(:, :, end + 1) = filt_sig;

                    cur_trace_noise = trial_data.spike_info.trace_noise;
                    cur_fov_tracenoises = horzcat_pad(cur_fov_tracenoises, cur_trace_noise);

                    cur_spike_amp = trial_data.spike_info.spike_amplitude{1};
                    cur_fov_spike_amp = horzcat_pad(cur_fov_spike_amp, cur_spike_amp');

                    % Calculate spikerate with a window size of 100
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    % Uses center window moving average with halved window on the ends
                    %
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_100, 1, 1);
                        
                    %Use preceding window points to average and asign to the last point
                    % |<-------->|  averaging
                    % ...........*  point assigned
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_100);
                    
                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_100 = horzcat_pad(cur_fov_srate_100, cur_spikerate');            
                    
                    
                    % Calculate spikerate with a window size of 50
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    % Uses center window moving average with halved window on the ends
                    %
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_50, 1, 1);
                        
                    %Use preceding window points to average and asign to the last point
                    % |<-------->|  averaging
                    % ...........*  point assigned
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_50);
                    
                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_50 = horzcat_pad(cur_fov_srate_50, cur_spikerate');            
        


                    % Calculate spikerate with window size of 20
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_20);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_20, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_20 = horzcat_pad(cur_fov_srate_20, cur_spikerate');            


                    % Calculate spikerate with a window size of 10
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_10);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_10, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_10 = horzcat_pad(cur_fov_srate_10, cur_spikerate');            

                    

                    % Calculate spikerate with a window size of 3
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_3);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_3, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx), 'omitnan');
                    cur_fov_srate_3 = horzcat_pad(cur_fov_srate_3, cur_spikerate');            


                    % Store the raster plot
                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');

                    % Grab transient, sustained and stimulation period for spikes and Vm
                    base_idx = find(cur_trace_time >= base_ped(1)./1000 & cur_trace_time < base_ped(2)./1000);
                    trans_idx = find(cur_trace_time > trans_ped(1)./1000 & cur_trace_time <= trans_ped(2)./1000);
                    sus_idx = find(cur_trace_time > sus_ped(1)./1000 & cur_trace_time <= cur_stim_time(end));
                    stim_idx = find(cur_trace_time > stim_ped(1)./1000 & cur_trace_time <= stim_ped(2)./1000);

                    cur_fov_base_Vm = horzcat_pad(cur_fov_base_Vm, mean(detrend_subVm(base_idx)', 'omitnan'));
                    cur_fov_trans_Vm = horzcat_pad(cur_fov_trans_Vm, mean(detrend_subVm(trans_idx)', 'omitnan'));
                    cur_fov_sus_Vm = horzcat_pad(cur_fov_sus_Vm, mean(detrend_subVm(sus_idx)', 'omitnan'));
                    cur_fov_stim_Vm = horzcat_pad(cur_fov_stim_Vm, mean(detrend_subVm(stim_idx)', 'omitnan'));

                    cur_fov_base_FR = horzcat_pad(cur_fov_base_FR, sum(cur_raster(base_idx))');
                    cur_fov_trans_FR = horzcat_pad(cur_fov_trans_FR, sum(cur_raster(trans_idx))./(diff(trans_ped)./1000)');
                    cur_fov_sus_FR = horzcat_pad(cur_fov_sus_FR, sum(cur_raster(sus_idx))./(diff(sus_ped)./1000)');
                    cur_fov_stim_FR = horzcat_pad(cur_fov_stim_FR, sum(cur_raster(stim_idx))./(diff(stim_ped)./1000)');

                    % Calculate the power spectra for each trial of subthreshold vm
                    [trial_wt, trial_f] = get_power_spec(detrend_subVm', nanmean(cur_fov_Fs));
                    cur_fov_subvm_wt = cat(3, cur_fov_subvm_wt, abs(trial_wt));
                    cur_fov_subvm_f = cat(3, cur_fov_subvm_f, trial_f);

                    % Calculate the power spectra for each trial of raw vm
                    [trial_wt, trial_f] = get_power_spec(detrend_Vm', nanmean(cur_fov_Fs));
                    cur_fov_rawvm_wt = cat(3, cur_fov_rawvm_wt, abs(trial_wt));
                    cur_fov_rawvm_f = cat(3, cur_fov_rawvm_f, trial_f);
                end 
            end
            % End looping through each neuron

            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            %Add neuron name
            data_bystim.(f_stim).neuron_name{end + 1} = cur_mat;
                
            % Add neuron trial number
            data_bystim.(f_stim).trial_num(end + 1) = size(cur_fov_subVm, 2);

            % Save all the subthreshold Vm
            data_bystim.(f_stim).all_trial_SubVm{end + 1} = cur_fov_subVm;

            % Save each neuron's average sub vm
            temp = data_bystim.(f_stim).neuron_SubVm;
            data_bystim.(f_stim).neuron_SubVm = horzcat_pad(temp, mean(cur_fov_subVm, 2, 'omitnan'));
 
            % Save each neuron's average raw vm
            temp = data_bystim.(f_stim).neuron_RawVm;
            data_bystim.(f_stim).neuron_RawVm = horzcat_pad(temp, mean(cur_fov_rawtraces, 2, 'omitnan'));

            % Store average spike rate for each neuron
            temp = data_bystim.(f_stim).neuron_srate_100;
            data_bystim.(f_stim).neuron_srate_100 = horzcat_pad(temp, nanmean(cur_fov_srate_100, 2));

            temp = data_bystim.(f_stim).neuron_srate_50;
            data_bystim.(f_stim).neuron_srate_50 = horzcat_pad(temp, nanmean(cur_fov_srate_50, 2));

            temp = data_bystim.(f_stim).neuron_srate_20;
            data_bystim.(f_stim).neuron_srate_20 = horzcat_pad(temp, nanmean(cur_fov_srate_20, 2));

            temp = data_bystim.(f_stim).neuron_srate_10;
            data_bystim.(f_stim).neuron_srate_10 = horzcat_pad(temp, nanmean(cur_fov_srate_10, 2));

            temp = data_bystim.(f_stim).neuron_srate_3;
            data_bystim.(f_stim).neuron_srate_3 = horzcat_pad(temp, nanmean(cur_fov_srate_3, 2));

            % Store the spike counts across trials in a raster
            temp = data_bystim.(f_stim).neuron_spikecounts_raster;
            data_bystim.(f_stim).neuron_spikecounts_raster = horzcat_pad(temp, sum(cur_fov_raster, 2));

            % Store the spike rasters for all trials
            temp = data_bystim.(f_stim).all_trial_spike_rasters;
            data_bystim.(f_stim).all_trial_spike_rasters{end + 1} = cur_fov_raster;

            % Save all trial raw traces
            temp = data_bystim.(f_stim).all_trial_rawVm;
            data_bystim.(f_stim).all_trial_rawVm{end + 1} = cur_fov_rawtraces;
            
            % Save all the trial trace noises
            temp = data_bystim.(f_stim).all_trial_trace_noise;
            data_bystim.(f_stim).all_trial_trace_noise{end + 1} = cur_fov_tracenoises;

            % Save all the trial spike amplitudes
            temp = data_bystim.(f_stim).all_trial_spike_amp;
            data_bystim.(f_stim).all_trial_spike_amp{end + 1} = cur_fov_spike_amp;

            % Average the spike amplitudes for each neuron
            temp = data_bystim.(f_stim).neuron_spike_amp;
            data_bystim.(f_stim).neuron_spike_amp = horzcat_pad(temp, nanmean(cur_fov_spike_amp, 'all'));

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
            %temp = data_bystim.(f_stim).neuron_trans_Vm;
            %data_bystim.(f_stim).neuron_trans_Vm = horzcat_pad(temp, mean(cur_fov_trans_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
            %temp = data_bystim.(f_stim).neuron_sus_Vm;
            %data_bystim.(f_stim).neuron_sus_Vm = horzcat_pad(temp, mean(cur_fov_sus_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
            %temp = data_bystim.(f_stim).neuron_stim_Vm;
            %data_bystim.(f_stim).neuron_stim_Vm = horzcat_pad(temp, mean(cur_fov_stim_Vm, 'omitnan') - mean(cur_fov_base_Vm, 'omitnan'));
            %
            %temp = data_bystim.(f_stim).neuron_trans_FR;
            %data_bystim.(f_stim).neuron_trans_FR = horzcat_pad(temp, mean(cur_fov_trans_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));
            %temp = data_bystim.(f_stim).neuron_sus_FR;
            %data_bystim.(f_stim).neuron_sus_FR = horzcat_pad(temp, mean(cur_fov_sus_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));
            %temp = data_bystim.(f_stim).neuron_stim_FR;
            %data_bystim.(f_stim).neuron_stim_FR = horzcat_pad(temp, mean(cur_fov_stim_FR, 'omitnan') - mean(cur_fov_base_FR, 'omitnan'));

            % Calculate and save frequency data
            % Old way of just using the neuron average Vm to start the Spectra calculation
            %[wt, f] = get_power_spec(nanmean(cur_fov_subVm, 2)', nanmean(cur_fov_Fs));

            % Calculate subthreshold vm
            temp = data_bystim.(f_stim).neuron_subvm_spec_power;
            data_bystim.(f_stim).neuron_subvm_spec_power = cat(3, temp, mean(cur_fov_subvm_wt, 3, 'omitnan'));
            temp = data_bystim.(f_stim).neuron_subvm_spec_freq;   
            data_bystim.(f_stim).neuron_subvm_spec_freq = cat(3, temp, mean(cur_fov_subvm_f, 3, 'omitnan'));
            
            % Calculate all Vm
            temp = data_bystim.(f_stim).neuron_rawvm_spec_power;
            data_bystim.(f_stim).neuron_rawvm_spec_power = cat(3, temp, mean(cur_fov_rawvm_wt, 3, 'omitnan'));
            temp = data_bystim.(f_stim).neuron_rawvm_spec_freq;   
            data_bystim.(f_stim).neuron_rawvm_spec_freq = cat(3, temp, mean(cur_fov_rawvm_f, 3, 'omitnan'));
            
            % Save individual trial subthreshold Vm
            data_bystim.(f_stim).all_trial_subvm_power_spec{end + 1} = cur_fov_subvm_wt;
            data_bystim.(f_stim).all_trial_subvm_spec_freq{end + 1} = cur_fov_subvm_f;

            % Save individual trial all Vm
            data_bystim.(f_stim).all_trial_rawvm_power_spec{end + 1} = cur_fov_rawvm_wt;
            data_bystim.(f_stim).all_trial_rawvm_spec_freq{end + 1} = cur_fov_rawvm_f;

            % Save the hilbert frequency data
            cur_fov_subvm_hilbfilt(:, :, 1) = []; % Remove empty first element
            data_bystim.(f_stim).neuron_subvm_hilbfilt{end + 1} = cur_fov_subvm_hilbfilt;

            cur_fov_rawvm_hilbfilt(:, :, 1) = []; % Remove empty first element
            data_bystim.(f_stim).neuron_rawvm_hilbfilt{end + 1} = cur_fov_rawvm_hilbfilt;

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region) = data_bystim;
end

%% Save CA1 data to matfile

% Check first if the matfile exists
if isfile(save_all_data_file)
    data = load(save_all_data_file);
    
    %Save all of the region data from the data file to this structure
    region_data.r_M1 = data.region_data.r_M1;
    region_data.r_V1 = data.region_data.r_V1;
    save([save_all_data_file], 'region_data', '-v7.3');
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'ca1_data.mat'];
end
%save([save_all_data_file], 'region_data', '-v7.3');

%%

% Return matfiles by stimulation condition
function [cond_struct] = stim_cond(matfile_names)
    cond_struct = struct();
    
    % Loop through each matfilename and group by stimulation conditions
    for i=1:length(matfile_names)
        file_parts = split(matfile_names{i}, '_');
            
        if length(file_parts) < 4
            continue;
        end

        stim = file_parts{4};
        
        % Create stimulation field if it does not exist
        if ~isfield(cond_struct, ['f_' stim])
            cond_struct.(['f_' stim]).names = {};
        end

        cond_struct.(['f_' stim]).names{end+1} = matfile_names{i};
    end
end

% Calculate cwt for input signal and 
function [wt, f] = get_power_spec(signal, samp_freq)
    freqLimits = [0 200];
    fb = cwtfilterbank(SignalLength=length(signal),...
                       SamplingFrequency=samp_freq,...
                       FrequencyLimits=freqLimits);
    [wt, f] = cwt(signal, FilterBank=fb);
end

% Filter data with hilbert transform
function  [filt_sig]=filt_data(sig,frs, FS)
    Fn = FS/2;
    for steps=frs;
        FB=[ frs(steps)*0.8 frs(steps)*1.2];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        filt_sig(steps, :)= hilbert(filtfilt(B,A,sig));
    end
end
