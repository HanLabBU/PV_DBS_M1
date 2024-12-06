%% Save all analysis data to a matfile that can be easily read for analysis
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

exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;


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

% Data path for all of the intermediate anaylsis data
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Check that the server path exists
if ~isfolder(local_root_path) || ~isfolder(pv_data_path)
    disp('Server rootpath does not exist!!!');
    return;
end

%% Structure that contains the surgery date based on mouse ID
% V1 mice surgery dates
surg_date_dict.m_109557 = '20231113';
surg_date_dict.m_109558 = '20231113';
surg_date_dict.m_23072 = '20211102';
surg_date_dict.m_611284 = '20210730';
surg_date_dict.m_96334 = '20231026';
surg_date_dict.m_109567 = '20231026';

% M1 mice surgery dates
surg_date_dict.m_31556eartag = '20220817';
surg_date_dict.m_31556noeartag = '20220819';
surg_date_dict.m_50373 = '20230621';
surg_date_dict.m_50464 = '20230308';
surg_date_dict.m_60430 = '20220901'; % Should be excluded
surg_date_dict.m_81631 = '20230508'; % Should be excluded
surg_date_dict.m_617100 = ''; % Do not have cage card unfortunately


%% Setup
ses = dir([pv_data_path '*.mat']);
all_matfiles = {ses.name};

% Select matfiles by brain region
[region_matfiles] = Multi_func.find_region(all_matfiles);
region_data = struct();
all_Fs = [];
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    % Select matfiles by stim specific conditions for all regions
    %[matfile_stim] = stim_cond(all_matfiles); 
    % Select matfiles by stim condition for given region
    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);

    % Loop through each field of the struct and concatenate everything together
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
        data_bystim.(f_stim).surgery_date = struct;
        data_bystim.(f_stim).trial_num = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            matfile = matfile{1};
            
            tokens = regexp(matfile, '_', 'split');

            % Determine the mouse ID and set the surgery date
            % Check if the mouse ID is already listed
            mouse_id_f = ['m_' tokens{1}];
            if ~isfield(data_bystim.(f_stim).surgery_date, mouse_id_f)
                data_bystim.(f_stim).surgery_date.(mouse_id_f) = surg_date_dict.(mouse_id_f);
            end

            % Read in the mat file of the current condition
            data = load([pv_data_path matfile]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            
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
                
                % Collect all data for individual ROIs
                cur_roi_Fs = [];
                cur_roi_subVm = [];
                cur_roi_rawtraces = [];
                cur_roi_tracenoises = [];
                cur_roi_spikeidx = [];
                cur_roi_spike_amp = [];

                cur_roi_srate_100 = [];
                cur_roi_srate_50 = [];
                cur_roi_srate_20 = [];
                cur_roi_srate_10 = [];
                cur_roi_srate_3 = [];
                cur_roi_srate_inst = [];
                cur_roi_raster = [];

                cur_roi_trans_Vm = [];      
                cur_roi_sus_Vm = []; 
                cur_roi_stim_Vm = [];
                cur_roi_trans_FR = [];     
                cur_roi_sus_FR = [];
                cur_roi_stim_FR = [];

                cur_roi_base_Vm = [];
                cur_roi_base_FR = [];

                cur_roi_stim_time = [];
                cur_roi_trace_time = [];

                cur_roi_subvm_wt = [];
                cur_roi_subvm_f = [];
            
                cur_roi_rawvm_wt = [];
                cur_roi_rawvm_f = [];

                cur_roi_subvm_hilbfilt = [];
                cur_roi_rawvm_hilbfilt = [];

                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    % Store the camera framerate
                    all_Fs(end+1) = trial_data.camera_framerate;
                    cur_roi_Fs(end + 1) = trial_data.camera_framerate;

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    
                    % Check if there is an interpolated trace (indicates a TICO recording)
                    if isfield(trial_data, 'interp_time')
                        cur_trace_time = trial_data.interp_time(front_frame_drop:back_frame_drop) - stim_start;
                    else
                        cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
                    end

                    % Perform sanity check for number of timestamps collected
                    if abs(str2num(ri{5}) - length(raw_trial_data.raw_stimulation_time)) > 10
                        disp(['Length of stim timestamps' length(raw_trial_data.raw_stimulation_time)]);
                        disp(['Labeled stim freq ' str2num(ri{5})]);
                        pause();
                    end

                    cur_stim_time = raw_trial_data.raw_stimulation_time(1:str2num(ri{5})) - stim_start;
                    cur_roi_stim_time = horzcat_pad(cur_roi_stim_time, cur_stim_time);
                    cur_roi_trace_time = horzcat_pad(cur_roi_trace_time, cur_trace_time(:));

                    % Grab the subthreshold Vm
                    % 'trace_ws' uses the orignal trace and replaces the spikes with the neighboring interpolation
                    cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, :);
                    
                    if isfield(trial_data, 'interp_time')
                        % Interpolating the subthreshold trace 
                        cur_trace_ws = interp1(trial_data.camera_frame_time, cur_trace_ws, trial_data.interp_time);
                    end
                    cur_trace_ws = cur_trace_ws(front_frame_drop:back_frame_drop);

                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));
                    detrend_subVm = cur_trace_ws - baseline;
                    cur_roi_subVm = horzcat_pad(cur_roi_subVm, detrend_subVm');

                    % Grab hilbert transform coefficients of subthreshold Vm
                    [filt_sig] = Multi_func.filt_data(detrend_subVm', [1:1:200], mean(cur_roi_Fs));
                    cur_roi_subvm_hilbfilt(:, :, end + 1) = filt_sig;


                    % Grab the spike idxs
                    cur_spike_idx = trial_data.spike_info375.spike_idx{roi_idx};
                    cur_raster = trial_data.spike_info375.roaster(roi_idx, :);
                    if isfield(trial_data, 'interp_time')
                        % adjust when using the spike raster
                        raster_interp = interp1(trial_data.camera_frame_time, cur_raster, trial_data.interp_time);
                        [~, cur_spike_idx] = findpeaks(raster_interp);
                        cur_raster = zeros(size(trial_data.interp_time));
                        cur_raster(cur_spike_idx) = 1;

                    end
                    cur_raster = cur_raster(front_frame_drop:back_frame_drop);
            
                    cur_spike_idx(find(cur_spike_idx < front_frame_drop | cur_spike_idx > back_frame_drop)) = NaN;
                    cur_spike_idx = cur_spike_idx - front_frame_drop + 1;
                    % Add if spikes were detected
                    if length(cur_spike_idx)  == 0
                        cur_spike_idx = [NaN];
                    end
                    cur_roi_spikeidx = horzcat_pad(cur_roi_spikeidx, cur_spike_idx(:));
                    
                    % Grab the raw traces
                    % Check if there is an interpolated trace (indicates a TICO recording)
                    if isfield(raw_trial_data, 'interp_raw_traces')
                        cur_raw_trace = raw_trial_data.interp_raw_traces(front_frame_drop:back_frame_drop, roi_idx);
                    else
                        cur_raw_trace = raw_trial_data.raw_traces(front_frame_drop:back_frame_drop, roi_idx);
                    end

                    [baseline, coeff] = Multi_func.exp_fit_Fx(cur_raw_trace, round(trial_data.camera_framerate));
                    detrend_Vm = cur_raw_trace - baseline';
                    cur_roi_rawtraces = horzcat_pad(cur_roi_rawtraces, detrend_Vm);
                    
                    % Grab hilbert transform coefficients of all Vm
                    [filt_sig] = Multi_func.filt_data(detrend_Vm', [1:1:200], mean(cur_roi_Fs));
                    cur_roi_rawvm_hilbfilt(:, :, end + 1) = filt_sig;

                    %if isfield(trial_data, 'interp_time')
                    %    figure;
                    %    plot(cur_raster);
                    %    figure;
                    %    plot(detrend_Vm, '-b');
                    %    hold on;
                    %    plot(cur_raster.*detrend_Vm', '.r');
                    %    pause();
                    %end
                    
                    cur_trace_noise = trial_data.spike_info375.trace_noise(roi_idx);
                    cur_roi_tracenoises = horzcat_pad(cur_roi_tracenoises, cur_trace_noise);

                    cur_spike_amp = trial_data.spike_info375.spike_amplitude{roi_idx};
                    
                    % set to Nan if there is no spikes for this trial
                    if isempty(cur_spike_amp)
                        cur_spike_amp = [NaN];
                    end

                    cur_roi_spike_amp = horzcat_pad(cur_roi_spike_amp, cur_spike_amp');

                    % Calculate spikerate with a window size of 100
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
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx));
                    cur_roi_srate_100 = horzcat_pad(cur_roi_srate_100, cur_spikerate');            
                    
                    
                    % Calculate spikerate with a window size of 50
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
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx));
                    cur_roi_srate_50 = horzcat_pad(cur_roi_srate_50, cur_spikerate');            
        
                    % Calculate spikerate with window size of 20
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_20);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_20, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx));
                    cur_roi_srate_20 = horzcat_pad(cur_roi_srate_20, cur_spikerate');            

                    % Calculate spikerate with a window size of 10
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_10);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_10, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx));
                    cur_roi_srate_10 = horzcat_pad(cur_roi_srate_10, cur_spikerate');            

                    % Calculate spikerate with a window size of 3
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;

                    %Use preceding window points to average and asign to the last point
                    %cur_spikerate = Multi_func.estimate_spikerate(cur_raster, trial_data.camera_framerate, srate_win_3);

                    % Uses center window moving average   
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win_3, 1, 1);

                    % Baseline subtract the mean 
                    baseline_idx = find(cur_trace_time < cur_stim_time(1));
                    cur_spikerate = cur_spikerate - mean(cur_spikerate(baseline_idx));
                    cur_roi_srate_3 = horzcat_pad(cur_roi_srate_3, cur_spikerate');            


                    % Store the raster plot
                    cur_roi_raster = horzcat_pad(cur_roi_raster, cur_raster');

                    % Grab transient, sustained and stimulation period for spikes and Vm
                    base_idx = find(cur_trace_time >= base_ped(1)./1000 & cur_trace_time < base_ped(2)./1000);
                    trans_idx = find(cur_trace_time > trans_ped(1)./1000 & cur_trace_time <= trans_ped(2)./1000);
                    sus_idx = find(cur_trace_time > sus_ped(1)./1000 & cur_trace_time <= cur_stim_time(end));
                    stim_idx = find(cur_trace_time > stim_ped(1)./1000 & cur_trace_time <= stim_ped(2)./1000);

                    cur_roi_base_Vm = horzcat_pad(cur_roi_base_Vm, mean(detrend_subVm(base_idx)'));
                    cur_roi_trans_Vm = horzcat_pad(cur_roi_trans_Vm, mean(detrend_subVm(trans_idx)'));
                    cur_roi_sus_Vm = horzcat_pad(cur_roi_sus_Vm, mean(detrend_subVm(sus_idx)'));
                    cur_roi_stim_Vm = horzcat_pad(cur_roi_stim_Vm, mean(detrend_subVm(stim_idx)'));

                    cur_roi_base_FR = horzcat_pad(cur_roi_base_FR, sum(cur_raster(base_idx))');
                    cur_roi_trans_FR = horzcat_pad(cur_roi_trans_FR, sum(cur_raster(trans_idx))./(diff(trans_ped)./1000)');
                    cur_roi_sus_FR = horzcat_pad(cur_roi_sus_FR, sum(cur_raster(sus_idx))./(diff(sus_ped)./1000)');
                    cur_roi_stim_FR = horzcat_pad(cur_roi_stim_FR, sum(cur_raster(stim_idx))./(diff(stim_ped)./1000)');

                    % Calculate the power spectra for each trial of subthreshold vm
                    [trial_wt, trial_f] = get_power_spec(detrend_subVm', nanmean(cur_roi_Fs));
                    cur_roi_subvm_wt = cat(3, cur_roi_subvm_wt, abs(trial_wt));
                    cur_roi_subvm_f = cat(3, cur_roi_subvm_f, trial_f);

                    % Calculate the power spectra for each trial of raw vm
                    [trial_wt, trial_f] = get_power_spec(detrend_Vm', nanmean(cur_roi_Fs));
                    cur_roi_rawvm_wt = cat(3, cur_roi_rawvm_wt, abs(trial_wt));
                    cur_roi_rawvm_f = cat(3, cur_roi_rawvm_f, trial_f);

                end % End looping through trials

                % Skip rest of the calculations if the subthreshold Vm is nan
                if sum(isnan(cur_roi_subVm(:))) > 0 || isempty(cur_roi_subVm)
                    continue;
                end

                %Add neuron name
                data_bystim.(f_stim).neuron_name{end + 1} = [matfile '_' num2str(roi_idx)];
                    
                % Add neuron trial number
                data_bystim.(f_stim).trial_num(end + 1) = size(cur_roi_subVm, 2);

                % Save all the subthreshold Vm
                data_bystim.(f_stim).all_trial_SubVm{end + 1} = cur_roi_subVm;

                % Save each neuron's average sub vm
                temp = data_bystim.(f_stim).neuron_SubVm;
                data_bystim.(f_stim).neuron_SubVm = horzcat_pad(temp, mean(cur_roi_subVm, 2));
                
                % Save all trial raw traces
                temp = data_bystim.(f_stim).all_trial_rawVm;
                data_bystim.(f_stim).all_trial_rawVm{end + 1} = cur_roi_rawtraces;

                % Save each neuron's average raw vm
                temp = data_bystim.(f_stim).neuron_RawVm;
                data_bystim.(f_stim).neuron_RawVm = horzcat_pad(temp, mean(cur_roi_rawtraces, 2));

                % Store average spike rate for each neuron
                temp = data_bystim.(f_stim).neuron_srate_100;
                data_bystim.(f_stim).neuron_srate_100 = horzcat_pad(temp, nanmean(cur_roi_srate_100, 2));

                temp = data_bystim.(f_stim).neuron_srate_50;
                data_bystim.(f_stim).neuron_srate_50 = horzcat_pad(temp, nanmean(cur_roi_srate_50, 2));

                temp = data_bystim.(f_stim).neuron_srate_20;
                data_bystim.(f_stim).neuron_srate_20 = horzcat_pad(temp, nanmean(cur_roi_srate_20, 2));

                temp = data_bystim.(f_stim).neuron_srate_10;
                data_bystim.(f_stim).neuron_srate_10 = horzcat_pad(temp, nanmean(cur_roi_srate_10, 2));

                temp = data_bystim.(f_stim).neuron_srate_3;
                data_bystim.(f_stim).neuron_srate_3 = horzcat_pad(temp, nanmean(cur_roi_srate_3, 2));

                % Store the spike counts across trials in a raster
                temp = data_bystim.(f_stim).neuron_spikecounts_raster;
                data_bystim.(f_stim).neuron_spikecounts_raster = horzcat_pad(temp, sum(cur_roi_raster, 2));

                % Store the spike rasters for all trials
                temp = data_bystim.(f_stim).all_trial_spike_rasters;
                data_bystim.(f_stim).all_trial_spike_rasters{end + 1} = cur_roi_raster;

                % Save all the trial trace noises
                temp = data_bystim.(f_stim).all_trial_trace_noise;
                data_bystim.(f_stim).all_trial_trace_noise{end + 1} = cur_roi_tracenoises;

                % Save all the trial spike amplitudes
                temp = data_bystim.(f_stim).all_trial_spike_amp;
                data_bystim.(f_stim).all_trial_spike_amp{end + 1} = cur_roi_spike_amp;

                % Average the spike amplitudes for each neuron
                temp = data_bystim.(f_stim).neuron_spike_amp;
                data_bystim.(f_stim).neuron_spike_amp = horzcat_pad(temp, nanmean(cur_roi_spike_amp, 'all'));

                % Save all spike idxs
                temp = data_bystim.(f_stim).all_trial_spikeidx;
                data_bystim.(f_stim).all_trial_spikeidx{end + 1} = cur_roi_spikeidx;

                % Store the timestamp data
                temp = data_bystim.(f_stim).stim_timestamps;
                data_bystim.(f_stim).stim_timestamps = horzcat_pad(temp, nanmean(cur_roi_stim_time, 2));
                temp = data_bystim.(f_stim).trace_timestamps;
                data_bystim.(f_stim).trace_timestamps = horzcat_pad(temp, nanmean(cur_roi_trace_time, 2));
                
                % Save framerate
                data_bystim.(f_stim).framerate(end + 1) = mean(cur_roi_Fs);
                    
                % Save all of the transient and sustained information
                %temp = data_bystim.(f_stim).neuron_trans_Vm;
                %data_bystim.(f_stim).neuron_trans_Vm = horzcat_pad(temp, mean(cur_roi_trans_Vm) - mean(cur_roi_base_Vm));
                %temp = data_bystim.(f_stim).neuron_sus_Vm;
                %data_bystim.(f_stim).neuron_sus_Vm = horzcat_pad(temp, mean(cur_roi_sus_Vm) - mean(cur_roi_base_Vm));
                %temp = data_bystim.(f_stim).neuron_stim_Vm;
                %data_bystim.(f_stim).neuron_stim_Vm = horzcat_pad(temp, mean(cur_roi_stim_Vm) - mean(cur_roi_base_Vm));
                
                %temp = data_bystim.(f_stim).neuron_trans_FR;
                %data_bystim.(f_stim).neuron_trans_FR = horzcat_pad(temp, mean(cur_roi_trans_FR) - mean(cur_roi_base_FR));
                %temp = data_bystim.(f_stim).neuron_sus_FR;
                %data_bystim.(f_stim).neuron_sus_FR = horzcat_pad(temp, mean(cur_roi_sus_FR) - mean(cur_roi_base_FR));
                %temp = data_bystim.(f_stim).neuron_stim_FR;
                %data_bystim.(f_stim).neuron_stim_FR = horzcat_pad(temp, mean(cur_roi_stim_FR) - mean(cur_roi_base_FR));

                % Calculate and save frequency data
                % Old way of just using the neuron average Vm to start the Spectra calculation
                %[wt, f] = get_power_spec(nanmean(cur_roi_subVm, 2)', nanmean(cur_roi_Fs));

                % Save power spectrum of subthreshold vm
                temp = data_bystim.(f_stim).neuron_subvm_spec_power;
                data_bystim.(f_stim).neuron_subvm_spec_power = cat(3, temp, mean(cur_roi_subvm_wt, 3));
                temp = data_bystim.(f_stim).neuron_subvm_spec_freq;   
                data_bystim.(f_stim).neuron_subvm_spec_freq = cat(3, temp, mean(cur_roi_subvm_f, 3));

                % Save the power spectrum of all Vm
                temp = data_bystim.(f_stim).neuron_rawvm_spec_power;
                data_bystim.(f_stim).neuron_rawvm_spec_power = cat(3, temp, mean(cur_roi_rawvm_wt, 3));
                temp = data_bystim.(f_stim).neuron_rawvm_spec_freq;   
                data_bystim.(f_stim).neuron_rawvm_spec_freq = cat(3, temp, mean(cur_roi_rawvm_f, 3));

                % Individual trial subthreshold power spectra
                data_bystim.(f_stim).all_trial_subvm_power_spec{end + 1} = cur_roi_subvm_wt;
                data_bystim.(f_stim).all_trial_subvm_spec_freq{end + 1} = cur_roi_subvm_f;

                % Individual trial all Vm power spectra
                data_bystim.(f_stim).all_trial_rawvm_power_spec{end + 1} = cur_roi_rawvm_wt;
                data_bystim.(f_stim).all_trial_rawvm_spec_freq{end + 1} = cur_roi_rawvm_f;

                % Save the hilbert frequency data
                cur_roi_subvm_hilbfilt(:, :, 1) = []; % Remove empty first element
                data_bystim.(f_stim).neuron_subvm_hilbfilt{end + 1} = cur_roi_subvm_hilbfilt;

                cur_roi_rawvm_hilbfilt(:, :, 1) = []; % Remove empty first element
                data_bystim.(f_stim).neuron_rawvm_hilbfilt{end + 1} = cur_roi_rawvm_hilbfilt;

            end % End looping through each ROI
        end % End looping through FOVs 
    end % End looping through

    % Save the VM to the specific region
    region_data.(f_region) = data_bystim;
end

% Check if there is already an interm_data file
% If there is, then just replace the regions that were created
if isfile(save_all_data_file)
    data = load(save_all_data_file);
    
    % Loop through all of the regions that already existed
    for field = fieldnames(data.region_data)'
        field = field{1};

        % Only save field if it was not constructed in this consolidation script
        if ~isfield(region_data, field)
            disp(['Preserved region ' field]);
            region_data.(field) = data.region_data.(field);
        end
    end
end

save([save_all_data_file], 'region_data', '-v7.3');

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
