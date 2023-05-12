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
all_regions = 1;

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
for f_region = {'r_M1'} %fieldnames(region_matfiles)'
    f_region = f_region{1};

    %% Select matfiles by stim specific conditions for all regions
    %[matfile_stim] = stim_cond(all_matfiles); 
    %% Select matfiles by stim condition for given region
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);

    %% Loop through each field of the struct and concatenate everything together
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies

    % Loop through 40Hz stimulated neurons
    for f_stim = {'f_40'}
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    
    
        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();
        data_bystim.(f_stim).neuron_Vm = [];
        data_bystim.(f_stim).neuron_spec_power = [];
        data_bystim.(f_stim).neuron_spec_freq = [];
        data_bystim.(f_stim).stim_timestamps = [];
        data_bystim.(f_stim).trace_timestamps = [];

        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles(17) % DEBUG just trying the one neuron first
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    
            cur_fov_Fs = [];
            cur_fov_subVm = [];
            cur_fov_phases = [];
            cur_fov_fourcoeff = [];
            cur_fov_stim_time = [];
            cur_fov_trace_time = [];
            cur_fov_filtsig = [];

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
                    if isempty(trial_data) || sum(isnan(cur_fov_subVm(:))) > 0
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

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
                    
                    % Save the phases for all trials
                    [wt, f] = get_power_spec(detrend_subVm', mean(cur_fov_Fs, 'omitnan'));                   
                    cur_fov_fourcoeff(:, :, end + 1) = wt;
                    
                    % Save the hilbert transform signal

                    [filt_sig] = filt_data(detrend_subVm', [1:1:150] , mean(cur_fov_Fs, 'omitnan'));
                    cur_fov_filtsig(:, :, end + 1) = filt_sig;

                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Plot the subthreshold heatmap, and then the coherence
            figure('Renderer', 'Painters', 'Position', [200 200 900 700]);
            tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            nexttile;
            surface(nanmean(cur_fov_trace_time, 2), 1:size(cur_fov_subVm, 2), cur_fov_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            Multi_func.set_default_axis(gca);
            xlabel('Time from onset(s)');
            ylabel('Trial #');
            xlim([-1, 2.05]);
            title('Sub Vm');
            
            nexttile;
            itc_map = abs(mean(cur_fov_fourcoeff./abs(cur_fov_fourcoeff), 3, 'omitnan'));
            %itc_map = exp(1i.*angle(cur_fov_fourcoeff));
            %itc_map = abs(mean(itc_map, 3, 'omitnan'));
            surface(nanmean(cur_fov_trace_time, 2), f, itc_map, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
            xlim([-1, 2.05]);
            xlabel('Time from onset(s)');
            ylabel('Frequency (Hz)');
            a = colorbar;    
            a.Label.String = 'ITC';
            Multi_func.set_default_axis(gca);
            title('ITC');

            nexttile;
            % One method has NaNs and the other has zeroes
            % Either one seems very close
            %itc_map_hil = abs(mean(cur_fov_filtsig./abs(cur_fov_filtsig), 3, 'omitnan'));
            itc_map_hil = exp(1i.*angle(cur_fov_filtsig));
            itc_map_hil = abs(mean(itc_map_hil, 3, 'omitnan'));

            surface(nanmean(cur_fov_trace_time, 2), [1:1:150], itc_map_hil, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
            xlim([-1, 2.05]);
            xlabel('Time from onset(s)');
            ylabel('Frequency (Hz)');
            a = colorbar;    
            a.Label.String = 'ITC';
            Multi_func.set_default_axis(gca);
            title('ITC from Hilbert');
            sgtitle(matfile, 'Interpreter', 'none');
            
            % Look at rose plot for 2.2 sec time point, 1817 idx
            %data_line = cur_fov_fourcoeff(35, 1817, :);
            %figure;
            %plot(squeeze(angle(data_line)));
            

            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_Vm;
            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
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

% Check if combining all of the regions or not
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

% Calculate the sampling frequency from all of the 
avg_Fs = nanmean(all_Fs);

% Filter data with hilbert transform
function  [filt_sig]=filt_data(sig,frs, FS)
    Fn = FS/2;

    for steps=1:length(frs);
        FB=[ frs(steps)*0.8 frs(steps)*1.2];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        filt_sig(steps, :)= hilbert(filtfilt(B,A,sig));
    end
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

