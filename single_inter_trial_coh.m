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

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

exclude_200ms = 1;

%%% END Modification

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else
    front_frame_drop = 15 + round((828*.200));
end

back_frame_drop = 2496;

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
    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);

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
        for matfile = matfiles([18, 19]) % use 18, 19, 22, or  % DEBUG just trying the one neuron first
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
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(matfile{1}, '_');
                try
                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                catch
                    trial_ignr_list = [];
                end
        
                % Remove ignored trials from trial idx list
                trial_idxs = setdiff(trial_idxs, trial_ignr_list);

                % Skip ROI if there are at most 2 trials
                if length(trial_idxs) <= 2
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

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);
                    cur_fov_stim_time = horzcat_pad(cur_fov_stim_time, raw_trial_data.raw_stimulation_time - stim_start);
                    cur_fov_trace_time = horzcat_pad(cur_fov_trace_time, trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start);
                    
                    % Save the phases for all trials
                    [wt, freq] = Multi_func.get_power_spec(detrend_subVm', mean(cur_fov_Fs, 'omitnan'));                   
                    cur_fov_fourcoeff(:, :, end + 1) = wt;
                    
                    % Save the hilbert transform signal

                    [filt_sig] = Multi_func.filt_data(detrend_subVm', [1:1:150] , mean(cur_fov_Fs, 'omitnan'));
                    cur_fov_filtsig(:, :, end + 1) = filt_sig;

                end % End looping through each neuron
            end
            
            % Skip rest of the calculations if the subthreshold Vm is nan
            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
                continue;
            end

            % Plot the subthreshold heatmap, and then the coherence
            figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 20.59 27.94]);
            tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 5.0, 6.0]);
            ax1 = nexttile;
            surface(nanmean(cur_fov_trace_time, 2), 1:size(cur_fov_subVm, 2), cur_fov_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            %colormap(ax1, Multi_func.green_warm_cold_color);
            Multi_func.set_default_axis(gca);
            y = gca; y = y.YAxis;
            Multi_func.set_spacing_axis(y, 3, 1);
            a = colorbar;    
            a.Label.String = 'Vm';
            xlabel('Time from onset(s)');
            ylabel('Trial #');
            xlim([-0.8, 2.05]);
            title('Sub Vm');
            
            ax2 = nexttile;
            % Creates smooth ITC
            itc_map = abs(mean(cur_fov_fourcoeff./abs(cur_fov_fourcoeff), 3, 'omitnan'));
            
            % Creates ridged ITC plot
            %itc_map = exp(1i.*angle(cur_fov_fourcoeff));
            %itc_map = abs(mean(itc_map, 3, 'omitnan'));

            surface(nanmean(cur_fov_trace_time, 2), freq, itc_map, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
            colormap(ax2, Multi_func.warm_cold_color);
            xlim([-0.8, 2.05]);
            xlabel('Time from onset(s)');
            ylabel('Frequency (Hz)');
            ylim([0 20]);
            a = colorbar;    
            a.Label.String = 'ITC';
            Multi_func.set_default_axis(gca);
            y = gca; y = y.YAxis;
            Multi_func.set_spacing_axis(y, 5, 1);
            title('ITC');
            fontsize(gcf, 7, "points");

            %nexttile;
            %% One method has NaNs and the other has zeroes
            %% Either one seems very close
            %%itc_map_hil = abs(mean(cur_fov_filtsig./abs(cur_fov_filtsig), 3, 'omitnan'));
            %itc_map_hil = exp(1i.*angle(cur_fov_filtsig));
            %itc_map_hil = abs(mean(itc_map_hil, 3, 'omitnan'));

            %surface(nanmean(cur_fov_trace_time, 2), [1:1:150], itc_map_hil, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
            %xlim([-1, 2.05]);
            %xlabel('Time from onset(s)');
            %ylabel('Frequency (Hz)');
            %a = colorbar;    
            %a.Label.String = 'ITC';
            %Multi_func.set_default_axis(gca);
            %title('ITC from Hilbert');
            sgtitle(matfile, 'Interpreter', 'none');
            saveas(gcf, [figure_path 'ITC' f matfile{1} '_ITC.pdf']);
            saveas(gcf, [figure_path 'ITC' f matfile{1} '_ITC.png']);

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
            [wt, freq] = Multi_func.get_power_spec(nanmean(cur_fov_subVm, 2)', nanmean(cur_fov_Fs));
            temp = data_bystim.(f_stim).neuron_spec_power;
            data_bystim.(f_stim).neuron_spec_power = cat(3, temp, wt);
            temp = data_bystim.(f_stim).neuron_spec_freq;   
            data_bystim.(f_stim).neuron_spec_freq = cat(3, temp, freq);

        end % End looping through FOVs of a condition
    end

    % Save the VM to the specific region
    region_data.(f_region).data_bystim = data_bystim;
end

% Check if combining all of the regions or not
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

%% Calculate the sampling frequency from all of the 
%field1 = fieldnames(region_data);
%field1 = field1{1};
%avg_Fs = mean(region_data.(field1).f_40.framerate, 'omitnan');
