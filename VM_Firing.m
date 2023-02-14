clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
server_root_path = '~/Projects/';
% Windows server
%server_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Smoothing parameter for spike rate
srate_win = 50;

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

%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%    %% Select matfiles by stim specific conditions for all regions
    [matfile_stim] = stim_cond(all_matfiles); %stim_cond(region_matfiles.(f_region).names);
    %% Loop through each field of the struct and concatenate everything together
    
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies
    all_Fs = [];

    % Loop through each stimulation condition
    for f_stim = fieldnames(matfile_stim)'
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    
    
        % Initialize field subthreshold array
        data_bystim.(f_stim) = struct();
        data_bystim.(f_stim).neuron_Vm = [];
        data_bystim.(f_stim).neuron_srate = [];
    
        % Loop through each matfile of the current stimulation condition
        for matfile = matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                cur_fov_subVm = [];
                cur_fov_srate = [];
                cur_fov_raster = [];

                % Store the camera framerate
                all_Fs(end+1) = trial_data.camera_framerate;
                
                % Loop through each trial                
                for tr_idx=trial_idxs        
                    trial_data = data.align.trial{tr_idx};
                    
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
    
                    % Grab the subthreshold Vm
                    % Chop the respective frames
                    cur_trace_ws = trial_data.spike_info.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                    cur_fov_subVm = horzcat_pad(cur_fov_subVm, cur_trace_ws');
                    
                    % Calculate the spike rate
                    cur_raster = trial_data.spike_info.roaster(roi_idx, front_frame_drop:back_frame_drop);
                    cur_spikerate = cur_raster.*trial_data.camera_framerate;
                    cur_spikerate = nanfastsmooth(cur_spikerate, srate_win, 1, 1);
                    cur_fov_srate = horzcat_pad(cur_fov_srate, cur_spikerate');            
                    
                    % Save the raster plot
                    cur_fov_raster = horzcat_pad(cur_fov_raster, cur_raster');

                    % Keep track of the figure being used
                    %if strcmp(ri{5}, '1000') == 1
                    %    figure('visible', 'off');
                    %    plot(cur_trace_ws)
                    %    title([matfile{1} ' ROI#' num2str(roi_idx) ' Trial ' num2str(tr_idx)], 'Interpreter', 'none');
                    %    saveas(gcf, ['debug/' matfile{1}(1:end-4) '.png']);
                    %    close gcf;
                    %end
                end

                % Plot each neuron's raster
                %figure('visible', 'off');
                %for i = 1:size(cur_fov_raster, 2)
                %    plot(cur_fov_raster(:, i)*i*0.2, '.');
                %    hold on;
                %end
                %ylim([0.1, (size(cur_fov_raster, 2)*0.2) + 1]);
                %title(['Raster of ' matfile{1}(1:end-4)], 'Interpreter', 'none');
                %saveas(gcf, [figure_path 'Raster plots' f matfile{1}(1:end-4) '.png']);
            end
            
            %EBUG
            %figure;
            %plot(nanmean(cur_fov_subVm, 2));
    
            % Average for each neuron and save the subthreshold Vm
            temp = data_bystim.(f_stim).neuron_Vm;
            data_bystim.(f_stim).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
            temp = data_bystim.(f_stim).neuron_srate;
            data_bystim.(f_stim).neuron_srate = horzcat_pad(temp, nanmean(cur_fov_srate, 2));

        end % End looping through FOVs of a condition
    end

%    % Save the VM to the specific region
%    region_data.(f_region) = data_bystim;

%end

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

%% Full collective subthreshold spectra
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
freqLimits = [0 150];

% Loop through each stimulation parameter
for stim=stims'
    
    % Loop through each neuron
    all_abs_wt = [];
    for i=1:size(data_bystim.(stim{1}).neuron_Vm, 2)
        cur_subVm = data_bystim.(stim{1}).neuron_Vm(:, i);
        
        % Check if there are nans in the trace
        if isnan(cur_subVm)
            continue;
        end
        fb = cwtfilterbank(SignalLength=length(cur_subVm),...
                           SamplingFrequency=avg_Fs,...
                           FrequencyLimits=freqLimits);
        [wt, f] = cwt(cur_subVm, FilterBank=fb);
        all_abs_wt = cat(3, all_abs_wt, abs(wt));
    end
    nexttile;
    surface(timeline, f, nanmean(all_abs_wt, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
    colorbar;
    title(stim{1}(3:end), 'Interpreter', 'none');
    %nexttile;
    %imagesc(timeline, f, abs(wt));
    %
    %%TODO change scale for showing the frequencies
    %yticks(flip(f([1, 12, 24, 64])));
    %yticklabels(string(flip(f([1, 12, 24, 64] ))));
    %title(stim, 'Interpreter', 'none');
end
sgtitle('Spectra from averaged Sub Vm Trace');

% Full collective spike rate over time
stims = fieldnames(data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for stim=stims'
    
    cur_srate = nanmean(data_bystim.(stim{1}).neuron_srate, 2);
    nexttile;
    plot(cur_srate);
    %nexttile;
    %imagesc(timeline, f, abs(wt));
    %
    %%TODO change scale for showing the frequencies
    %yticks(flip(f([1, 12, 24, 64])));
    %yticklabels(string(flip(f([1, 12, 24, 64] ))));
    title(stim{1}(3:end), 'Interpreter', 'none');
end
sgtitle('Average Spike rate');

figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for stim=stims'
    stim{1}
    cur_Vm = nanmean(data_bystim.(stim{1}).neuron_Vm, 2);
    nexttile;
    plot(cur_Vm);
    %nexttile;
    %imagesc(timeline, f, abs(wt));
    %
    %%TODO change scale for showing the frequencies
    %yticks(flip(f([1, 12, 24, 64])));
    %yticklabels(string(flip(f([1, 12, 24, 64] ))));
    title(stim{1}(3:end), 'Interpreter', 'none');
end
sgtitle('Average subthreshold Vm');


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
