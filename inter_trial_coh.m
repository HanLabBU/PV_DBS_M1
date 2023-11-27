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

set(0,'DefaultFigureVisible','off');

%%% END Modification

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
%for f_region = {'r_M1'} %fieldnames(region_matfiles)'
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
%    % Loop through 40Hz stimulated neurons
%    for f_stim = {'f_40'}
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
%        for matfile = matfiles(17) % DEBUG just trying the one neuron first
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile{1}]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%            cur_fov_Fs = [];
%            cur_fov_subVm = [];
%            cur_fov_phases = [];
%            cur_fov_fourcoeff = [];
%            cur_fov_stim_time = [];
%            cur_fov_trace_time = [];
%            cur_fov_filtsig = [];
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
%                    
%                    % Save the phases for all trials
%                    [wt, f] = get_power_spec(detrend_subVm', mean(cur_fov_Fs, 'omitnan'));                   
%                    cur_fov_fourcoeff(:, :, end + 1) = wt;
%                    
%                    % Save the hilbert transform signal
%
%                    [filt_sig] = filt_data(detrend_subVm', [1:1:150] , mean(cur_fov_Fs, 'omitnan'));
%                    cur_fov_filtsig(:, :, end + 1) = filt_sig;
%
%                end % End looping through each neuron
%            end
%            
%            % Skip rest of the calculations if the subthreshold Vm is nan
%            if sum(isnan(cur_fov_subVm(:))) > 0 || isempty(cur_fov_subVm)
%                continue;
%            end
%
%            % Plot the subthreshold heatmap, and then the coherence
%            figure('Renderer', 'Painters', 'Position', [200 200 900 700]);
%            tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
%            nexttile;
%            surface(nanmean(cur_fov_trace_time, 2), 1:size(cur_fov_subVm, 2), cur_fov_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%            Multi_func.set_default_axis(gca);
%            xlabel('Time from onset(s)');
%            ylabel('Trial #');
%            xlim([-1, 2.05]);
%            title('Sub Vm');
%            
%            nexttile;
%            itc_map = abs(mean(cur_fov_fourcoeff./abs(cur_fov_fourcoeff), 3, 'omitnan'));
%            %itc_map = exp(1i.*angle(cur_fov_fourcoeff));
%            %itc_map = abs(mean(itc_map, 3, 'omitnan'));
%            surface(nanmean(cur_fov_trace_time, 2), f, itc_map, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
%            xlim([-1, 2.05]);
%            xlabel('Time from onset(s)');
%            ylabel('Frequency (Hz)');
%            ylim([0 20]);
%            a = colorbar;    
%            a.Label.String = 'ITC';
%            Multi_func.set_default_axis(gca);
%            title('ITC');
%
%            %nexttile;
%            %% One method has NaNs and the other has zeroes
%            %% Either one seems very close
%            %%itc_map_hil = abs(mean(cur_fov_filtsig./abs(cur_fov_filtsig), 3, 'omitnan'));
%            %itc_map_hil = exp(1i.*angle(cur_fov_filtsig));
%            %itc_map_hil = abs(mean(itc_map_hil, 3, 'omitnan'));
%
%            %surface(nanmean(cur_fov_trace_time, 2), [1:1:150], itc_map_hil, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
%            %xlim([-1, 2.05]);
%            %xlabel('Time from onset(s)');
%            %ylabel('Frequency (Hz)');
%            %a = colorbar;    
%            %a.Label.String = 'ITC';
%            %Multi_func.set_default_axis(gca);
%            %title('ITC from Hilbert');
%            sgtitle(matfile, 'Interpreter', 'none');
%            
%            % Look at rose plot for 2.2 sec time point, 1817 idx
%            %data_line = cur_fov_fourcoeff(35, 1817, :);
%            %figure;
%            %plot(squeeze(angle(data_line)));
%            
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
if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

disp('Finished Loading');

% Set to look for 
f_region = 'r_M1';
f_stim = 'f_40';
timeline = nanmean(region_data.(f_region).(f_stim).trace_timestamps, 2);
data_bystim = region_data.(f_region);
% Save the index of the neurons with the ITC phenomena
itc_neuron = [];
post_ran_idx = find(timeline >= 1 & timeline <= 1.5);
base_idx = find(timeline < 0);
freq_idx = 4:12; % Check the theta band %TODO maybe change to 2:12? Need to determine a good theta range

% Store each neuron's base and post stimulation ITC for comparison
nr_base_itc = [];
nr_post_itc = [];

neuron_trial_map = {};
% Sort all of the neurons based on the post ITC value
for neuron = 1:length(data_bystim.(f_stim).neuron_hilbfilt)
    cur_fourcoeff = data_bystim.(f_stim).neuron_hilbfilt{neuron};
    itc_map = abs(mean(cur_fourcoeff./abs(cur_fourcoeff), 3, 'omitnan'));
 
    % Append neuron's ITC values to population wide
    nr_base_itc(end + 1) = mean( itc_map(freq_idx, base_idx), 'all');
    nr_post_itc(end + 1) = mean( itc_map(freq_idx, post_ran_idx), 'all');
    
    % Get the trial's Vm for current neuron
    neuron_trial_map{end + 1} = data_bystim.(f_stim).all_trial_rawVm{neuron};
end

[~, sort_i] = sort(nr_post_itc);
sort_i = flip(sort_i);
sorted_trials = neuron_trial_map(sort_i);
nr_base_itc = nr_base_itc(sort_i);
nr_post_itc = nr_post_itc(sort_i);

neuron_bound = [0];
%Loop through the trials to create the bounds
for i=1:length(sorted_trials)
    neuron_bound(end + 1) = neuron_bound(end) + size(sorted_trials{i}, 2);
end

trial_map = [sorted_trials{:}];
%trial_map = reshape(trial_map, length(timeline), []);

% Plot the sorted neurons
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.09 27.94]);
surface(timeline, 1:size(trial_map, 2), trial_map', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
hold on;
yline(neuron_bound + 0.5);

text_y = 0.5 + [neuron_bound(1:end - 1) + neuron_bound(2:end)]/2;
% Plot the Post stim ITC values
for i=1:length(nr_post_itc)
    text(2.2, text_y(i), num2str(nr_post_itc(i)));
end
title('Sorted Trial Map');
Multi_func.set_default_axis(gca);
xlim([-1, 2.5]);

% Plot violins to compare the base vs post stim ITC values for all neurons
figure('Renderer', 'Painters', 'Position', [200 200 700 700]);
data = [nr_base_itc, nr_post_itc];
labels = [repmat({'Base'}, 1, length(nr_base_itc)), repmat({'Post'}, 1, length(nr_post_itc))];
violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Post'});
hold on;
% Conditionally color lines between violinplots
for i=1:length(nr_base_itc)
    color = [0 0 0 0.4];

    % Neuron increased their ITC (considered ITC neuron)
    if nr_post_itc(i) > nr_base_itc(i)*1.20
        color = [[30, 2, 237]/255, 0.4];
        
    % Neuron had a decrease in their ITC
    elseif nr_post_itc(i) < nr_base_itc(i)*0.8
        color = [[235, 5, 28]/255, 0.4];
    end

    plot([1 2], [nr_base_itc(i), nr_post_itc(i)], '-', 'Color', color);
    hold on;
end
Multi_func.set_default_axis(gca);
ylabel('ITC Value');
title('All 40Hz Neurons ITC Base vs Stim comparison');
saveas(gcf, [figure_path 'ITC/40Hz_base_vs_stim.png']);

% Perform statistical test
[p, h, stats] = signrank(nr_base_itc, nr_post_itc)
signrank_itc_base_post_stats = struct();
signrank_itc_base_post_stats.p = p
signrank_itc_base_post_stats.h = h;
signrank_itc_base_post_stats.stats = stats;

% Loop through each itc neuron and create a large ITC map
itc_four = [];
itc_trials = []; 
itc_neuron = find(nr_post_itc > 0.5); % Finding criteria where the coherence is >0.5
itc_neuron = sort_i(itc_neuron);
for neuron = itc_neuron
    cur_fourcoeff = data_bystim.(f_stim).neuron_hilbfilt{neuron};
    itc_four = cat(3, itc_four, cur_fourcoeff);
    itc_trials = [itc_trials; data_bystim.(f_stim).all_trial_SubVm{neuron}'];

    % Print out the ITC neurons source name
    disp(data_bystim.(f_stim).neuron_name{neuron});
end

% Calculate the ITC for all trails across neurons that show an above 0.5 ITC
itc_map = abs(mean(itc_four./abs(itc_four), 3, 'omitnan'));

figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 7.0658, 6.36]);

% Display all of the trials
nexttile;
surface(timeline, 1:size(itc_trials, 1), itc_trials, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
Multi_func.set_default_axis(gca);
xlim([-1, 2.05]);
xlabel('Time from onset (s)');
ylabel('Trace #');
title('Collective Traces', 'Interpreter', 'none');
a = colorbar;    
a.Label.String = 'Vm';

ax2 = nexttile;
surface(timeline, [1:1:20], itc_map(1:20, :), 'CDataMapping',  'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
colormap(ax2, Multi_func.warm_cold_color);
Multi_func.set_default_axis(gca);
xlim([-1, 2.05]);
xlabel('Time from onset(s)');
ylabel('Frequency (Hz)');
%ylim([0 20]);
a = colorbar;    
a.Label.String = 'ITC';
Multi_func.set_default_axis(gca);
title('Collective ITC', 'Interpreter', 'none');
saveas(gcf, [figure_path 'ITC/40Hz_above_half_itc.png']);
saveas(gcf, [figure_path 'ITC/40Hz_above_half_itc.pdf']);

%saveas(gcf, [figure_path 'ITC/40Hz_all_rest_itc.pdf']);

% Using individual criteria for selecting neurons to do the final ITC
% Loop through each neuron
itc_neuron = [];
for neuron = 1:length(data_bystim.(f_stim).neuron_hilbfilt)
    cur_fourcoeff = data_bystim.(f_stim).neuron_hilbfilt{neuron};
    itc_map = abs(mean(cur_fourcoeff./abs(cur_fourcoeff), 3, 'omitnan'));
 
    % Append neuron's ITC values to population wide
    nr_base_itc(end + 1) = mean( itc_map(freq_idx, base_idx), 'all');
    nr_post_itc(end + 1) = mean( itc_map(freq_idx, post_ran_idx), 'all');

    % Check if the in phase theta reset is high for this neuron
    if mean( itc_map(freq_idx, post_ran_idx), 'all') > 1.20*mean( itc_map(freq_idx, base_idx), 'all') % & mean( itc_map(freq_idx, post_ran_idx), 'all') > 0.7
        itc_neuron(end + 1) = neuron;
        
        %figure('Renderer', 'Painters', 'Position', [200 200 700 700]);
        %tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        %
        %% Display all of the trials
        %nexttile;
        %cur_fov_subVm = data_bystim.(f_stim).all_trial_rawVm{neuron};
        %surface(timeline, 1:size(cur_fov_subVm, 2), cur_fov_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        %nexttile;
        %surface(timeline, [1:1:150], itc_map, 'CDataMapping',  'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
        %xlim([-1, 2.05]);
        %xlabel('Time from onset(s)');
        %ylabel('Frequency (Hz)');
        %%ylim([0 20]);
        %a = colorbar;    
        %a.Label.String = 'ITC';
        %Multi_func.set_default_axis(gca);
        %title(data_bystim.(f_stim).neuron_name{neuron}, 'Interpreter', 'none');
        %% Save each ITC maps
        %disp('Saving Time');
        %tic
        %saveas(gcf, [figure_path 'ITC/' data_bystim.(f_stim).neuron_name{neuron} '_ITC.png']);
        %toc
    end

%saveas(gcf, [figure_path 'ITC/' f_region '_' f_stim '_neuronwise_ITC.pdf']);
end

% Plot the neuron's trials that meet the 1.20*base criteria
crit_trials = neuron_trial_map(itc_neuron);
trial_map = [crit_trials{:}];
figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
surface(timeline, 1:size(trial_map, 2), trial_map', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
hold on;
%%yline(neuron_bound + 0.5);
%%
%%text_y = 0.5 + [neuron_bound(1:end - 1) + neuron_bound(2:end)]/2;
%%% Plot the Post stim ITC values
%%for i=1:length(nr_post_itc)
%%    text(2.2, text_y(i), num2str(nr_post_itc(i)));
%%end
title('Trials that showed criteria');
Multi_func.set_default_axis(gca);
xlim([-1, 2.5]);

% 
%% Each region ITC
%for f_region = fieldnames(region_data)'
%    f_region = f_region{1};
%    data_bystim = region_data.(f_region);
%    stims = fieldnames(data_bystim);
%    
%    % Loop through stimulation
%    for f_stim=stims'
%        f_stim = f_stim{1};
%        timeline = nanmean(data_bystim.(f_stim).trace_timestamps, 2);
%        
%        % Loop through each neuron
%        for neuron = 1:length(data_bystim.(f_stim).neuron_hilbfilt)
%            cur_fourcoeff = data_bystim.(f_stim).neuron_hilbfilt{neuron};
%            itc_map = abs(mean(cur_fourcoeff./abs(cur_fourcoeff), 3, 'omitnan'));
%            figure('Renderer', 'Painters', 'Position', [200 200 700 700]);
%            disp('n');
%            surface(timeline, [1:1:150], itc_map, 'CDataMapping',  'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none' );
%            xlim([-1, 2.05]);
%            xlabel('Time from onset(s)');
%            ylabel('Frequency (Hz)');
%            %ylim([0 20]);
%            a = colorbar;    
%            a.Label.String = 'ITC';
%            Multi_func.set_default_axis(gca);
%            title(data_bystim.(f_stim).neuron_name{neuron}, 'Interpreter', 'none');
%        % Save each ITC maps
%        saveas(gcf, [figure_path 'ITC/' data_bystim.(f_stim).neuron_name{neuron} '_ITC.png']);
%        %saveas(gcf, [figure_path 'ITC/' f_region '_' f_stim '_neuronwise_ITC.pdf']);
%        end
%        %sgtitle([f_region '_' f_stim '_ITC'], 'Interpreter', 'none');
%    end
%end

% Filter data with hilbert transform
function  [filt_sig]=filt_data(sig,frs, FS)
    Fn = FS/2;

    for steps=1:length(frs);
        FB=[ frs(steps)*0.8 frs(steps)*1.2];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        filt_sig(steps, :)= hilbert(filtfilt(B,A,sig));
    end
end
