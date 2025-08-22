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

figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

exclude_200ms = 1;

%% Read in the saved pv data and perform analysis
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
    %save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'ca1_data.mat'];
end
%Load the data
load(save_all_data_file);

% -- Depracated --
% Check if combining all of the regions or not
%if all_regions == 1
%    region_data = Multi_func.combine_regions_old(region_data);
%end

if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

field1 = fieldnames(region_data);
field1 = field1(1);
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');

%% Time to peak SubVm before and after stimulation comparison for each trial
% Find the highest peak within 500ms before stimulation and compare to post stimulation
t_to_stim = 250 ; % in ms
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    stims = fieldnames(data_bystim);


    % Loop through stimulation frequencies
    for f_stim=stims'
        f_stim = f_stim{1};
        
        % Specify timepoints within trial
        timeline = nanmean(region_data.(f_region).(f_stim).trace_timestamps, 2);
        post_ran_idx = find(timeline >= 1 & timeline <= 1.5);
        base_idx = find(timeline < 0);

        % Build a filter

        neuron_trial_map = {};
        
        % Store each neuron's pre and post coefficient of variation
        pre_coef_var = [];
        post_coef_var = [];

        % Loop through each neuron
        for neuron=1:length(data_bystim.(f_stim).all_trial_SubVm)
            % Get the raw trials
            nr_subvm = data_bystim.(f_stim).all_trial_SubVm{neuron};
            
            % Time data for separating periods
            timeline = data_bystim.(f_stim).trace_timestamps(:, neuron);
            pulse_diff = nanmean(diff(data_bystim.(f_stim).stim_timestamps(:, neuron)));
            stim_offset_t = data_bystim.(f_stim).stim_timestamps(end, neuron) + pulse_diff;

            % Baseline to stim time 
            base_idx = find(timeline < 0 & timeline > 0 - t_to_stim/1000);
            % Post from offset time
            post_idx = find(timeline > stim_offset_t & timeline < stim_offset_t + t_to_stim/1000);
            
            % Store the times
            pre_stim_times = [];
            post_stim_times = [];

            %DEBUG print figure with peak times
            figure();
            tiledlayout(size(nr_subvm, 2), 3, 'TileSpacing', 'compact', 'Padding', 'compact');
            % Loop through each trial
            for tr=1:size(nr_subvm, 2)
                % Filter the subvm
                %TODO Fn = 

                % Peak pre stim
                [pks, locs, w, pre_p] = findpeaks(nr_subvm(base_idx, tr));
                pre_pk_loc = base_idx(locs(find(pre_p == max(pre_p))));
            
                % Peak post stim
                [pks, locs, w, post_p] = findpeaks(nr_subvm(post_idx, tr));
                post_pk_loc = post_idx(locs(find(post_p == max(post_p))));

                %DEBUG plot the trial Vm along with the closest peak height
                nexttile;
                plot(timeline, nr_subvm(:, tr));
                hold on;
                plot(timeline(pre_pk_loc), nr_subvm(pre_pk_loc, tr), 'or');
                hold on;
                plot(repmat(timeline(pre_pk_loc), 2, 1), ...
                     [nr_subvm(pre_pk_loc, tr), ...
                    nr_subvm(pre_pk_loc, tr) - pre_p(find(pre_p == max(pre_p)))], '-g');
                hold on;
                plot(timeline(post_pk_loc), nr_subvm(post_pk_loc, tr), 'or');
                hold on;
                plot(repmat(timeline(post_pk_loc), 2, 1), ...
                     [nr_subvm(post_pk_loc, tr), ...
                    nr_subvm(post_pk_loc, tr) - post_p(find(post_p == max(post_p)))], '-g');
                
                % Look at the distribution of prominensces
                nexttile;
                histogram(pre_p, 1000);
                nexttile;
                histogram(post_p, 1000);
                
                % END DEBUG

                pre_stim_times(end + 1) = 0 - timeline(pre_pk_loc);
                post_stim_times(end + 1) = timeline(post_pk_loc) - stim_offset_t ;
            end
            
            % Store the coefficient of variation for pre and post
            pre_coef_var(end + 1) = nanstd(pre_stim_times)/nanmean(pre_stim_times);
            post_coef_var(end + 1) = nanstd(post_stim_times)/nanmean(post_stim_times);
            
            legend(['Pre CV:' num2str(pre_coef_var(end)) ...
                    ' Post CV:' num2str(post_coef_var(end))]);

        end % End-- looping through neuron
        
        % Plot the violin
        figure();
        data = [pre_coef_var, post_coef_var];
        labels = [repmat({'Pre'}, 1, length(pre_coef_var)),...
            repmat({'Post'}, 1, length(post_coef_var))];
            
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Pre', 'Post'}, ViolinOpts);
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.post_color};
        Multi_func.set_default_axis(gca);
        hold on;
        plot([1, 2], [pre_coef_var(:), post_coef_var(:)]);

        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        
            
    end % End-- looping through stimulation
end % End-- looping through region
