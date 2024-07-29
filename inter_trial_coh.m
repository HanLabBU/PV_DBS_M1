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

% Flag for neurons that have the low frequency oscillations
sep_neurons_by_freq = 0;
use_low_freq = 1; % 1 to calculate itc for low freq neurons

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

%% Identify low frequency oscillations
if sep_neurons_by_freq == 1
    % Determine neurons that do or do not have low frequency oscillations
    for f_region = fieldnames(region_data)'
        f_region = f_region{1};
        for f_stim = fieldnames(region_data.(f_region))'
            f_stim = f_stim{1};
            region_data.(f_region).(f_stim).has_low_freq = Multi_func.get_low_freq_neurons(region_data.(f_region).(f_stim));
        end
    end
end


%% Calculate the ITC for all regions and stimulation frequencies
freq_idx = Multi_func.theta_range(1):Multi_func.theta_range(2); 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    stims = fieldnames(data_bystim);

    % Loop through stimulation frequencies
    for f_stim=stims'
        f_stim = f_stim{1};
        
        % Specify timepoints within trial
        timeline = nanmean(region_data.(f_region).(f_stim).trace_timestamps, 2);
        post_ran_idx = find(timeline >= 1 & timeline <= 1.5); %TODO should use a shared variable across whole script
        base_idx = find(timeline < 0);

        % Store each neuron's base and post stimulation ITC for comparison
        nr_base_itc = [];
        nr_post_itc = [];
        
        neuron_trial_map = {};

        % Loop through each neuron
        for neuron=1:length(data_bystim.(f_stim).neuron_hilbfilt)
            
            if sep_neurons_by_freq == 1 && data_bystim.(f_stim).has_low_freq(neuron) == ~use_low_freq
                disp('Skipped');
                continue;                
            end

            cur_fourcoeff = data_bystim.(f_stim).neuron_hilbfilt{neuron};
            itc_map = abs(mean(cur_fourcoeff./abs(cur_fourcoeff), 3, 'omitnan'));
 
            % Append neuron's ITC values to population wide
            nr_base_itc(end + 1) = mean( itc_map(freq_idx, base_idx), 'all');
            nr_post_itc(end + 1) = mean( itc_map(freq_idx, post_ran_idx), 'all');
            
            % Get the trial's Vm for current neuron
            neuron_trial_map{end + 1} = data_bystim.(f_stim).all_trial_rawVm{neuron};
        end        

        % Sort neurons based on the ITC value after stimulation period
        [~, sort_i] = sort(nr_post_itc);
        sort_i = flip(sort_i);
        sorted_trials = neuron_trial_map(sort_i);
        nr_base_itc = nr_base_itc(sort_i);
        nr_post_itc = nr_post_itc(sort_i);

        % Store the trial maps back into the structures
        region_data.(f_region).(f_stim).itc_sorted_neuron_trials_map = sorted_trials;
        region_data.(f_region).(f_stim).itc_neuron_baseline = nr_base_itc;
        region_data.(f_region).(f_stim).itc_neuron_post = nr_post_itc;
    end
end

%% Plot the trial heatmap for all brain regions and stimulation frequencies
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region); %

    stims = fieldnames(data_bystim);
    % Loop through stimulation frequencies
    for f_stim=stims'
        f_stim = f_stim{1};
        
        neuron_bound = [0.5];
        sorted_trials = data_bystim.(f_stim).itc_sorted_neuron_trials_map;
        %Loop through the trials to create the bounds
        for i=1:length(sorted_trials)
            neuron_bound(end + 1) = neuron_bound(end) + size(sorted_trials{i}, 2);
        end
        
        trial_map = [sorted_trials{:}];

        % Plot all of the neurons Vm mapping sorted by ITC
        figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.09 27.94]);
        
        imagesc('XData', timeline, 'YData', 1:size(trial_map, 2) , 'CData', trial_map');
        hold on;
        yline(neuron_bound);
        
        text_y = 0.5 + [neuron_bound(1:end - 1) + neuron_bound(2:end)]/2;
        % Plot the Post stim ITC values
        for i=1:length(data_bystim.(f_stim).itc_neuron_post)
            text(2.2, text_y(i), num2str(data_bystim.(f_stim).itc_neuron_post(i)));
        end
     
        % Indicate if plotting for those with low frequency neuron
        title_ext = '';
        if sep_neurons_by_freq == 1 && use_low_freq == 1
            title_ext = 'Low Frequency Neurons';
        elseif sep_neurons_by_freq == 1 && use_low_freq == 0
            title_ext = 'non-low freq neurons';
        end

        Multi_func.set_default_axis(gca);
        
        title([f_region ' ' f_stim ' ' title_ext], 'Interpreter', 'none');
        xlim([-1, 2.5]);
        
        saveas(gcf, [figure_path 'ITC/Vm_heatmap_itc_sorted_' title_ext '_' f_region '_' f_stim '.png']);
        saveas(gcf, [figure_path 'ITC/Vm_heatmap_itc_sorted_' title_ext '_' f_region '_' f_stim '.pdf']);

    end
end

%% Plot violins to comapre baseline and post across all regions and stimulation frequency
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    stims = fieldnames(region_data.(f_region))
    % Start the figures for each brain region
    figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [1, 1, 7, 12]);
    
    for f_stim = stims'
        f_stim = f_stim{1};
        nexttile;
        popul_data = region_data.(f_region).(f_stim);
                
        data = [popul_data.itc_neuron_baseline, popul_data.itc_neuron_post];
        labels = [repmat({'Base'}, 1, length(popul_data.itc_neuron_baseline)), repmat({'Post'}, 1, length(popul_data.itc_neuron_post))];
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Post'}, ViolinOpts);
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.post_color};
        hold on;
        
        % Conditionally color lines between violinplots
        for i=1:length(popul_data.itc_neuron_post)
            color = [0 0 0 0.4];
        
            % Neuron increased their ITC (considered ITC neuron)
            if popul_data.itc_neuron_post(i) > popul_data.itc_neuron_baseline(i)*1.20
                color = [[30, 2, 237]/255, 0.4];
                
            % Neuron had a decrease in their ITC
            elseif popul_data.itc_neuron_post(i) < popul_data.itc_neuron_baseline(i)*0.8
                color = [[235, 5, 28]/255, 0.4];
            end
        
            plot([1 2], [popul_data.itc_neuron_baseline(i), popul_data.itc_neuron_post(i)], '-', 'Color', color);
            hold on
        end

        % Setting the axis default
        Multi_func.set_default_axis(gca);
    
        % Indicate if plotting for those with low frequency neuron
        title_ext = '';
        if sep_neurons_by_freq == 1 && use_low_freq == 1
            title_ext = 'Low Frequency Neurons';
        elseif sep_neurons_by_freq == 1 && use_low_freq == 0
            title_ext = 'non-low freq neurons';
        end
        
        title([f_region ' ' f_stim], 'Interpreter', 'none');
    end
    sgtitle(['Base to Post ITC Comparison' title_ext ' neurons']);

    saveas(gcf, [figure_path 'ITC/Base_vs_Post_ITC_for_' title_ext '_' f_region '.png']);
    saveas(gcf, [figure_path 'ITC/Base_vs_Post_ITC_for_' title_ext '_' f_region '.pdf']);

end

%TODO perform statistical test between each group for each condition

%% Calculate shuffled distribution of offset times with ITC
num_shuf = 500;
wind_oi = 500;
% Ensure to change the offset time between onset to end of trial-(itc period calculation)
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    stims = fieldnames(data_bystim);

    % Loop through stimulation frequencies
    for f_stim=stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        figure('Position', [0, 0, 750, 2000]);
        tiledlayout(length(popul_data.trial_num), 2, 'TileSpacing', 'tight', 'Padding', 'tight');

        % Loop through each neuron
        for nr=1:length(popul_data.trial_num)
            timeline = popul_data.trace_timestamps(:, nr);
            post_stim_idx = find(timeline >= 1 & timeline <= 1.5);
            
            nr_coeffs = popul_data.neuron_hilbfilt{nr};

            % Shuffled post-stim ITC
            shuf_itc = [];
            
            % loop through and shuffle the offset time for each trial and calculate the post_itc
            for i=1:num_shuf
                shuf_tr_coeff = [];
                %loop through each trial
                for tr=1:popul_data.trial_num(nr)
                    % Grab the stim onset to max post-stim offset period
                    stim_onset_range = [find(timeline == 0), length(timeline) - wind_oi ];
                    rand_start = randi(stim_onset_range, 1);
                    rand_range = rand_start:rand_start + wind_oi;

                    shuf_tr_coeff(:, :, end+1) = nr_coeffs(:, rand_range, tr);
                end
                
                % Remove the initial trial
                shuf_tr_coeff(:, :, 1) = [];
                
                % Calculate shuffled trial ITC
                itc_wind_map = abs(mean(shuf_tr_coeff./abs(shuf_tr_coeff), 3, 'omitnan'));
                shuf_itc(end + 1) = mean(itc_wind_map(freq_idx, :), 'all');
            end
            
            % Caluclate the acual post-stim ITC
            itc_map = abs(mean(nr_coeffs./abs(nr_coeffs), 3, 'omitnan'));
            obs_post_itc = mean(itc_map(freq_idx, post_stim_idx), 'all');

            % Plot the shuffled distribution
            nexttile;
            high_prc = prctile(shuf_itc(:), 97.5);
            low_prc = prctile(shuf_itc(:), 2.5);
            
            histogram(shuf_itc, 1000);
            hold on;
            xline([low_prc, high_prc], '-b');
            hold on;
            
            Multi_func.set_default_axis(gca);

            % If observed is higher make the line green, otherwise
            %within distribution just make it red
            if obs_post_itc > high_prc
                xline(obs_post_itc, '-g');
            else
                xline(obs_post_itc, '-r');
            end

            % Plot the subvm heatmap
            nexttile;
            imagesc('XData', timeline, 'YData', 1:popul_data.trial_num(nr), 'CData',...
                popul_data.all_trial_rawVm{nr}');

            Multi_func.set_default_axis(gca);
        end

        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
        saveas(gcf, [figure_path 'ITC/Shuffled_histogram_' f_region '_' f_stim '.png']);
    end
end



%% Calculate all of the ITCs for 40Hz M1 region neurons, and sort by the post ITRC  
f_region = 'r_M1';
f_stim = 'f_40';
timeline = nanmean(region_data.(f_region).(f_stim).trace_timestamps, 2);
data_bystim = region_data.(f_region);
% Save the index of the neurons with the ITC phenomena
itc_neuron = [];
post_ran_idx = find(timeline >= 1 & timeline <= 1.5);
base_idx = find(timeline < 0);
freq_idx = Multi_func.theta_range(1):Multi_func.theta_range(2); 

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

% Sort neurons based on the ITC value after stimulation period
[~, sort_i] = sort(nr_post_itc);
sort_i = flip(sort_i);
sorted_trials = neuron_trial_map(sort_i);
nr_base_itc = nr_base_itc(sort_i);
nr_post_itc = nr_post_itc(sort_i);


%% Plot all of the neuron's in sorting post-stim ITC order
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

%% Plot violins to compare the base vs post stim ITC values for 40Hz motor cortex
figure('Renderer', 'Painters', 'Position', [200 200 700 700]);
data = [nr_base_itc, nr_post_itc];
labels = [repmat({'Base'}, 1, length(nr_base_itc)), repmat({'Post'}, 1, length(nr_post_itc))];
ViolinOpts = Multi_func.get_default_violin();
violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Post'}, ViolinOpts);
violins(1).ViolinColor = {Multi_func.base_color};
violins(2).ViolinColor = {Multi_func.post_color};

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
ax = gca;
ax.Units = 'Centimeters';
ax.InnerPosition = [5 5 4.0 6.43];
Multi_func.set_default_axis(gca);
ylabel('ITC Value');
title('All 40Hz Neurons ITC Base vs Stim comparison');

saveas(gcf, [figure_path 'ITC' f '40Hz_base_vs_stim.pdf']);
saveas(gcf, [figure_path 'ITC' f '40Hz_base_vs_stim.png']);

%% Perform statistical test
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

%% Calculate all of the ITCs for 140Hz M1 region neurons, and sort by the post ITRC  
f_region = 'r_M1';
f_stim = 'f_140';
timeline = nanmean(region_data.(f_region).(f_stim).trace_timestamps, 2);
data_bystim = region_data.(f_region);
% Save the index of the neurons with the ITC phenomena
itc_neuron = [];
post_ran_idx = find(timeline >= 1 & timeline <= 1.5);
base_idx = find(timeline < 0);
freq_idx = Multi_func.theta_range(1):Multi_func.theta_range(2); 

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

% Sort neurons based on the ITC value after stimulation period
[~, sort_i] = sort(nr_post_itc);
sort_i = flip(sort_i);
sorted_trials = neuron_trial_map(sort_i);
nr_base_itc = nr_base_itc(sort_i);
nr_post_itc = nr_post_itc(sort_i);


%% Plot violins to compare the base vs post stim ITC values for all neurons
figure('Renderer', 'Painters', 'Position', [200 200 700 700]);
data = [nr_base_itc, nr_post_itc];
labels = [repmat({'Base'}, 1, length(nr_base_itc)), repmat({'Post'}, 1, length(nr_post_itc))];
ViolinOpts = Multi_func.get_default_violin();
violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Post'}, ViolinOpts);
violins(1).ViolinColor = {Multi_func.base_color};
violins(2).ViolinColor = {Multi_func.post_color};

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
ax = gca;
ax.Units = 'Centimeters';
ax.InnerPosition = [5 5 4.0 6.43];
Multi_func.set_default_axis(gca);
ylabel('ITC Value');
title('All 140Hz Neurons ITC Base vs Stim comparison');

saveas(gcf, [figure_path 'ITC' f '140Hz_base_vs_stim.pdf']);
saveas(gcf, [figure_path 'ITC' f '140Hz_base_vs_stim.png']);

%% Perform statistical test
[p, h, stats] = signrank(nr_base_itc, nr_post_itc)
signrank_itc_base_post_stats = struct();
signrank_itc_base_post_stats.p = p
signrank_itc_base_post_stats.h = h;
signrank_itc_base_post_stats.stats = stats;

%%
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
