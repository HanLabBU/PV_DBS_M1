% Plot the currents used in the experiments
clear all;
close all;
clc;

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
savepath = Multi_func.save_plot;

exclude_200ms = 1;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

%% Load all region data

% Use congregated data to save all of the currents
% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end
tic;
%Load the data
load(save_all_data_file);
toc
disp('Finished Loading Data');
%% TODO plot these as days from surgery
%% TODO plot currents of each brain region and frequency

%%TODO plot the currents for each mouse, and try to indicate brain region or frequency


%% Plotting current from days to surgery for each neuron
% Red: 140Hz
% Blue: 40Hz
% Circle: M1
% Plus: V1

%TODO have a map that keeps track of how many repetitions there are
%of number of days since surgery and amperage
days_amp_map = [];
mouse_color = struct();

figure('Position', [0 0 1000 800]);
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Keep track of neurons plotted
        neurons_plotted = [];

        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            if ismember(nr, neurons_plotted)
                continue;
            end

            tokens = regexp(popul_data.neuron_name{nr}, '_', 'split');

            % Find all of the neurons that come from the same mouse
            same_m_nrs = find(contains(popul_data.neuron_name, tokens{1}) == 1);
            same_m_nrs_pts = [];
            neurons_plotted = cat(2, neurons_plotted, same_m_nrs);
            for i=same_m_nrs

                tokens = regexp(popul_data.neuron_name{i}, '_', 'split');

                % Calculate the days between imaging and surgery date
                surg_date = popul_data.surgery_date.(['m_' tokens{1}]);
                if isempty(surg_date)
                    continue;
                end
                tokens{3} = erase(tokens{3}, 'rec');

                surg_date = datetime(surg_date, 'InputFormat', 'yyyyMMdd');
                img_date = datetime(tokens{3}, 'InputFormat', 'yyyyMMdd');

                days_from_surg = days(img_date - surg_date);

                % Parse out the current amperage
                fov_i = find(contains(tokens, 'FOV') == 1);
                amp_i = fov_i + 2;
                if contains(tokens{amp_i}, '.') == 1
                    amp_str = regexp(tokens{amp_i}, '\.', 'split');
                    amp_str = amp_str{1};
                else
                    amp_str = tokens{amp_i};
                end

                % Set labels based on brain region and stim frequency
                plot_str = '';
                if strcmp(f_region, 'r_M1') == 1
                    plot_str(end+1) = 'o';
                elseif strcmp(f_region, 'r_V1') == 1
                    plot_str(end+1) = '+';
                end

                % Set the frequency stuff
                if strcmp(f_stim, 'f_40') == 1
                    plot_str(end+1) = 'b';
                elseif strcmp(f_stim, 'f_140') == 1
                    plot_str(end+1) = 'r';
                end

                % Increment the marker size for repeated points
                if size(days_amp_map, 1) < days_from_surg | size(days_amp_map, 2) < str2num(amp_str)
                    days_amp_map(days_from_surg, str2num(amp_str)) = 1;
                else
                    days_amp_map(days_from_surg, str2num(amp_str)) = ...
                        1 + days_amp_map(days_from_surg, str2num(amp_str));
                end

                % Plot data point
                hold on;
                plot(days_from_surg, str2num(amp_str), plot_str, 'MarkerSize', 8 + days_amp_map(days_from_surg, str2num(amp_str)));
                hold on;

                % Keep track of points
                same_m_nrs_pts(:, end + 1) = [days_from_surg; str2num(amp_str)];
            end

            % Plot lines across all FOVs of the same mouse
            %TODO sort the x values from decreasing to increasing, makes figure easier to read
            if ~isempty(same_m_nrs_pts)
                
                [~, i] = sort(same_m_nrs_pts(1, :));
                if isfield(mouse_color, ['m_' tokens{1}])
                    color = mouse_color.(['m_' tokens{1}]);
                    plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i), '-', 'Color', color);
                else
                    plot_h = plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i),  '-');
                    mouse_color.(['m_' tokens{1}]) = get(plot_h, 'Color');
                end
            end
        end
    end
end


% Add legend label for each item
text(400, 600, '+: V1 region');
hold on;
text(400, 580, 'o: M1 region');
hold on;
text(400, 560, '140Hz', 'Color', 'r');
hold on;
text(400, 540, '40Hz', 'Color', 'b');
hold on;

y_val = 600;
for mouse_f=fieldnames(mouse_color)'
    mouse_f = mouse_f{1};
    text(450, y_val, '--', 'Color', mouse_color.(mouse_f), 'FontWeight', 'bold');
    y_val = y_val - 20;
end
% Label 
xlabel('Days since surgery (days)');
ylabel('Amperage (uamp)');

saveas(gcf, [savepath 'Current/' '2d_currents_by_region_frequency.png']);
saveas(gcf, [savepath 'Current/' '2d_currents_by_region_frequency.pdf']);

%TODO make a 3D plot so that each mouse's data is in a different axis
% Try to reuse as much as possible from the top figure plotting
%% Generates a 3D plot of the mouse current stuff

%%

all_currents = [];
all_Fs = [];

% Loop through and grab each neuron's current amplitude
for f_region = {'r_M1'}
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    for f_stim = fieldnames(data_bystim)'
        f_stim = f_stim{1};

        % Loop through each neuron
        for nr = 1:length(data_bystim.(f_stim).neuron_name)
            cur_nr_name = strsplit(data_bystim.(f_stim).neuron_name{nr}, '_');
            
            % Get the current from the matfile
            f_mouse = ['m_' cur_nr_name{1}];
            f_rec = ['r_' cur_nr_name{3}];               
            f_neuron= ['n_' cur_nr_name{4}];
            f_stim = ['f_' cur_nr_name{5}];

            neuron_num = str2num(erase(cur_nr_name{4}, 'FOV'));
            
            [region_data] = create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
            
            % Saving the current at neuron number index
            region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(cur_nr_name{6});
            
        end
    end
end

%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%
%    % Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    % Select matfiles by stim condition for given region
%    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);
%
%    % Loop through each stimulation condition
%    for f_stim = fieldnames(matfile_stim)'
%        f_stim = f_stim{1};
%        matfiles = matfile_stim.(f_stim).names;    
%    
%        % Initialize field subthreshold array
%        data_bystim.(f_stim) = struct();
%
%        % Loop through each matfile of the current stimulation condition
%        for matfile = matfiles
%            matfile = matfile{1};
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%         
%            
%            % Loop through each ROI
%            for roi_idx=1:size(trial_data.detrend_traces, 2)
%                %Determine whether this roi is to be ignored for this particular trial
%                ri = strsplit(matfile, '_');
%                try
%                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
%                catch
%                    trial_ignr_list = [];
%                end
%                
%                % Remove trials from trial idx list
%                trial_idxs = setdiff(trial_idxs, trial_ignr_list);
%
%                % Skip if there is at most 2 trials
%                if length(trial_idxs) <= 2
%                    continue;
%                end
%                
%                % Get the current from the matfile
%                f_mouse = ['m_' ri{1}];
%                f_rec = ['r_' ri{3}];               
%                f_neuron= ['n_' ri{4}];
%                f_stim = ['s_' ri{5}];
%                
%                % convert FOV number to single number
%                neuron_num = str2num(erase(ri{4}, 'FOV'));
%
%                [region_data] = create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
%                
%                % Saving the current at neuron number index
%                region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(ri{6});
%
%                % Save all of the currents
%                all_currents(end + 1) = str2num(ri{6});
%
%            end
%        end
%    end
%end



% Clear all zeros in the current array
all_currents(find(all_currents == 0)) = [];

disp(['Average currents ' num2str(mean(all_currents))]);
disp(['Std currents ' num2str(std(all_currents))]);

% Plot all of the currents as a single plot
figure;
tick_labels = {};
tick_spots = [];
i = 1;
for f_region = {'r_M1'} %fieldnames(region_data)'
    f_region = f_region{1};

    for f_mouse = fieldnames(region_data.(f_region))'
        f_mouse = f_mouse{1};
        
        % Fields with leading m indicate a mouse entry
        if strcmp(f_mouse(1), 'm') ~= 1
            continue;
        end
    
        for f_rec = fieldnames(region_data.(f_region).(f_mouse))'
            f_rec = f_rec{1};
 
            % Check if each stimulation condition exists for the neuron
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_40')
                current = region_data.(f_region).(f_mouse).(f_rec).f_40.currents;
                current(find(current == 0)) = [];
                
                plot(i + rand(1, length(current)), current, '.r');
                hold on;
            end
            i = i + 1;
            
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_140')
                current = region_data.(f_region).(f_mouse).(f_rec).f_140.currents;
                current(find(current == 0)) = [];
                plot(i + rand(1, length(current)), current, '.b');
                hold on;
            end

            % Add the mouse and recording date for the plot
            tick_labels{end + 1} = [erase(f_mouse, 'm_')]; % ' ' erase(f_rec, 'r_')
            tick_spots(end + 1) = i;

            i = i + 2;
        end 

        i = i + 3;
    end
end
xticks(tick_spots);
xticklabels(tick_labels);
ylabel('uamp');

saveas(gcf, [savepath 'Current/' 'all_mice_currents.png']);
saveas(gcf, [savepath 'Current/' 'all_mice_currents.pdf']);

% Plot the Normalized Vm Change vs Current amplitude

% Create fields in structure if they do not exist
function [result] = create_struct(src, f_region, f_mouse, f_rec, f_neuron, f_stim)
    result = src;
    % Create region field
    if isfield(result, f_region) == 0
        result.(f_region) = struct();
    end
    
    % Create mouse field
    if isfield(result.(f_region), f_mouse) == 0
        result.(f_region).(f_mouse) = struct();
    end

    % Create recording field
    if isfield(result.(f_region).(f_mouse), f_rec) == 0
        result.(f_region).(f_mouse).(f_rec) = struct();
    end

    % Create stim field
    if isfield(result.(f_region).(f_mouse).(f_rec), f_stim) == 0
        result.(f_region).(f_mouse).(f_rec).(f_stim) = struct();
        result.(f_region).(f_mouse).(f_rec).(f_stim).currents = [];
    end

    % Create neuron field
    % Old code for saving everything into the neuron field's
    %if isfield(result.(f_region).(f_mouse).(f_rec).(f_stim), f_neuron) == 0
    %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron) = struct();
    %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron).currents = [];
    %    
    %end
end
