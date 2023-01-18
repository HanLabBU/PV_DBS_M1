clear all;
close all;
clc;
%% Script to loop through all D2 data that is aligned and saved into Matlab
% Possibly only need to be ran on linux machine
f = filesep;

%USER modify everything below here
%DMD
%server_rootpath = 'Z:\';
%Dual scope
%server_rootpath = '\\engnas.bu.edu\research\eng_research_handata\';
% Maingear office computer
server_rootpath = '/home/pierfier/handata_server/';

% Folder location of all the saved and aligned data
%data_path = [server_rootpath 'Jenn' f 'SomArchon' f 'A2A_Cre' f 'D2 Data auto' f];
% Local maingear linux folder of PV data
data_path = ['~/Projects' f 'PV DBS Project' f 'PV_Data' f];

%USER END modification

matfile_names = dir(data_path);
isfile_idx = ~[matfile_names.isdir];
matfile_names = {matfile_names.name};
matfile_names = matfile_names(isfile_idx);

ignored_trials_dict = [];
total_num_trials = 0;
ignored_num_trials = 0;
all_trial_snr = [];

snr_threshold_trace = 3.5;

% Loop through each FOV matfile
for matfile_idx=1:length(matfile_names)
    
    % Only open matfiles
    if ~contains(matfile_names{matfile_idx}, '.mat')
        continue
    end
    
    % Load the FOV data
    data = load([data_path matfile_names{matfile_idx}]);
    mat_parts = strsplit(matfile_names{matfile_idx}, '_');
    cur_mouse = mat_parts{1};
    cur_rec = erase(mat_parts{3}, 'rec');
    cur_fov = mat_parts{4};
    cur_cond = mat_parts{5};

    % Loop through each trial
    for j=1:length(data.align.trial)
        trial_data = data.align.trial{j};
        if isempty(trial_data)
            continue
        end

        % Loop through each ROI
        for k=1:size(trial_data.detrend_traces, 2)
            ave_snr = mean(trial_data.spike_info.spike_snr{k}); 
            all_trial_snr = [all_trial_snr ave_snr];
            total_num_trials = total_num_trials + 1;

            %DEBUG
            if ave_snr > 5
                matfile_names{matfile_idx}
                j
                k
            end

            % Ignore based on the SNR threshold
            if ave_snr < snr_threshold_trace
                ignored_num_trials = ignored_num_trials + 1;
                
                % Add this particular trace to the dictionary to ignore
                mouse_key = ['mouse_' cur_mouse];
                rec_key = ['rec_' cur_rec];
                FOV_key = [cur_fov];
                ROI_key = ['ROI' num2str(k)];
                cond_key = ['f_' cur_cond];

                % Append trial to list of ignored trials for this ROI, otherwise just initialize array
                try 
                    roi_trial_array = ignored_trials_dict.(mouse_key).(rec_key).(FOV_key).(cond_key).(ROI_key);
                    ignored_trials_dict.(mouse_key).(rec_key).(FOV_key).(cond_key).(ROI_key) = [roi_trial_array j];
                catch
                    ignored_trials_dict.(mouse_key).(rec_key).(FOV_key).(cond_key).(ROI_key) = [j];
                end
            end
        end

    end
end

%Plot the distribution of all the SNRs
figure;
boxplot(all_trial_snr);
hold on;
plot(repmat(1, size(all_trial_snr, 1), size(all_trial_snr, 2)), all_trial_snr, '.b');
title('All SNRs of every single trace');

% Display ignoring info
disp(['Ignored ' num2str(ignored_num_trials) ' trials of ' num2str(total_num_trials)]);

% Create CSV of all the ignored trials
struct_to_csv(ignored_trials_dict, ['thres_' num2str(snr_threshold_trace) '_ignore.csv']);

% Grab a string that contains the list of the filenames
function [selected_folder, selected_files] = str_to_list(paths_string)
    f = filesep;

    path = paths_string{1};
    path = erase(path, '[');
    path = erase(path, ']');
    path = erase(path, "'");
    path_list = split(path, ', ');
    
    % Find the directory of where all of the files are stored
    last_slash_index = find(path_list{1} == f, 1, 'last');
    first_path = path_list{1};
    selected_folder = first_path(1:last_slash_index);

    for i=1:length(path_list)
        last_slash_index = find(path_list{i} == f, 1, 'last');
        cur_path = path_list{i};
        path_list{i} = cur_path(last_slash_index+1:end);
    end

    selected_files = path_list;
end

% Convert struct to csv, makes a new column until individual ROIs are specified. At that point there is a string representation of a dictionary
function [error] = struct_to_csv(arg_dict, csv_filepathname)
    error = [];

    % Intialize string representation of the csv file with the headers (standard for our data)
    str_array = {'', 'mouse', 'rec_date', 'FOV', 'cond', 'roi_trial_dict', 'comments'};
    index = 0;

    % Loop through the mice
    mice = fieldnames(arg_dict);
    for i=1:length(mice)
        
        % Loop through the sessions
        sessions = fieldnames(arg_dict.(mice{i}));
        for j=1:length(sessions)
            
            % Loop through the FOVs
            fovs = fieldnames(arg_dict.(mice{i}).(sessions{j}));
            for k=1:length(fovs)
                
                % Loop through each stimulation condition
                conds = fieldnames(arg_dict.(mice{i}).(sessions{j}).(fovs{k}));
                for l=1:length(conds)
                    % Take the rest of the dictionary and convert it to a string representation
                    roi_dict = arg_dict.(mice{i}).(sessions{j}).(fovs{k}).(conds{l});
                    [roi_trials_str] = dict_to_string(roi_dict);
                    row = {num2str(index), erase(mice{i}, 'mouse_'), erase(sessions{j}, 'rec_'), erase(fovs{k}, 'FOV'), erase(conds{l}, 'f_'), roi_trials_str, ''};
                    str_array = cat(1, str_array, row);
                    index = index + 1;
                end
            end
        end
    end

    fileID = fopen(csv_filepathname, 'w');
    for i=1:size(str_array, 1)
        row = str_array(i, :);
        for j=1:length(row) - 1
            fprintf(fileID, [row{j} ';']);
        end
        fprintf(fileID, [row{end} '\n']);
    end
    fclose(fileID);
end

% Convert ROI with trials to ignore dictionary 
function [str_rep] = dict_to_string(arg_dict)
    str_rep = '{';
    
    % Loop through ROIs
    rois = fieldnames(arg_dict);
    for i=1:length(rois)
        str_rep = [ str_rep, erase(rois{i}, 'ROI'), ':['];

        % Loop through individual trials
        trials = arg_dict.(rois{i});
        for j=1:length(trials)
            str_rep = [str_rep, num2str(trials(j)), ','];
        end
        str_rep = [str_rep(1:end-1), '],'];
    end
    str_rep = [str_rep(1:end-1), '}'];

end

% Read the data within the saved matfile filename
function [mouse_str, rec_date, FOV_str] = get_session_info(matfile_name)
    % Split elements of filename
    elem = split(matfile_name, '_');
    
    mouse_str = elem{1};
    rec_date = erase(elem{2}, 'rec');
    FOV_str = erase(elem{3}, 'FOV');
    FOV_str = erase(FOV_str, '.mat');
end
