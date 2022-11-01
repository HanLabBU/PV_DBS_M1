%% Script to loop through a dataset that is aligned and saved into Matlab
% Possibly only need to be ran on linux machine

f = filesep;

%USER modify everything below here
%% Specific to the computer being used
%DMD
server_rootpath = 'Z:\';

%Dual scope
%server_rootpath = '\\engnas.bu.edu\research\eng_research_handata\';

% Maingear office computer
%server_rootpath = '/home/pierfier/handata_server/';

% Folder location of all the saved and aligned data
data_path = [server_rootpath 'Jenn' f 'SomArchon' f 'A2A_Cre' f 'D2 Data auto' f];

%USER END modification

% Scripts folder for performing data alignment and saving relevant data to mat files
addpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'LFP_Trace_Alignment' ]);

matfile_names = dir(data_path);
isfile_idx = ~[matfile_names.isdir];
matfile_names = {matfile_names.name};
matfile_names = matfile_names(isfile_idx);

% Loop through each FOV matfile
for i=1:length(matfile_names)
    
    % Only open matfiles
    if ~contains(matfile_names{i}, '.mat')
        continue
    end
    
    % Check which trials to ignore
    fov_dict = [];
    fov_dict.a = [];
    [cur_mouse, cur_rec, cur_fov] = get_session_info(matfile_names{matfile_idx});
    if isfield(ignore_trial_dict, ['mouse_' cur_mouse]) && ...
        isfield(ignore_trial_dict.(['mouse_' cur_mouse]), ['rec_' cur_rec]) && ...
        isfield(ignore_trial_dict.(['mouse_' cur_mouse]).(['rec_' cur_rec]), ['FOV' cur_fov])
        
        fov_dict = ignore_trial_dict.(['mouse_' cur_mouse]).(['rec_' cur_rec]).(['FOV' cur_fov]);
        
        disp(['DEBUG ignoring ' cur_mouse ' ' cur_rec ' ' cur_fov]);
        
        if strcmp(fov_dict, 'all')
            continue;
        end
    end

    % Load the FOV data
    data = load([data_path matfile_names{matfile_idx}]);

    % Loop through each trial
    for j=1:length(data.align.trial)
        trial_data = data.align.trial{j};
            
        % Check if it is an empty trial and move on 
        if isempty(trial_data) || strcmp(fov_dict.([rois{l}]), 'all')
            continue
        end

        % Ignore specific ROI traces
        rois = fieldnames(fov_dict);
        %trace_idx = zeros(size(data.align.trial{j}.detrend_traces));
        for l=1:length(rois)
            if  ismember(j, fov_dict.([rois{l}]))
                %DEBUG
                disp(['NaNs for ' rois{l} ' of trial ' num2str(j)]);
                trial_data.detrend_traces(:, str2num(erase('ROI', rois{l}))) = NaN;
            end
        end
        
        % TODO here is where individual trial functions need to be implemented
        

    end
end

% 

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
