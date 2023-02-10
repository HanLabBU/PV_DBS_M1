
clear all;
close all;

f = filesep;

%%%---- USER Modification--------
% Linux server
server_root_path = '~/Projects/';

% Windows server
%server_root_path = 'Z:\';


% Data on server
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

% Specify which ignore trials to use
ignore_trial_dict = Multi_func.csv_to_struct([server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Recordings' f 'Data_Config' f 'byvis_ignore.csv']);
%%%----- END Modification---------

% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);
all_matfiles = {ses.name};

% Loop through each brain region
[region_matfiles] = find_region(all_matfiles);
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    % Create brain region struct
    region_stats.(f_region) = struct();
    
    %% Select matfiles by stim condition by this iteration's brain region
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);
    
    %% Loop through each field of the struct and concatenate everything together
    % Keep track the stats of the number of trials and neurons per mouse
    cond_stats = struct();
    % Loop through each stimulation condition
    for field = fieldnames(matfile_stim)'
        field = field{1};
    
        matfiles = matfile_stim.(field).names;    
     
        % Initialize field to keep track of neurons and trials
        cond_stats.(field) = struct();
        cond_stats.(field).trial_nums = [];
        cond_stats.(field).num_fovs = 0;
    
    
        % Loop through each matfile of the current stimulation condition
        for matfile =matfiles
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile{1}]);
     
            %TODO need to loop through ROIs
    
            %% OLD CODE
            % Grab the number of trials to ignore from this neuron
            %ri = strsplit(matfile{1}, '_');
            %if isfield(ignore_trial_dict.(['mouse_' ri{1}]), (['rec_' erase(ri{3}, 'rec')])) && ...
            %    isfield(ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]), (ri{4})) && ...
            %    isfield(ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}), (['f_' ri{5}]))
    
            %    num_ignore = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).ROI1;
            %else
            %    num_ignore = 0;
            %end
    
            %TODO I should just encapsulate this into a try/catch block and then catch with num_ignore = 0
            ri = strsplit(matfile{1}, '_');
            try 
                % Try to grab ROI1's ignoring, if doesn't exist look at ROI2 
                try
                    num_ignore = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).ROI1;
                catch
                    num_ignore = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).ROI2;
                end
                
                num_ignore = length(num_ignore);
            catch
                num_ignore = 0;
            end
    
            % Add the number of trials for this specific neuron and check if at least one trial is included in calculation
            cond_stats.(field).trial_nums(end + 1) = length(data.align.trial) - num_ignore;
            
            if length(data.align.trial) - num_ignore > 2
                cond_stats.(field).num_fovs = cond_stats.(field).num_fovs + 1;
            end
    
        end % End looping through FOVs of a condition
    end

    % Save the stim conditions for this iteration's brain region
    region_stats.(f_region).cond_stats = cond_stats;
end

% Loop through and show each regions conditions
regions = fieldnames(region_stats);
for f_region = regions'
    f_region = f_region{1};
    disp(f_region);

    stims = fieldnames(region_stats.(f_region).cond_stats);
    cond_stats = region_stats.(f_region).cond_stats;
    for stim=stims'
        stim = stim{1};
        disp(stim);
        disp(['Average trial per neuron ' num2str(nanmean(cond_stats.(stim).trial_nums)) ]);
        disp(['Number of neurons ' num2str(cond_stats.(stim).num_fovs) ]);
    end

    fprintf('\n\n');
end

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


%% Specific functions for determining which FOVs to look at
function [region_struct] = find_region(matfile_names)
    region_struct = struct();

    for i=1:length(matfile_names)
        file_parts = split(matfile_names{i}, '_');
        reg = file_parts{2};
        
        % Create region field if it does not exist
        if ~isfield(region_struct, ['r_' reg])
            region_struct.(['r_' reg]).names = {};
        end

        region_struct.(['r_' reg]).names{end+1} = matfile_names{i};
    end
end
