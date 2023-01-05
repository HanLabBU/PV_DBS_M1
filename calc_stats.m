
clear all;
close all;

f = filesep;

% USER Modification
% Linux server
server_root_path = '/home/pierfier/handata_server/';

% Windows server
%server_root_path = 'Z:\';
% END Modification

% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

pv_data_path = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

ses = dir([pv_data_path '*.mat']);

all_matfiles = {ses.name};

%% Select matfiles by specific conditions

[matfile_stim] = stim_cond(all_matfiles);

%% Loop through each field of the struct and concatenate everything together

% TODO maybe try to calculate the number of neurons and trials per
%mouse

% Keep track the stats of the number of trials and neurons per mouse
cond_stats = struct();

% Loop through each stimulation condition
for field = fieldnames(matfile_stim)'
    field = field{1};

    matfiles = matfile_stim.(field).names;    
 
    % Initialize field to keep track of neurons and trials
    cond_stats.(field) = struct();
    cond_stats.(field).trial_nums = [];
    cond_stats.(field).num_neurons = 0;


    % Loop through each matfile of the current stimulation condition
    for matfile =matfiles
        % Read in the mat file of the current condition
        data = load([pv_data_path matfile{1}]);
 
        % Add the number of trials for this specific neuron 
        cond_stats.(field).trial_nums(end + 1) = length(data.align.trial);

        % Increment number of neurons for condition
        cond_stats.(field).num_neurons = cond_stats.(field).num_neurons + 1;

    end % End looping through FOVs of a condition
end 

% Plot all of the subthreshold Vm for each stimulation condition
stims = fieldnames(cond_stats);
for stim=stims'
    field = stim{1};
    disp(field);
    disp(['Average trial per neuron for condition ' num2str(nanmean(cond_stats.(field).trial_nums)) ]);
    disp(['All trials of condition ' num2str(nansum(cond_stats.(field).trial_nums)) ]);
    
    disp(['Number of neurons ' num2str(cond_stats.(field).num_neurons)]);
    fprintf('\n');

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
