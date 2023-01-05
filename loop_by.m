
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

% Store Vm data by each individual stimulation condition
Vm_data_bystim = struct();

% Loop through each stimulation condition
for field = fieldnames(matfile_stim)'
    field = field{1};

    % Initialize field subthreshold array
    Vm_data_bystim.(field) = struct();
    Vm_data_bystim.(field).neuron_Vm = [];
    matfiles = matfile_stim.(field).names;    
    

    % Loop through each matfile of the current stimulation condition
    for matfile =matfiles
        % Read in the mat file of the current condition
        data = load([pv_data_path matfile{1}]);
        
        cur_fov_subVm = [];
        % Loop through each trial
        for trial_data=data.align.trial
            trial_data = trial_data{1};
            
            % Check if there is only one neuron
            if size(trial_data.detrend_traces, 2) == 1
                cur_fov_subVm = horzcat_pad(cur_fov_subVm, trial_data.spike_info.trace_ws(:));
            else
                cur_fov_subVm = horzcat_pad(cur_fov_subVm, trial_data.spike_info.trace_ws);
            end
        end
        
        %DEBUG
        figure;
        plot(nanmean(cur_fov_subVm, 2));

        % Average for each neuron and save the subthreshold Vm
        temp = Vm_data_bystim.(field).neuron_Vm;
        Vm_data_bystim.(field).neuron_Vm = horzcat_pad(temp, nanmean(cur_fov_subVm, 2));
    end % End looping through FOVs of a condition
end 

% Plot all of the subthreshold Vm for each stimulation condition
stims = fieldnames(Vm_data_bystim);
figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
for stim=stims'
    plot(nanmean(Vm_data_bystim.(stim{1}).neuron_Vm, 2));
    hold on;
end
legend(stims, 'Interpreter', 'none');
title('Average Sub Vm by Stimulation condition', 'Interpreter', 'none');

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
