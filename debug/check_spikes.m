clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
 
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

%% Loop through and plot the spike amplitude for each condition
for f_region = fieldnames(region_data)' % {'r_V1'} %
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    for f_stim=stims'
        f_stim = f_stim{1};
        
        popul_data = data_bystim.(f_stim);

        % Grab only modulated neurons
        nr_idxs = find(sum(popul_data.mod_matrix, 2) > 0);
           
        figure;
    
        % Loop through each neuron
        for nr=nr_idxs'
            num_s = length(popul_data.all_trial_spike_amp{nr});
            plot(repmat(nr, 1, num_s), popul_data.all_trial_spike_amp{nr}, '.');
            hold on;
        end


        % Title of condition
        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
    end
end
