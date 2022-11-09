clear all;
% Linux server
server_root_path = '/home/pierfier/handata_server/';

% Windows server
%server_root_path = 'Z:\';

f = filesep;

% Specify the datapath, at the momement all of the data is in Eric's
pv_data_path = [server_root_path 'EricLowet' f 'DBS' f 'PV' f '40' f];

% Grab all 40Hz data
ses=dir('*.mat');

% Store the FOVs that were not being used
not_used = {};

% Specify used imaging parameters
Fs = 828;

% Loop through all trial data
for idx=1:length(ses)
    % Load the data
    data = load([server_root_path ses(ind).name]);
    
    % Use some filtering condition
    if length(find(result.trial_vec==trial_numb(1) )) < 3500 & length(trial_numb) > 2
        
    else
        % Store the files that do not meet the trial length conditions
        not_used{end+1} = ses(ind).name;
    end
end

% Print out all of the unused matfiles
disp('Not used trials');
for i=1:length(not_used)
    disp(not_used{i});
end
