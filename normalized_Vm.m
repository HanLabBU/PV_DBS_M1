clear all;
% Linux server
server_root_path = '/home/pierfier/handata_server/';

% Windows server
%server_root_path = 'Z:\';

f = filesep;

% Specify the datapath, at the momement all of the data is in Eric's
pv_data_path = [server_root_path 'EricLowet' f 'DBS' f 'PV' f ];


% Store the FOVs that were not being used
not_used = {};

% Specify analysis parameters
Fs = 828;
init_frame_drop = 15;


% Store all of the traces, use the string as a key for the data
stim_traces = struct();

for cond={'40', '140'}
    
    %Grab the respective frequency stimulation data
    ses=dir([pv_data_path cond{1} f '*.mat']);
    
    % Loop through all 40Hz FOV data
    for idx=1:length(ses)
    
        % Load the data
        data = load([pv_data_path ses(idx).name]);
        
        % Use some filtering condition
        if length(find(data.result.trial_vec==1 )) < 3500 & length(unique(data.result.trial_vec)) > 2
            
            % Grab the full voltage trace
            for trial=unique(data.result.trial_vec)
                trace = data.result.traces(trial == data.result.trial_vec);
                stim_traces.(['stim' cond{1}])(:, trial) = trace;
            end
    
        else
            % Store the files that do not meet the trial length conditions
            not_used{end+1} = ses(ind).name;
        end
    end
end

% Print out all of the unused matfiles
disp('Not used trials');
for i=1:length(not_used)
    disp(not_used{i});
end
