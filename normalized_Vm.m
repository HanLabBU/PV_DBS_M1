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
front_frame_drop = 15;
back_frame_drop = 2500;
stim_conditions = {'40', '140'};

% Store all of the traces, use the string as a key for the data
fields = [strcat('stim', stim_conditions); {[], []}];
trial_avg_norm_stim_traces = struct(fields{:});


for cond=stim_conditions
    
    %Grab the respective frequency stimulation data
    ses=dir([pv_data_path cond{1} f '*.mat']);
    
    % Loop through all 40Hz FOV data
    for idx=1:length(ses)
    
        % Load the data
        data = load([pv_data_path ses(idx).name]);
        
        % Use some filtering condition
        if length(find(data.result.trial_vec==1 )) < 3500 & length(unique(data.result.trial_vec)) > 2 & isfield(data.result, 'resultS')
           
            cur_stim_traces = [];
            cur_spike_amp = [];
            % Loop through all trials
            for trial=unique(data.result.trial_vec)
                
                % Grab the full voltage trace
                trace = data.result.traces(trial == data.result.trial_vec);
                %TODO got to detrend with the exponential fit
                cur_stim_traces(:, end+1) = trace(front_frame_drop:back_frame_drop);
                
                % Only use spikes within truncated region
                spx = data.result.resultS{trial}.spike_idx{1};
                samp = data.result.resultS{trial}.spike_amplitude{1};
                samp(spx <= front_frame_drop | spx >= back_frame_drop) = NaN;
                
                % Save all of the spike amplitudes in one matrix
                cur_spike_amp = [cur_spike_amp; samp(:)];
            end
            
            % Saving trial averaged trace normalized by average spike amplitude for all FOVs
            trial_avg_norm_stim_traces.(['stim' cond{1}])(:, end+1) = nanmean(cur_stim_traces, 2)./nanmean(cur_spike_amp)
        else
            % Store the files that do not meet the trial length conditions
            not_used{end+1} = ses(idx).name;
        end
    end
end

% Print out median of voltage trace


% Print out all of the unused matfiles
disp('Not used trials');
for i=1:length(not_used)
    disp(not_used{i});
end
