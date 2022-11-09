clear all;
% Linux server
%server_root_path = '/home/pierfier/handata_server/';

% Windows server
server_root_path = 'Z:\';

f = filesep;

% Specify the datapath, at the momement all of the data is in Eric's
pv_data_path = [server_root_path 'EricLowet' f 'DBS' f 'PV' f ];

savepath = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f];

% Specify analysis parameters
Fs = 828;
front_frame_drop = 15;
back_frame_drop = 2500;
stim_conditions = {'40', '140'};

% Store all of the traces, use the string as a key for the data
fields = [strcat('stim', stim_conditions); {[], []}];
trial_avg_norm_stim_traces = struct(fields{:});

% Store the FOVs that were not being used
not_used = {};

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
                
                [x, y] = exp_fit(cur_stim_traces(:, end), Fs);
                cur_stim_traces(:, end) = cur_stim_traces(:, end) - y';
                
                %DEBUG for the exponential fit
                %figure;
                %tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
                %nexttile;
                %plot(trace);
                %nexttile;
                %plot(cur_stim_traces(:, end));

                % Only use spikes within truncated region
                spx = data.result.resultS{trial}.spike_idx{1};
                samp = data.result.resultS{trial}.spike_amplitude{1};
                samp(spx <= front_frame_drop | spx >= back_frame_drop) = NaN;
                
                % Save all of the spike amplitudes in one matrix
                cur_spike_amp = [cur_spike_amp; samp(:)];
            end
            
            % Saving trial averaged trace normalized by average spike amplitude for all FOVs
            trial_avg_norm_stim_traces.(['stim' cond{1}])(:, end+1) = nanmean(cur_stim_traces, 2)./nanmean(cur_spike_amp);
        else
            % Store the files that do not meet the trial length conditions
            not_used{end+1} = ses(idx).name;
        end
    end
end

% Calculate the time axis
tim_axis = (front_frame_drop:back_frame_drop)./Fs;
stimoffset_tim_axis = tim_axis - 1; % Offset in seconds

% Print out median of voltage trace
figure('Renderer', 'painters');
stim40_trace_median = nanmedian(trial_avg_norm_stim_traces.stim40, 2);
stim140_trace_median = nanmedian(trial_avg_norm_stim_traces.stim140, 2);
plot(stimoffset_tim_axis, stim40_trace_median,'b','Linewidth',1.5);
hold on;
plot(stimoffset_tim_axis, stim140_trace_median,'r','Linewidth',1.5)
axis tight;
xlim([-0.7 1.7]);
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Nonsmooth_trace_median_40_140t.pdf']);


% Print out all of the unused matfiles
disp('Not used trials');
for i=1:length(not_used)
    disp(not_used{i});
end

% Function that will do exponential fit and return the x and y values for the fit line
% This uses the Curve Fitting Toolbox to get the values
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
