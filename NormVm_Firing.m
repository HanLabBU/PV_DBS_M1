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


% Specify the datapath, at the momement all of the data is in Eric's
pv_data_path = [server_root_path 'EricLowet' f 'DBS' f 'PV' f ];

savepath = [server_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f];

%% Specify analysis parameters
Fs = 828;
front_frame_drop = 15;
back_frame_drop = 2500;
stim_conditions = {'40', '140'};

%% Spike rate stuff
% Window for spike rate calculation
srate_window = 26;

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
            fov_stim_traces = [];
            fov_spike_amp = [];
            fov_cont_spike_rate = [];
            fov_sraster = [];

            % Loop through all trials
            for trial=unique(data.result.trial_vec)
                
                % Grab and detrend the full voltage trace
                trace = data.result.traces(trial == data.result.trial_vec);
                fov_stim_traces(:, end+1) = trace(front_frame_drop:back_frame_drop);
                [x, y] = exp_fit(fov_stim_traces(:, end), Fs);
                fov_stim_traces(:, end) = fov_stim_traces(:, end) - y';

                % Grab spike amplitudes within truncated region
                % Not sure if the indexing is shifted as previous spots that use front_frame_drop:back_frame_drop
                spike_idx = data.result.resultS{trial}.spike_idx{1};
                samp = data.result.resultS{trial}.spike_amplitude{1};
                samp(spike_idx <= front_frame_drop | spike_idx >= back_frame_drop) = NaN;
                
                % Save all of the spike amplitudes in one matrix
                fov_spike_amp = [fov_spike_amp; samp(:)];

                %% Spike Data
                % Grab spike raster
                spike_raster = data.result.resultS{trial}.roaster(front_frame_drop:back_frame_drop);
                % Calculate the spike rate over time
                [t cont_spike_rate] = calc_cont_spike_rate(spike_raster, srate_window, Fs);
                
                fov_cont_spike_rate(:, end+1) = cont_spike_rate;
                fov_sraster(:, end+1) = spike_raster;
            end
            
            %%Vm data
            % Saving trial averaged trace normalized by average spike amplitude for all FOVs
            trial_avg_norm_stim_traces.(['stim' cond{1}])(:, end+1) = nanmean(fov_stim_traces, 2)./nanmean(fov_spike_amp);
            
            %%Spike data
            
            %Plotting individual spike rate with rasters aligned
            figure('Renderer', 'Painters', 'Position', [200 200 1000 1000]);
            title(ses(idx).name);
            tiledlayout(2, 1, 'TileSpacing','compact', 'Padding','compact');
            nexttile;
            for i=1:size(fov_sraster, 2)
                spike_idx = find(fov_sraster(:, i) == 1);
                plot(spike_idx, repmat(i, length(spike_idx), 1), 'b.', 'MarkerSize', 6);
                hold on;
            end
            nexttile;
            plot(fov_cont_spike_rate);
            sgtitle(ses(idx).name, 'Interpreter', 'none');
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
% Calculates Firing rate using a gaussian shifting window with specified gaussian sigma (in ms)
function [t cont_spike_rate] = calc_cont_spike_rate(raster, gauss_sigma, Fs)
    cont_spike_rate = [];

    % Formula taken from gausswin matlab documentation
    L = 6*gauss_sigma +1;
    alpha = (L - 1)./(2*gauss_sigma);
    gw = gausswin(L, alpha);
    % Normalize the gaussian window
    gw = gw ./ sum(gw);
    
    % Apply gaussian filter across data
    cont_spike_rate = filtfilt(gw, 1, raster);
    t = (1:length(cont_spike_rate))./Fs;
end

% Function that will do exponential fit and return the x and y values for the fit line
% This uses the Curve Fitting Toolbox to get the values
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
