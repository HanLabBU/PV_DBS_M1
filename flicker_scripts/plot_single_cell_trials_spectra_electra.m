clear all;
close all;
clc;

f = filesep;
local_rootpath = '~/Projects/';
server_rootpath = '/home/pierfier/handata_server/eng_research_handata3/';


% Use fixed framerate
Fs = 1/(1.5 * 10^-3);

% Modifying which datapath to use
data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f];
%data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f '109567_Vb_male' f];

%savepath = [Multi_func.save_plot 'Quick_Analysis' f 'Random_Blocks' f];
%savepath = [Multi_func.save_plot 'Quick_Analysis' f 'PhaseWidthMapping' f];

% if 1: Have user select folder to use for the quick check
% if 0: Recursively search root data_path to find area folders
select_folder = 1;

% TODO change to flicker experiment matfiles
if select_folder == 1
    [fname,fdir] = uigetfile([data_path '*.mat'],'MultiSelect','on');
else
    %Looking for diretories that contain an 'area' in the subfolder item
    list = dir([data_path '**/*flicker*.mat']);
    list_dirs = unique({list.folder});
end

data = load([fdir f fname]);

traces = data.roi_list;

%%
% Construct filterbank
freqLimits = [0 200];

time_t = [1:length(traces.traces(traces.trial_vec == 1))]./Fs;

% Plot trails with spikes detected
figure;
all_f = [];
all_wt = [];
for tr_n=unique(traces.trial_vec)
    sp_info = data.spike_detect_SA_v4_info{tr_n};
    sp_idx = sp_info.spike_idx{1};
    raw_trace = traces.traces(traces.trial_vec == tr_n, 1);
    trace_mov = movmean(raw_trace, Fs);
    
    detrend_trace = -1*(raw_trace - trace_mov);
    max_t = max(detrend_trace);
    min_t = min(detrend_trace);
    
    norm_trace = (detrend_trace - min_t)./(max_t - min_t);
    
    plot(time_t, (1.5*tr_n) + norm_trace, 'k');
    hold on;
    %TODO plot the spikes as well
    plot(time_t(sp_idx), (1.5*tr_n) + norm_trace(sp_idx), '.r');
    hold on;
    
    fb = cwtfilterbank(SignalLength=length(detrend_trace),...
                   SamplingFrequency=Fs,...
                   FrequencyLimits=freqLimits);

    [wt, f] = cwt(detrend_trace, FilterBank=fb);

    % Append to array
    all_f = [all_f, f];

    all_wt = cat(3, all_wt, wt);;
end

% Perform trial-averaged spectra
avg_f = mean(all_f, 2, 'omitnan');
avg_wt = mean(abs(all_wt), 3, 'omitnan');

% Plot the spectra
figure;
surface(time_t, avg_f, avg_wt, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
