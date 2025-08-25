%% Close and clear all
clc;
close all;
clear all;

%%
f= filesep;
server_rootpath = '/home/pierfier/handata_server/eng_research_handata3/';

% Use fixed framerate
Fs = 500;

% Modifying which datapath to use
data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f];
%data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f '109567_Vb_male' f];

% Pick matfile
[file, path] = uigetfile('*.mat', 'Select a file', data_path);

data = load([path file]);

%% Loop through and calculate the average spectra for each trial
base_idx = 1:Fs;
stim_idx = Fs:2*Fs;

all_trial_spec = [];
all_trial_f = [];
freqLimits = [0 200];
for i=unique(data.roi_list.trial_vec)
    cur_trace = data.roi_list.traces(find(data.roi_list.trial_vec == i));

    % Create filterbank
    fb = cwtfilterbank(SignalLength=length(cur_trace),...
                       SamplingFrequency=Fs,...
                       FrequencyLimits=freqLimits);
    
    [wt, f] = cwt(cur_trace, FilterBank = fb);
    cur_pow = abs(wt);

    % Normalize the power
    base_pow = nanmean(cur_pow(:, base_idx), 2);
    stim_pow = nanmean(cur_pow(:, stim_idx), 2);
    cur_pow = (cur_pow - base_pow)./(base_pow + stim_pow);

    % Add trial to massive array
    all_trial_spec = cat(3, all_trial_spec, cur_pow);
    all_trial_f = cat(3, all_trial_f, f);

    % Plot the individual spectra trials
    figure;
    surface(1:size(wt, 2), ...
            f, ...
            cur_pow, 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        colormap(jet*.8);
        colorbar;
    title(['Inidividual trial ' num2str(i)]);
end

% Average across all trials
figure;
surface(1:size(all_trial_spec, 2), nanmean(all_trial_f, 3), ...
        nanmean(all_trial_spec, 3), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
colormap(jet*.8);
colorbar;
title('All Trial Average');
