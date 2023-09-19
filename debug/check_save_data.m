clc;
clear all;
close all;

f = filesep;

local_root_path = '~/Projects/';

% Read in the interim data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];

load(save_all_data_file);

data = region_data.r_V1;

for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    trace_time = data.(f_stim).trace_timestamps*1000;
    stim_time = data.(f_stim).stim_timestamps*1000;

    figure;
    % Loop through each trial
    plot(trace_time, repmat([1:size(trace_time, 2)], size(trace_time, 1), 1), '|r');
    hold on;
    plot(stim_time, repmat([1:size(stim_time, 2)] + 0.5, size(stim_time, 1), 1), '|b');
    
    %TODO Check before and after chopping of frames
    
end
