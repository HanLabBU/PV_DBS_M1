% Possible sub Vm's
% '50373_M1_rec20230718_FOV1_140_100' trial 6 or 7?
% '50373_M1_rec20230801_FOV2_40_200' trial 6 or 8? or 13
% '617100_M1_rec20211110_FOV3_40_60' trial 4 <---- TODO  7
clc;
clear all;
close all;
f = filesep;

% Option for excluding or including first 200ms
exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;

% Maingear office computer
local_root_path = '/home/pierfier/Projects/';
server_root_path = '~/handata_server/';
% Local windows machine
%local_root_path = 'Z:\Local_Data\';
%server_root_path = 'X:\';
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
pv_data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Path to save the figures
%savefig_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Figures' f 'Exemplary' f];
% handata3 server path
%savefig_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'Plots' f];
savefig_path = Multi_func.save_plot();

% Check if the figure path exists
if ~exist(savefig_path)
    disp('Figure path not found');
    return;
end

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end

tic;
%Load the data
load(save_all_data_file);
toc

% Get only M1 data
data_bystim = region_data.r_M1;

avg_Fs = mean(region_data.r_M1.f_40.framerate, 'omitnan');

%% Select neuron and trial
neuron = '617100_M1_rec20211110_FOV3_40_60';
trial_num = 4;
neuron_idx = find(contains(data_bystim.f_40.neuron_name, neuron));

%% Plot the Vm with overlay
figure('Renderer', 'painters', 'Position', [0 0 1000 1000]);
ax = gca;
ax.Units = 'centimeters';
ax.InnerPosition = [2 2 15.37 4.0];

trace_noise = data_bystim.f_40.all_trial_trace_noise{neuron_idx}(trial_num);
rawvm = data_bystim.f_40.all_trial_rawVm{neuron_idx}(:, trial_num)./trace_noise;
subvm = data_bystim.f_40.all_trial_SubVm{neuron_idx}(:, trial_num)./trace_noise;
plot(rawvm, 'k');
hold on;

% Remove the unfiltered subVm plot
%plot(subvm, 'b');
%hold on;
spikes = data_bystim.f_40.all_trial_spikeidx{neuron_idx}(:, trial_num);
plot(spikes, 14*ones(length(spikes)), '.g');
hold on;
filt_vm = filt_range(subvm, [2 10], avg_Fs);
plot(filt_vm, 'r');
hold on;
time_idx = [0 round(avg_Fs)/4];
plot(time_idx, [-10, -10], 'k', 'LineWidth', 4);
hold on;
text(time_idx(1), -12, [num2str(range(time_idx)/round(avg_Fs)) 'S']);
hold on;
snr_scale = 5;
plot([-1, -1], [0 snr_scale], 'k', 'LineWidth', 4);
hold on;
ht = text(- 5, 0, ['SBR ' num2str(snr_scale)]);
set(ht, 'rotation', 90);
Multi_func.set_default_axis(gca);
axis off;

saveas(gcf, [savefig_path 'Exemplary' f neuron '_tr_' num2str(trial_num) '_subthresh_overlay.pdf']);

%%

% Function to filter out range
function [filt_sig] = filt_range(sig, range, FS)
    Fn = FS/2;
    FB = [0.8 1.2].*range;
    
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    filt_sig = filtfilt(B,A,sig);
end

function [filt_sig] = hilb_filt_range(sig, range, FS)
    Fn = FS/2;
    FB = [0.8 1.2].*range;
    
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    filt_sig = hilbert(filtfilt(B,A,sig));
end
