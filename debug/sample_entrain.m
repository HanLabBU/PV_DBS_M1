close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% Parameters for frames to chop off
front_frame_drop = 15;
back_frame_drop = 2496;


figure_path = Multi_func.save_plot;

% Read in the saved pv data and perform analysis
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
%Load the data
load(save_all_data_file);

% Only M1 data
f_region = 'r_M1';
data_bystim = region_data.(f_region);
stims = fieldnames(data_bystim);

%% Grab Vm and Low Frequency Entrainment
% Vm and low frequency power
for f_stim=stims'
    f_stim = f_stim{1};
    cur_vm_change = [];
    cur_lf_change = [];
    for i = 1:size(data_bystim.(f_stim).trace_timestamps, 2)
        base_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) < data_bystim.(f_stim).stim_timestamps(1, i));
        stim_idx = find(data_bystim.(f_stim).trace_timestamps(:, i) >= data_bystim.(f_stim).stim_timestamps(1, i) & ...
                    data_bystim.(f_stim).trace_timestamps(:, i) <= data_bystim.(f_stim).stim_timestamps(end, i));
        Vm_base = nanmean(data_bystim.(f_stim).neuron_Vm(base_idx, i));
        Vm_stim = nanmean(data_bystim.(f_stim).neuron_Vm(stim_idx, i));
        cur_vm_change = [cur_vm_change, (Vm_stim - Vm_base)./(Vm_base + Vm_stim)];
               
        LF_pow_base = nanmean(abs(data_bystim.(f_stim).neuron_hilbfilt{i}(2:10, base_idx, :)), 'all');
        LF_pow_stim = nanmean(abs(data_bystim.(f_stim).neuron_hilbfilt{i}(2:10, stim_idx, :)), 'all');
        cur_lf_change = [cur_lf_change, (LF_pow_stim - LF_pow_base)./(LF_pow_base + LF_pow_stim)];
    end
    data_bystim.(f_stim).Vm_change = cur_vm_change;
    data_bystim.(f_stim).Low_Freq_change = cur_lf_change;
end

%% Plot the data
figure;
i = 0;
for f_stim=stims'
    f_stim = f_stim{1};
    col = [0 i 1];
    cur_Vm = data_bystim.(f_stim).Vm_change;
    cur_LF = data_bystim.(f_stim).Low_Freq_change;
    plot(cur_Vm, cur_LF, '.', 'color', col, 'MarkerSize', 10, 'DisplayName', [f_stim]);
    hold on;
    i = i + 1;
    b = polyfit(cur_Vm, cur_LF, 1);
    xfit = min(cur_Vm):0.1:max(cur_Vm);
    fit = polyval(b, xfit);
    plot(xfit, fit, 'color', col);
    hold on;
end
xlabel('Vm');
ylabel('2-10Hz Power');
legend('Interpreter', 'none');
