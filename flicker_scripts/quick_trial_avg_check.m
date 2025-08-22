clear all;
close all;
clc;

f = filesep;
local_rootpath = '~/Projects/';
server_rootpath = '/home/pierfier/handata_server/eng_research_handata3/';


% Use fixed framerate
Fs = 500;

% Modifying which datapath to use
data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f];
%data_path = [server_rootpath 'Yangyang_Wang' f 'PV_V1_LED_SomArchon' f '109567_Vb_male' f];

%savepath = [Multi_func.save_plot 'Quick_Analysis' f 'Random_Blocks' f];
%savepath = [Multi_func.save_plot 'Quick_Analysis' f 'PhaseWidthMapping' f];

% if 1: Have user select folder to use for the quick check
% if 0: Recursively search root data_path to find area folders
select_folder = 0;

% TODO change to flicker experiment matfiles
if select_folder == 1
    [path] = uigetdir(data_path, 'Select the root folder that holds trial folders');
else
    %Looking for diretories that contain an 'area' in the subfolder item
    list = dir([data_path '**/*flicker*.mat']);
    list_dirs = unique({list.folder});
end

% Pulse times
pulse = [ (1:10), (1.1:10.1) ] - 1;

%%
for path_i = 1:length(list_dirs)
    path = list_dirs{path_i};
%%
    % Find matfiles in current path
    matfiles = dir([path f '*flicker*.mat']);
    matfiles = {matfiles.name};

    for matfile_i = 1:length(matfiles)
        matfile = matfiles{matfile_i};
    
        % Load the current areas traces data
        data = load([path f matfile]);
        data = data.roi_list;    

        % Store the spike detection info into the matfile 
        % The structure is {neuron, traces}
        spike_detect_SA_v4_info = {};

        % Loop through each roi
        for roi = 1:size(data.traces, 2)
            roi_vm = [];
            roi_raster = [];
            % Loop through each trial

            for tr = unique(data.trial_vec) %1:length(data.traces(roi).file)
                raw_trace = data.traces(data.trial_vec == tr, roi);
                trace_mov = movmean(raw_trace, Fs);
        
                % Flips the trace as well
                detrend_trace = (raw_trace - trace_mov)./trace_mov;
        
                roi_vm = horzcat_pad(roi_vm, detrend_trace);
                
                % Perform spike detection SomArchon version on detrended trace
                [spike_info] = spike_detect_SNR_v4(detrend_trace, Fs, 3); % Was originally 3.75
                spike_detect_SA_v4_info{roi, tr} = spike_info;
                
                roi_raster = horzcat_pad(roi_raster, spike_info.roaster');
            end            
        end

        % Save the spike detection back into matfiles
        save([path f matfile], 'spike_detect_SA_v4_info', '-append');
    end
end
    
    %%% Plot all of the roi traces in each area
    %for area_i = 1:length(area_data)
    %    if isempty(area_data{area_i})
    %        continue;
    %    end

    %    figure('Position', [100 100 500 200*length(area_data{area_i}.roi)]);
    %    tiledlayout(length(area_data{area_i}.roi), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    %    % Loop through each roi and plot the traces
    %    for roi_i = 1:length(area_data{area_i}.roi)
    %        roi_data = area_data{area_i}.roi{roi_i};
    %        nexttile;
    %        timeline = ( 4 + (1:size(roi_data.traces, 1) ))./Fs - 1;
    %        plot(timeline, mean(roi_data.traces, 2));

    %        if length(pulse) ~= 0
    %            hold on;
    %            xline(pulse);
    %        end
  
    %        %[mi, ma] = bounds(mean(roi_data.traces, 2));
    %        %ylim([0.9 1.1].*[mi ma]);
    %    end
    %    
    %    sgtitle([area_data{area_i}.name ' area' num2str(area_i) ] , 'Interpreter', 'none');
    %    savefig(gcf, [savepath area_data{area_i}.name ' area' num2str(area_i) ' _avg_vm.fig']);
    %    
    %end
    %
    %%% Plot the power spectra
    %avg_Fs = 500;
    %for area_i = 1:length(area_data)
    %    if isempty(area_data{area_i})
    %        continue;
    %    end
    %    figure('Position', [100 100 500 200*length(area_data{area_i}.roi)]);
    %    tiledlayout(length(area_data{area_i}.roi), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    %    % Loop through each roi and plot the traces
    %    for roi_i = 1:length(area_data{area_i}.roi)
    %        roi_data = area_data{area_i}.roi{roi_i};
    %
    %        all_pow = [];
    %        all_freq = [];
    %        % Loop through each trial
    %        for tr = 1:size(roi_data.traces, 2)
    %            % Check if trace values are Nans
    %            if sum(~isnan(roi_data.traces), 'all') ~= 0
    %                [wt, freq] = Multi_func.get_power_spec(roi_data.traces(:, tr), Fs);
    %                base_idx = 1:avg_Fs;
    %                stim_idx = avg_Fs:2*avg_Fs;
    %                base_avg = nanmean(abs(wt(base_idx)), 2);
    %                stim_avg = nanmean(abs(wt(stim_idx)), 2);
    %                wt = (abs(wt) - base_avg)/(base_avg + stim_avg);
    %                all_pow = cat(3, all_pow, abs(wt));
    %                all_freq = cat(3, all_freq, freq);
    %            end
    %        end
    %
    %        nexttile;
    %        Multi_func.set_default_axis(gca);
    %        timeline = ( 4 + (1:size(all_pow, 2) ))./Fs - 1;
    %        avg_spec = mean(all_pow, 3);
    %        surface(timeline, ...
    %                mean(all_freq, 3)', ...
    %                avg_spec, ...
    %                'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
    %            colormap(jet*.8);
    %
    %        %[mi, ma] = bounds(mean(roi_data.traces, 2));
    %        %ylim([0.9 1.1].*[mi ma]);
    %        
    %        if length(pulse) ~= 0
    %            hold on;
    %            xline(pulse);
    %        end
    %            
    %    end
    %    sgtitle(['A+B normalization' area_data{area_i}.name ' area' num2str(area_i)], 'Interpreter', 'none');
    %    savefig(gcf, [savepath area_data{area_i}.name ' area' num2str(area_i) ' _avg_spec.fig']);
    %end
    %
    %%% Plot each individual trials
    %for area_i = 1:length(area_data)
    %    
    %    if isempty(area_data{area_i})
    %        continue;
    %    end
    %    
    %    figure('Position', [100 100 500 200*length(area_data{area_i}.roi)]);
    %    tiledlayout(length(area_data{area_i}.roi), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    %    % Loop through each roi and plot the traces
    %    for roi_i = 1:length(area_data{area_i}.roi)
    %        roi_data = area_data{area_i}.roi{roi_i};
    %
    %        
    %        nexttile;
    %
    %        Multi_func.set_default_axis(gca);
    %        % Loop through each trial
    %        for tr = 1:size(roi_data.traces, 2)
    %            timeline = ( 4 + (1:size(roi_data.traces, 1) ))./Fs - 1;
    %            cur_trace = roi_data.traces(:, tr);
    %            min_trace = min(roi_data.traces(:, tr));
    %            max_trace = max(roi_data.traces(:, tr));
    %
    %            norm_trace = (cur_trace - min_trace)./(max_trace - min_trace);
    %            plot(timeline, norm_trace + (tr*1.5) + .5);
    %            hold on;
    %        end

    %        if length(pulse) ~= 0
    %            hold on;
    %            xline(pulse);
    %        end
    %
    %    end
    %    sgtitle([area_data{area_i}.name ' area' num2str(area_i) ], 'Interpreter', 'none');
    %    savefig(gcf, [savepath area_data{area_i}.name ' area' num2str(area_i) ' _traces.fig']);
    %end
    %
    %
    %%% Plot trials heatmap
    %for area_i = 1:length(area_data)
    %    
    %    if isempty(area_data{area_i})
    %        continue;
    %    end
    %    figure('Position', [100 100 500 200*length(area_data{area_i}.roi)]);
    %    tiledlayout(length(area_data{area_i}.roi), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    %    % Loop through each roi and plot the traces
    %    for roi_i = 1:length(area_data{area_i}.roi)
    %        roi_data = area_data{area_i}.roi{roi_i};
    %
    %        
    %        nexttile;
    %        Multi_func.set_default_axis(gca);
    %        timeline = ( 4 + (1:size(roi_data.traces, 1) ))./Fs - 1;
    %        imagesc('XData', timeline, 'YData', 1:size(roi_data.traces, 2), 'CData', roi_data.traces');
    %        
    %    end
    %    if length(pulse) ~= 0
    %        hold on;
    %        xline(pulse);
    %    end
    %    sgtitle([area_data{area_i}.name ' area' num2str(area_i) ], 'Interpreter', 'none');
    %    savefig(gcf, [savepath area_data{area_i}.name ' area' num2str(area_i) ' _heatmap.fig']);
    %end
    %
    %%% Plot the raster plots for each trial
    %for area_i = 1:length(area_data)
    %    
    %    if isempty(area_data{area_i})
    %        continue;
    %    end
    %    figure('Position', [100 100 500 200*length(area_data{area_i}.roi)]);
    %    tiledlayout(length(area_data{area_i}.roi), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    %    % Loop through each roi and plot the traces
    %    for roi_i = 1:length(area_data{area_i}.roi)
    %        roi_data = area_data{area_i}.roi{roi_i};
    %        nexttile;
    %        Multi_func.set_default_axis(gca);
    %        timeline = ( 4 + (1:size(roi_data.traces, 1) ))./Fs - 1;
    %        
    %        % loop through each trial
    %        for tr = 1:size(roi_data.rasters, 2)
    %            spike_idx = logical(roi_data.rasters(:, tr));
    %            plot(timeline(spike_idx), tr*roi_data.rasters(spike_idx, tr), '.');
    %            hold on;
    %            %ylim([0.5 size(roi_data.traces, 2) + 1]);
    %            %[mi, ma] = bounds(mean(roi_data.traces, 2));
    %            %ylim([0.9 1.1].*[mi ma]);
    %        end
    %    end
    %    
    %    if length(pulse) ~= 0
    %        hold on;
    %        xline(pulse);
    %    end

    %    sgtitle([area_data{area_i}.name ' area' num2str(area_i)], 'Interpreter', 'none');
    %    savefig(gcf, [savepath area_data{area_i}.name ' area' num2str(area_i) ' _raster.fig']);
    %end
%%
% This end is needed to recurse through all of the create list paths
%end 

%all_subVm = [];
%all_rasters = [];
%
%for trial=trial_nums
%    trial
%    cur_trace = data.result.traces(data.result.trial_vec == trial);    
%    spike_info = spike_detect_SNR_sim3(cur_trace, 3.75, 4, 7);
%
%    % Store the subthreshold Vm
%    [baseline, coeff] = Multi_func.exp_fit_Fx(spike_info.trace_ws', round(Fs));
%    cur_subVm = spike_info.trace_ws - baseline;
%    %cur_subVm = cur_trace';
%    all_subVm = horzcat_pad(all_subVm, cur_subVm');
%
%    % Store the spike raster
%    cur_raster = spike_info.roaster;
%    all_rasters = horzcat_pad(all_rasters, cur_raster');
%end
%
%% Show the average sub Threshold
%figure;
%timeline = (1:size(all_subVm, 1))./Fs - 1;
%plot(timeline, mean(all_subVm, 2, 'omitnan'));
%title('Plotting average subthreshold Vm');
%savefig(gcf, [path files(1:end - 5) '_subVmAvg.fig']);
%
%% Plot the spike rasters
%figure;
%timeline = (1:size(all_subVm, 1))./Fs - 1;
%for i=1:size(all_rasters, 2)
%    %cur_color = [rand, rand, rand];
%    cur_color = [0 0 0];
%    spike_idx = find(all_rasters(:, i) == 1);
%
%    if length(spike_idx) == 0
%        continue;
%    end
%
%    %plot(timeline(spike_idx), repmat(i*5, length(spike_idx), 1), '.', 'Color', cur_color);
%    plot(timeline(spike_idx), i*5, '.', 'Color', cur_color);
%    hold on;
%end
%title('Plotting spike raster');
%savefig(gcf, [path files(1:end - 5) '_spike_raster.fig']);
%
%% Plot all of the traces
%figure;
%surface(timeline, 1:size(all_subVm, 2), all_subVm', 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%title('All Traces Heatmap');
%savefig(gcf, [path files(1:end - 5) '_trials_heatmap.fig']);
%
%%TODO Show the subthreshold spectra
%Fs = 828;
%[wt, f] = get_power_spec(mean(all_subVm, 2, 'omitnan'), Fs);
%figure;
%surface(timeline, ... 
%        f, ...
%        abs(wt), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
%title('Power spectra of average signal');
%savefig(gcf, [path files(1:end - 5) '_SpecAvg.fig']);
%
%% Calculate cwt for input signal and 
%function [wt, f] = get_power_spec(signal, samp_freq)
%    freqLimits = [0 150];
%    fb = cwtfilterbank(SignalLength=length(signal),...
%                       SamplingFrequency=samp_freq,...
%                       FrequencyLimits=freqLimits);
%    [wt, f] = cwt(signal, FilterBank=fb);
%end
