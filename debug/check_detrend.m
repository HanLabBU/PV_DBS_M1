%% Save all analysis data to a matfile that can be easily read for analysis
clear all;
close all;

f = filesep;

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

exclude_200ms = 1;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else
    front_frame_drop = 15 + round((828*.200));
end
back_frame_drop = 2496;

% Time periods for comparison of firing rate and sub Vm
base_ped = Multi_func.base_ped;
trans_ped =Multi_func.trans_ped;
sus_ped = Multi_func.sus_ped;
stim_ped = Multi_func.stim_ped;
offset_trans_ped = Multi_func.offset_trans_ped;

% Data path for all of the intermediate anaylsis data
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Check that the server path exists
if ~isfolder(local_root_path) || ~isfolder(pv_data_path)
    disp('Server rootpath does not exist!!!');
    return;
end

save_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Debug_detrend' f];

%% Setup
ses = dir([pv_data_path '*.mat']); %DEBUG for debugging just try V1 data
all_matfiles = {ses.name};

% Select matfiles by brain region
[region_matfiles] = Multi_func.find_region(all_matfiles);
region_data = struct();
all_Fs = [];
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    % Select matfiles by stim specific conditions for all regions
    %[matfile_stim] = stim_cond(all_matfiles); 
    % Select matfiles by stim condition for given region
    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);

    % Loop through each field of the struct and concatenate everything together
    % Store trace aspect data by each individual stimulation condition
    data_bystim = struct();
    % Store all of the calculated sampling frequencies

    % Loop through each stimulation condition
    for f_stim = fieldnames(matfile_stim)'
        f_stim = f_stim{1};
        matfiles = matfile_stim.(f_stim).names;    

        all_base_only_exp_vm = [];
        all_base_offset_exp_vm = [];
        all_filt_vm = [];
        all_base_off_filt_vm = [];
        all_base_lin_vm = [];
        
        all_base_only_exp_fit = [];
        all_base_offset_exp_fit = [];
        all_filter_fit = [];
        all_base_offset_filter_fit = [];
        all_base_lin_fit = [];

        for matfile = matfiles %(1:2) % Testing only a few
            matfile = matfile{1};
            % Read in the mat file of the current condition
            data = load([pv_data_path matfile]);
            trial_idxs = find(~cellfun(@isempty, data.align.trial));
            trial_data = data.align.trial{trial_idxs(1)};    

            % Loop through each ROI
            for roi_idx=1:size(trial_data.detrend_traces, 2)
                %Determine whether this roi is to be ignored for this particular trial
                ri = strsplit(matfile, '_');
                try
                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                catch
                    trial_ignr_list = [];
                end
                
                % Remove trials from trial idx list
                trial_idxs = setdiff(trial_idxs, trial_ignr_list);

                % Skip if there is at most 2 trials
                if length(trial_idxs) <= 2
                    continue;
                end
                
                % Loop through each trial                
                for tr_idx=trial_idxs %(1:2)
                    trial_data = data.align.trial{tr_idx};
                    raw_trial_data = data.raw.trial{tr_idx};

                    % Store all of the timestamp info
                    stim_start = raw_trial_data.raw_stimulation_time(1);

                    % Test different detrending methods
                    if isfield(trial_data, 'interp_time')
                        cur_trace_time = trial_data.interp_time(front_frame_drop:back_frame_drop) - stim_start;
                    else
                        cur_trace_time = trial_data.camera_frame_time(front_frame_drop:back_frame_drop) - stim_start;
                    end

                    % Check the raw traces for the detrending methods
                    if isfield(trial_data, 'interp_time')
                        % Interpolate the raw traces
                        cur_trace_ws = raw_trial_data.interp_raw_traces(front_frame_drop:back_frame_drop, roi_idx)';
                    else
                        cur_trace_ws = raw_trial_data.raw_traces(front_frame_drop:back_frame_drop, roi_idx)';
                    end
                    
                    %cur_trace_ws = trial_data.spike_info375.trace_ws(roi_idx, front_frame_drop:back_frame_drop);
                        
                    [base_offset_exp, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));

                    [base_only_exp, coeff] = exp_fit_new(cur_trace_ws', round(trial_data.camera_framerate));
                    
                    [base_only_lin_fit] = lin_fit(cur_trace_ws, round(trial_data.camera_framerate));

                    trace_filter = fastsmooth(cur_trace_ws, 2000, 1, 1);
                    
                    base_offset_filter = mov_avg_adj(cur_trace_ws, round(trial_data.camera_framerate));

                    % Test different detrending methods
                    figure('Position', [1, 1, 1000, 900]);
                    tiledlayout(5, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters');
                    nexttile;
                    plot(cur_trace_time, cur_trace_ws); % Raw trial
                    hold on;
                    
                    plot(cur_trace_time, base_offset_exp, 'DisplayName', 'Base and Offset exponential');
                    hold on;
                    plot(cur_trace_time, base_only_lin_fit, 'DisplayName', 'Base only linear fit');

                    % Plotting the exponential fit from just the baseline
                    %plot(cur_trace_time, base_only_exp, 'DisplayName', 'Baseline only exponential');
                    hold on;
                    plot(cur_trace_time, trace_filter, 'g', 'DisplayName', 'Filter detrending method');
                    hold on;
                    plot(cur_trace_time, base_offset_filter, 'm', 'DisplayName', 'Filtering around points method');
                    legend();
                    title([f_region ' ' f_stim ' ' matfile], 'Interpreter', 'none');
                    

                    nexttile;
                    base_offset_exp_detrend_trace = cur_trace_ws - base_offset_exp;
                    plot(cur_trace_time, base_offset_exp_detrend_trace);
                    title('Old method (baseline and offset points)');
                    
                    nexttile;
                    base_only_lin_detrend_trace = cur_trace_ws - base_only_lin_fit;
                    plot(cur_trace_time, base_only_lin_detrend_trace);
                    title('Linear fit method (baseline points)');
 
                    % Baseline only exponential fit
                    base_only_exp_detrend_trace = cur_trace_ws - base_only_exp;
                    %plot(cur_trace_time, base_only_exp_detrend_trace);
                    %title('New method (baseline only points)');
                    
                    nexttile;
                    filt_detrend_trace = cur_trace_ws - trace_filter;
                    plot(cur_trace_time, filt_detrend_trace);
                    title('Filter method (all points)');

                    nexttile;
                    base_off_filt_detrend_trace = cur_trace_ws - base_offset_filter;
                    plot(cur_trace_time, base_offset_exp_detrend_trace);
                    title('Filter method (base and offset points)');

                    % Save to specific location to check on all matfiles
                    saveas(gcf, [save_path 'raw_' matfile '_roi' num2str(roi_idx) '_tr' num2str(tr_idx) '.png']);
                
                    % Save all of the Vms to respective arrays
                    all_base_only_exp_vm = horzcat_pad(all_base_only_exp_vm, base_only_exp_detrend_trace(:));
                    all_base_offset_exp_vm = horzcat_pad(all_base_offset_exp_vm, base_offset_exp_detrend_trace(:));
                    all_filt_vm = horzcat_pad(all_filt_vm, filt_detrend_trace(:));
                    all_base_off_filt_vm = horzcat_pad(all_base_off_filt_vm, base_off_filt_detrend_trace(:));
                    all_base_lin_vm = horzcat_pad(all_base_lin_vm, base_only_lin_detrend_trace(:));

                    all_base_only_exp_fit = horzcat_pad(all_base_only_exp_fit, base_only_exp(:));
                    all_base_offset_exp_fit = horzcat_pad(all_base_offset_exp_fit, base_offset_exp(:));
                    all_filter_fit = horzcat_pad(all_filter_fit, trace_filter(:));
                    all_base_offset_filter_fit = horzcat_pad(all_base_offset_filter_fit, base_offset_filter(:));

                    all_base_lin_fit = horzcat_pad(all_base_lin_fit, base_only_lin_fit(:));
			
                end % End looping through trial
            end % End looping through ROI
        end % End looping through individual FOVS
        
        data_bystim.(f_stim).all_base_only_exp_vm = all_base_only_exp_vm;
        data_bystim.(f_stim).all_base_offset_exp_vm = all_base_offset_exp_vm;
        data_bystim.(f_stim).all_filt_vm = all_filt_vm;
        data_bystim.(f_stim).all_base_off_filt_vm = all_base_off_filt_vm;
        data_bystim.(f_stim).all_base_lin_vm = all_base_lin_vm;

        data_bystim.(f_stim).all_base_only_exp_fit = all_base_only_exp_fit;
        data_bystim.(f_stim).all_base_offset_exp_fit = all_base_offset_exp_fit;
        data_bystim.(f_stim).all_filter_fit = all_filter_fit;
        data_bystim.(f_stim).all_base_offset_filter_fit = all_base_offset_filter_fit;
        data_bystim.(f_stim).all_base_lin_fit = all_base_lin_fit;

    end % End looping stim conditions
    region_data.(f_region) = data_bystim;
end % End looping through regions
save('detrend_debug.mat', 'region_data');

%% Plot all of the photobleach averages
figure;
tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters');
nexttile;
base_avg = mean(all_base_only_exp_vm, 2, 'omitnan');
base_std = std(all_base_only_exp_vm, 0, 2, 'omitnan');
num_vm = size(all_base_only_exp_vm, 2);
base_sem = base_std./sqrt(num_vm);
timeline = [1:size(all_base_only_exp_vm, 1)]';
fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;
plot(timeline, base_avg, 'k', 'LineWidth', 0.3);

nexttile;
base_avg = mean(all_base_offset_exp_vm, 2, 'omitnan');
base_std = std(all_base_offset_exp_vm, 0, 2, 'omitnan');
num_vm = size(all_base_offset_exp_vm, 2);
base_sem = base_std./sqrt(num_vm);
timeline = [1:size(all_base_offset_exp_vm, 1)]';
fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;
plot(timeline, base_avg, 'k', 'LineWidth', 0.3);


nexttile;
base_avg = mean(all_filt_vm, 2, 'omitnan');
base_std = std(all_filt_vm, 0, 2, 'omitnan');
num_vm = size(all_filt_vm, 2);
base_sem = base_std./sqrt(num_vm);
timeline = [1:size(all_filt_vm, 1)]';
fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
Multi_func.set_fill_properties(fill_h);
hold on;
plot(timeline, base_avg, 'k', 'LineWidth', 0.3);

saveas(gcf, 'Raw_Avg_SEM_Traces.png');
%%
function [ fitbaseline]=mov_avg_adj(v,FS)    
    v1=v;    
    v1([FS-10:FS*2+20])= mean([ v([FS-40:FS-20 2*FS+20:2*FS+40 ])]);    
    v1(1)=v1(2);

    figure('Position', [400, 400, 800, 500]);
    plot(v, 'DisplayName', 'Original Trace');
    hold on;
    plot(v1, 'DisplayName', 'Stimulation period Interpolated');
    legend();

    fitbaseline = fastsmooth(v1, 2000, 1, 1);
end

function [baseline] = lin_fit(v, Fs)
    t = 1:length(v);
    tsel = t(1:Fs);
    
    %
    p = polyfit(t(tsel), v(tsel), 1);
    slope = p(1);
    intercept = p(2);

    baseline = slope.*t + intercept;
end

function [ fitbaseline, coeff]=exp_fit_new(v,FS)    
    %% fits an expoential to estimate the photobleaching    
    v1=v;    
    v1([FS-10:FS*2+20])= mean([ v([FS-40:FS-20 2*FS+20:2*FS+40 ])]);    
    v1(1)=v1(2);    
    F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3)); %+ x(3)*exp(- xdata./x(4))  ;    
    x0 = [mean(v1) 40 1.5] ;    
    OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');    
    t = (1:length(v))./FS;
    % Whole trace
    %tsel=[1:length(t)];    
    %Base only, --- and offset period consideration
    tsel=[1:FS ];    
    [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);    
    fitbaseline=F(xunc, t);    
    %DEBUG
    %figure; tiledlayout(2, 1); nexttile; plot(v); hold on; plot(fitbaseline);
    %nexttile;
    %plot(v' - fitbaseline);
    coeff=xunc;    
end
