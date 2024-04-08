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
all_base_only_vm = [];
all_base_offset_vm = [];
all_filt_vm = [];
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
                        
                    [base_offset, coeff] = Multi_func.exp_fit_Fx(cur_trace_ws', round(trial_data.camera_framerate));

                    [base_only, coeff] = exp_fit_new(cur_trace_ws', round(trial_data.camera_framerate));
                    
                    base_filter = fastsmooth(cur_trace_ws, 2000, 1, 1);

                    % Test different detrending methods
                    figure('Position', [1, 1, 1000, 900]);
                    tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters');
                    nexttile;
                    plot(cur_trace_time, cur_trace_ws); % Raw trial
                    hold on;
                    plot(cur_trace_time, baseline, 'DisplayName', 'Old detrending method');
                    hold on;
                    plot(cur_trace_time, base_new, 'DisplayName', 'New detrending method');
                    hold on;
                    plot(cur_trace_time, base_filter, 'g', 'DisplayName', 'Filter detrending method');
                    legend();
                    title([f_region ' ' f_stim ' ' matfile], 'Interpreter', 'none');
                    nexttile;
                    base_offset_detrend_trace = cur_trace_ws - base_offset;
                    plot(cur_trace_time, base_offset_detrend_trace);
                    title('Old method (baseline and offset points)');
                    nexttile;
                    base_only_detrend_trace = cur_trace_ws - base_only;
                    plot(cur_trace_time, base_only_detrend_trace);
                    title('New method (baseline only points)');
                    nexttile;
                    filt_detrend_trace = cur_trace_ws - base_filter;
                    plot(cur_trace_time, filt_detrend_trace);
                    title('Filter method (all points)');

                    % Save to specific location to check on all matfiles
                    saveas(gcf, [save_path 'raw_' matfile '_roi' num2str(roi_idx) '_tr' num2str(tr_idx) '.png']);
                
                    % Save all of the Vms to respective arrays
                    all_base_only_vm = horzcat_pad(all_base_only_vm, base_only_detrend_trace(:));
                    all_base_offset_vm = horzcat_pad(all_base_only_vm, base_offset_detrend_trace(:));
                    all_filt_vm = horzcat_pad(all_filt_vm, filt_detrend_trace(:));
                end % End looping through trial
            end % End looping through ROI
        end % End looping through individual FOVS
    end % End looping stim conditions
end % End looping through regions



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
