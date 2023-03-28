% Housekeeping clear all stuff
clc;
clear all;
close all;

f = filesep;

%------------ USER modification
% Server root path
%server_rootpath = 'Z:\';

% Maingear office computer
local_rootpath = '/home/pierfier/Projects/';

server_rootpath = ['~/handata_server' f 'eng_research_handata3' f];

% Server folder location of the saved and aligned data
%data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Local linux machine
%data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data share on server
data_path = ['~/handata_server' f 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Determine whether to save the ignored traces (1) or not ignored traces (0)
show_ignored = 1;
% USER make sure this path changes based on the above line
%save_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'Ignored_Traces' f];

figure_path = [server_rootpath 'Pierre Fabris' f 'PV Project' f 'Figures' f 'Raster plots' f];

% Filepath name for ignoring individual trial csv
ignore_trial_csv = [local_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv'];


%------------- END modification

% Get all of the trials that are being ignored
ignore_trial_dict = Multi_func.csv_to_struct(ignore_trial_csv);

% Get all FOV matfiles
matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Grab only 140Hz FOVs
[mat_struct] = Multi_func.stim_cond(matfile_names);
matfile_names = mat_struct.f_40.names;

% Grab only the M1 FOVs from this list
[region_struct] = Multi_func.find_region(matfile_names);
matfile_names = region_struct.r_V1.names;

% Setup figure to show alignment data for all trials
figure('Position', [0, 0, 800, 1000]);
posy = 0;
% Loop through each matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);

    %DEBUG
    disp(matfile_names{i});
        
    % Set color for current FOV
    cur_color = [rand, rand, rand]*0.7;
    
    % Loop through each trial
    for tr=1:length(data.align.trial)
        if isempty(data.align.trial{tr})
            continue;
        end

        % Check if trial is in the ignore list
        ri = strsplit(matfile_names{i}, '_');
        try
            trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).('ROI1'); % Always assuming ROI 1
        catch
            trial_ignr_list = [];
        end

        % Check if current trial is in the ignore list
        if ismember(tr, trial_ignr_list) ~= show_ignored
            continue;
        end

        %tiledlayout((size(data.align.trial{tr}.detrend_traces, 2)*2) + 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

        % Save current trial
        trial_data = data.align.trial{tr};
        % Loop through each ROI
        for roi=1:size(trial_data.detrend_traces, 2)
            
            % Check if trace has NaNs or contains spikes
            if sum(~isnan(data.raw.trial{tr}.raw_traces(:, roi))) == 0 || isempty(trial_data.spike_info375.spike_idx{1})
                continue;
            end
            
            % Plotting each raster
            plot(trial_data.spike_info375.spike_idx{1}, posy, '.', 'Color', cur_color);
            hold on;
            posy = posy + 1;
        end
    end
    %yline(posy - 0.5, 'k--');
    %hold on;
end
%set(gca,'color','none');
title('40_V1_raster', 'Interpreter', 'none');
saveas(gcf, [figure_path '40_V1_raster.png']);

% Perform exponential fit
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
