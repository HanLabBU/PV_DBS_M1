clc;
clear all;
close all;
f = filesep;

%%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

exclude_200ms = 1;

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

% Parameters for frames to chop off
if ~exclude_200ms
    front_frame_drop = 15;
else 
    front_frame_drop = 15 + round((828*.200));
end

back_frame_drop = 2496;

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
%pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
%
%% CSV file to determine which trials to ignore
%ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
%                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% USER make sure this path changes based on the above line
save_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Traces' f];

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end

set(0,'DefaultFigureVisible','off');

%Load the data
load(save_all_data_file);

% -- Depracated --
% Check if combining all of the regions or not
%if all_regions == 1
%    region_data = Multi_func.combine_regions_old(region_data);
%end

if all_regions == 1
    region_data = Multi_func.combine_regions(region_data);
end

field1 = fieldnames(region_data);
field1 = field1(1);
avg_Fs = mean(region_data.(field1{1}).f_40.framerate, 'omitnan');
timeline = ( (4+(front_frame_drop:back_frame_drop) )./avg_Fs) - 1;

% Loop through each raw trace, detrended trace, and movement if it is there
for f_region = {'r_M1'}%fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    for f_stim=stims'
        f_stim = f_stim{1};

        for nr=1:length(data_bystim.(f_stim).neuron_hilbfilt)

            
            figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
            tiledlayout(size(data_bystim.(f_stim).all_trial_SubVm{nr}, 2), 1, 'TileSpacing', 'compact', 'Padding', 'compact'); %  'Units', 'centimeters', 'InnerPosition', [4, 20, 10.68, 3.0886]

            timeline = data_bystim.(f_stim).trace_timestamps(:, nr);
            
            for tr=1:size(data_bystim.(f_stim).all_trial_SubVm{nr}, 2)
                nexttile;
                % Plot the raw trace
                rawvm = data_bystim.(f_stim).all_trial_rawVm{nr}(:, tr);
                subvm = data_bystim.(f_stim).all_trial_SubVm{nr}(:, tr);
                plot(timeline, rawvm, 'k');
                hold on;
                filt_vm = hilb_filt_range(subvm, [2 10], avg_Fs);
                plot(timeline, subvm, 'r');
            end
            sgtitle(data_bystim.(f_stim).neuron_name{nr}, 'Interpreter', 'none');
            saveas(gcf, [save_path data_bystim.(f_stim).neuron_name{nr} ' tr: ' num2str(tr) '.png']);
        end
    end
end
   
function [filt_sig] = hilb_filt_range(sig, range, FS)
    Fn = FS/2;
    FB = [0.8 1.2].*range;
    
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    filt_sig = hilbert(filtfilt(B,A,sig));
end
