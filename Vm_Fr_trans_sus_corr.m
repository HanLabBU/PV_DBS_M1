clear all;
close all;
f = filesep;

%% USER Modification
% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

exclude_200ms = 1;

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
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

%figure_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'Plots' f];
figure_path = Multi_func.save_plot;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Parameter to determine whether to combine all regions as one data
all_regions = 0;

%% END Modification

%% Check that the server path exists
if ~isfolder(server_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end
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

%% Vm vs Firing Rate correlation for whole stimulation period
% Loop through each region 
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    figure('visible', 'on', 'Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
    tiledlayout(length(stims), 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 12, 12]);

    % Loop through stim condition
    for f_stim = stims'
        f_stim = f_stim{1};
        cond_data = data_bystim.(f_stim);

        % Stim depolarization zscored modulated to baseline subtracted
        % Could zscore to the mean of each factor
        %Vm_change = zscore(cond_data.neuron_stim_Vm);
        %Fr_change = zscore(cond_data.neuron_stim_FR);
        Vm_change = cond_data.neuron_stim_Vm;
        Fr_change = cond_data.neuron_stim_FR;


        nexttile;
        % Plot zscored change between firing rate and Vm
        plot(Vm_change, Fr_change, '.');
        hold on;
        
        % Calculate the correlation value
        fitresult = polyfit(Vm_change, Fr_change, 1);
        x_fit = Vm_change;
        y_fit = polyval(fitresult, x_fit);
        plot(x_fit, y_fit);
        SSR = sum((Fr_change - y_fit).^2);
        SST = sum((Fr_change - mean(y_fit)).^2);
        r2 = 1 - (SSR/SST);
        r = sqrt(r2);

        % TODO save correlation value into legend
        legend(num2str(r));
        xlabel('Vm Change');
        ylabel('Firing Rate Change');
        title([f_stim(3:end) 'Vm change vs FR change'], 'Interpreter', 'none');
    end
end
