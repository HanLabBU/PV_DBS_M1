% Plot the currents used in the experiments
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
savepath = Multi_func.save_plot;

exclude_200ms = 1;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

all_currents = [];
all_Fs = [];

% Use congregated data to save all of the currents
% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end
%Load the data
load(save_all_data_file);

% Loop through and grab each neuron's current amplitude
for f_region = {'r_M1'}
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

    for f_stim = fieldnames(data_bystim)'
        f_stim = f_stim{1};

        % Loop through each neuron
        for nr = 1:length(data_bystim.(f_stim).neuron_name)
            cur_nr_name = strsplit(data_bystim.(f_stim).neuron_name{nr}, '_');
            
            % Get the current from the matfile
            f_mouse = ['m_' cur_nr_name{1}];
            f_rec = ['r_' cur_nr_name{3}];               
            f_neuron= ['n_' cur_nr_name{4}];
            f_stim = ['f_' cur_nr_name{5}];

            neuron_num = str2num(erase(cur_nr_name{4}, 'FOV'));
            
            [region_data] = create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
            
            % Saving the current at neuron number index
            region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(cur_nr_name{6});
            
        end
    end
end

%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%
%    % Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    % Select matfiles by stim condition for given region
%    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);
%
%    % Loop through each stimulation condition
%    for f_stim = fieldnames(matfile_stim)'
%        f_stim = f_stim{1};
%        matfiles = matfile_stim.(f_stim).names;    
%    
%        % Initialize field subthreshold array
%        data_bystim.(f_stim) = struct();
%
%        % Loop through each matfile of the current stimulation condition
%        for matfile = matfiles
%            matfile = matfile{1};
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%         
%            
%            % Loop through each ROI
%            for roi_idx=1:size(trial_data.detrend_traces, 2)
%                %Determine whether this roi is to be ignored for this particular trial
%                ri = strsplit(matfile, '_');
%                try
%                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
%                catch
%                    trial_ignr_list = [];
%                end
%                
%                % Remove trials from trial idx list
%                trial_idxs = setdiff(trial_idxs, trial_ignr_list);
%
%                % Skip if there is at most 2 trials
%                if length(trial_idxs) <= 2
%                    continue;
%                end
%                
%                % Get the current from the matfile
%                f_mouse = ['m_' ri{1}];
%                f_rec = ['r_' ri{3}];               
%                f_neuron= ['n_' ri{4}];
%                f_stim = ['s_' ri{5}];
%                
%                % convert FOV number to single number
%                neuron_num = str2num(erase(ri{4}, 'FOV'));
%
%                [region_data] = create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
%                
%                % Saving the current at neuron number index
%                region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(ri{6});
%
%                % Save all of the currents
%                all_currents(end + 1) = str2num(ri{6});
%
%            end
%        end
%    end
%end



% Clear all zeros in the current array
all_currents(find(all_currents == 0)) = [];

disp(['Average currents ' num2str(mean(all_currents))]);
disp(['Std currents ' num2str(std(all_currents))]);

% Plot all of the currents as a single plot
figure;
tick_labels = {};
tick_spots = [];
i = 1;
for f_region = {'r_M1'} %fieldnames(region_data)'
    f_region = f_region{1};

    for f_mouse = fieldnames(region_data.(f_region))'
        f_mouse = f_mouse{1};
        
        % Fields with leading m indicate a mouse entry
        if strcmp(f_mouse(1), 'm') ~= 1
            continue;
        end
    
        for f_rec = fieldnames(region_data.(f_region).(f_mouse))'
            f_rec = f_rec{1};
 
            % Check if each stimulation condition exists for the neuron
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_40')
                current = region_data.(f_region).(f_mouse).(f_rec).f_40.currents;
                current(find(current == 0)) = [];
                
                plot(i + rand(1, length(current)), current, '.r');
                hold on;
            end
            i = i + 1;
            
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_140')
                current = region_data.(f_region).(f_mouse).(f_rec).f_140.currents;
                current(find(current == 0)) = [];
                plot(i + rand(1, length(current)), current, '.b');
                hold on;
            end

            % Add the mouse and recording date for the plot
            tick_labels{end + 1} = [erase(f_mouse, 'm_')]; % ' ' erase(f_rec, 'r_')
            tick_spots(end + 1) = i;

            i = i + 2;
        end 

        i = i + 3;
    end
end
xticks(tick_spots);
xticklabels(tick_labels);
ylabel('uamp');

saveas(gcf, [savepath 'Current/' 'all_mice_currents.png']);
saveas(gcf, [savepath 'Current/' 'all_mice_currents.pdf']);

% Plot the Normalized Vm Change vs Current amplitude

% Create fields in structure if they do not exist
function [result] = create_struct(src, f_region, f_mouse, f_rec, f_neuron, f_stim)
    result = src;
    % Create region field
    if isfield(result, f_region) == 0
        result.(f_region) = struct();
    end
    
    % Create mouse field
    if isfield(result.(f_region), f_mouse) == 0
        result.(f_region).(f_mouse) = struct();
    end

    % Create recording field
    if isfield(result.(f_region).(f_mouse), f_rec) == 0
        result.(f_region).(f_mouse).(f_rec) = struct();
    end

    % Create stim field
    if isfield(result.(f_region).(f_mouse).(f_rec), f_stim) == 0
        result.(f_region).(f_mouse).(f_rec).(f_stim) = struct();
        result.(f_region).(f_mouse).(f_rec).(f_stim).currents = [];
    end

    % Create neuron field
    % Old code for saving everything into the neuron field's
    %if isfield(result.(f_region).(f_mouse).(f_rec).(f_stim), f_neuron) == 0
    %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron) = struct();
    %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron).currents = [];
    %    
    %end
end
