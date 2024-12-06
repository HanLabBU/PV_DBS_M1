% This script will perform spike detection with the specified threshold and append it as a new spike_info<threshold>
f = filesep;

%USER modify everything below here
%% Specific to the computer being used
%DMD
%local_root_path = 'Z:\';

%Dual scope
%local_root_path = '\\engnas.bu.edu\research\eng_research_handata\';

% Maingear office computer
local_root_path = '~/Projects/';
server_root_path = '~/handata_server/';

% Eric's scripts path
%addpath(genpath([server_root_path 'EricLowet']));

% Folder location of all the saved and aligned data
data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Data on local computer
%data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

%USER END modification

matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Path for spike detection debug folder
debug_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Spikes_check' f];

%% Load all of the stim data
save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];

%Load the data
load(save_all_data_file);

%% TODO plot the SBRs for each neuron and compare across brain region and what not if possible

%% Plot all of the SBRs of each neuron for each condition

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Find all modulated neurons
        mod_nr = find(sum(popul_data.mod_matrix, 2) > 0);

        figure('Position', [0 0 2000 1500]);
        data = [];
        labels = [];
        % Loop through each neuron
        for nr = mod_nr' 
            sbrs = popul_data.all_trial_spike_amp{nr}./popul_data.all_trial_trace_noise{nr};
            sbrs(isnan(sbrs)) = [];
            data = [data, sbrs(:)'];
            labels = [labels, repmat({num2str(nr)}, 1, length(sbrs(:)'))];
        end

        % Make the violin plot
        ViolinOpts = Multi_func.get_default_violin();
        ViolinOpts.QuartileStyle = 'shadow';
        group_order = unique(labels, 'stable');
        violins = violinplot(data, labels, 'GroupOrder', group_order, ViolinOpts);
        xlabel('Neuron #');
        ylabel('SBR');

        sgtitle([f_region ' ' f_stim], 'Interpreter', 'none');
    end
end

%% Detect spikes for individual matfiles
% Loop through each FOV matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);
    
    % Loop through each trial
    for  tr=1:length(data.align.trial)
        if ~isempty(data.align.trial{tr})
            %% Use DMD spike detection script
            %[spike_info] = spike_detect_SNR_v3b(data.raw.trial{tr}.raw_traces, 3.75);
            %data.align.trial{tr}.spike_info375 = spike_info;
            
            % Use Eric's good spike detection
            [spike_info] = spike_detect_SNR_sim3(data.raw.trial{tr}.raw_traces, 3.75, 4, 7);
            data.align.trial{tr}.spike_info375 = spike_info;
        end
    end
    
    % Update the matfile
    align = data.align;
    save([data_path matfile_names{i}], 'align', '-append');
    
    % Set the title of the Figure
    %sgtitle(matfile_names{i}, 'Interpreter', 'none');
end

%% Sweep through each SBR threshold to determine what would be the best spike SBR for neurons
thresholds=3:0.1:4;

% TODO I think I should use the mod_matrix to sift through which neurons to actually use

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Find all modulated neurons
        mod_nr = find(sum(popul_data.mod_matrix, 2) > 0);

        % Loop through each neuron
        for nr = mod_nr' % TODO change to all neurons 
            raw_vm = popul_data.all_trial_rawVm{nr}; 
            norm_vm = Multi_func.norm_signals(popul_data.all_trial_rawVm{nr});
            
            sp_idx_mat = {};
            % Loop through each threshold
            for thres = thresholds %TODO remove the paranetheses
                [sp_info] = spike_detect_SNR_sim3(raw_vm, thres, 4, 7);
                sp_idx_mat = [sp_idx_mat, sp_info.spike_idx];
            end
            
            % Loop through each trial and plot for each figure
            for tr_i = 1:size(raw_vm, 2)
                figure('Position', [0 0 2000 1500]);
                
                %loop through each threshold
                for thres_i = 1:length(thresholds)
                    plot(norm_vm(:, tr_i) + thres_i, '-k');
                    hold on;
                    plot( sp_idx_mat{tr_i, thres_i}, ...
                        norm_vm(sp_idx_mat{tr_i, thres_i}, tr_i) + thres_i, 'or');
                    hold on;
                    text(min(-1), max(norm_vm(:, tr_i) + thres_i), ...
                        num2str(thresholds(thres_i)), ...
                        'FontSize', 14, 'Color', 'g', 'FontWeight', 'bold');
                end
                
                Multi_func.set_default_axis(gca);

                % Set title figure
                sgtitle([popul_data.neuron_name{nr} ' tr ' num2str(tr_i)], 'Interpreter', 'none');
                saveas(gcf, [debug_path popul_data.neuron_name{nr} ' tr ' num2str(tr_i) '.png']);
            end
        end
    end
end

%% TODO plot the SBRs for each neuron and compare across brain region and what not if possible


%% Loop through by each matfile
%Loop through each matfile
for i=1 % TODO include all of the other matfiles :length(matfile_names)
    data = load([data_path matfile_names{i}]);
    
    % Create figure for each neuron
    figure('Position', [0 0 2000 1500]);
    tiledlayout(1, length(thresholds));

    % Sweep through each spike detection threshold
    for thres=thresholds
        nexttile;
        % Loop through each trial
        for  tr=1:length(data.align.trial)
            if ~isempty(data.align.trial{tr})
                
            end
        end
    end
end
