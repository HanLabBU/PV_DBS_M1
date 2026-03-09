clear all;
close all;

f = filesep;

%%%---- USER Modification--------
% Linux server
local_root_path = '~/Projects/';
server_root_path = '~/handata_server/eng_research_handata3/';

% Windows server
%local_root_path = 'Z:\';

% Save the figure path
figure_path = Multi_func.save_plot;

% Data on server
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Specify which ignore trials to use
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);
%%%----- END Modification---------

%% Count the total number of mice used for each brain region
% This section requires 'plot_current.m' ran with the additional flicker experimental values

% Loop through
for f_region=fieldnames(region_data)'
    f_region = f_region{1};
    
    % Get all of the region neuron names
    all_nr_names = {};

    % Loop through each stimulation condition
    stims=fieldnames(region_data.(f_region))';
    for f_stim=stims
        f_stim = f_stim{1};
        popul_data = region_data.(f_region).(f_stim);
        
        all_nr_names = cat(1, all_nr_names, popul_data.neuron_name(:));

    end
    mouse_names = cellfun(@(x) strtok(x, '_'), all_nr_names, 'UniformOutput', false);

    % Display the unique number of mice for each brain region
    disp(['Region ' f_region]);
    disp(['Num mice: ' num2str(length(unique(mouse_names) ) )]);
    disp(['Num total recordings: ' num2str(length(mouse_names))]);
    
    %DEBUG displaying all of the mice names
    %sort(mouse_names)
end

%% Create table showing all estim-only neurons from the regular stimulation experiments
estim_nr_t = table();

% Have the relabels ready
mouse_rename = struct();
mice_names = fieldnames(Multi_func.mouse_color)';
for f_i = 1:length(mice_names)
    mouse_rename.(mice_names{f_i}) = ['Mouse ' num2str(f_i)];
end

nr_num = 1;
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequenciesd
    for f_stim = stims'

        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Loop through each neuron
        for nr_i = 1:length(popul_data.neuron_name(:)')
            nr_i;
            
            % Parse out recording info
            name_parts = strsplit(popul_data.neuron_name{nr_i}, "_");
            mouse_name = name_parts{1};

            % Find and grab the current amperage value
            idxs = contains(name_parts, 'mat');
            
            % Current amplitude is attached to the extension
            if ~isempty(regexp(name_parts{idxs}, '\d', 'once'))
                name_parts{idxs};
                current = strsplit(name_parts{idxs}, '.');
                current = current{1};
            else
                % Current is the previous entry from the "_" split
                idxs = find(idxs == 1);
                current = name_parts{idxs - 1};
            end

            % Enter the neuron identification info
            estim_nr_t{num2str(nr_num), 'Mouse_ID'} = string(mouse_rename.(['m' mouse_name]));
            %estim_nr_t(num2str(nr_num), 'Neuron_ID') = {popul_data.neuron_name{nr_i}};
            estim_nr_t{num2str(nr_num), 'Brain Region'} = string(f_region(3:end));
            
            % Enter the stimulation parameters
            estim_nr_t{num2str(nr_num), 'Current'} = str2num(current);
            estim_nr_t{num2str(nr_num), 'Stim Frequency'} = str2num(f_stim(3:end));

            % Enter the statistical significances for each neuron
            temp_plv = [popul_data.plv_mod_stats.obs_PLV2];
            temp_plv_mod = [popul_data.plv_mod_stats.mod];
            %estim_nr_t(num2str(nr_num), 'PLV') = {[num2str(temp_plv(nr_i)) ' (' num2str(temp_plv_mod(nr_i)) ')' ]};
            estim_nr_t{num2str(nr_num), 'PLV Stat'} = temp_plv(nr_i);
            estim_nr_t{num2str(nr_num), 'PLV Sign'} = temp_plv_mod(nr_i);

            temp_vm = [popul_data.Vm_trans_mod_stats.p_sign];
            temp_vm_mod = [popul_data.Vm_trans_mod_stats.mod];
            %estim_nr_t(num2str(nr_num), 'Vm Transient') = {[num2str(temp_vm(nr_i)) ' (' num2str(temp_vm_mod(nr_i)) ')' ]};
            estim_nr_t{num2str(nr_num), 'Vm Transient Stat'} = temp_vm(nr_i);
            estim_nr_t{num2str(nr_num), 'Vm Transient Sign'} = temp_vm_mod(nr_i);

            
            temp_vm = [popul_data.Vm_sus_mod_stats.p_sign];
            temp_vm_mod = [popul_data.Vm_sus_mod_stats.mod];
            %estim_nr_t(num2str(nr_num), 'Vm Sustained') = {[num2str(temp_vm(nr_i)) ' (' num2str(temp_vm_mod(nr_i)) ')' ]};
            estim_nr_t{num2str(nr_num), 'Vm Sustained Stat'} = temp_vm(nr_i);
            estim_nr_t{num2str(nr_num), 'Vm Sustained Sign'} = temp_vm_mod(nr_i);
            
            temp_fr = [popul_data.fr_trans_mod_stats.p_sign];
            temp_fr_mod = [popul_data.fr_trans_mod_stats.mod];
            %estim_nr_t(num2str(nr_num), 'fr Transient') = {[num2str(temp_fr(nr_i)) ' (' num2str(temp_fr_mod(nr_i)) ')' ]};
            estim_nr_t{num2str(nr_num), 'fr Transient Stat'} = temp_fr(nr_i);
            estim_nr_t{num2str(nr_num), 'fr Transient Sign'} = temp_fr_mod(nr_i);

            
            temp_fr = [popul_data.fr_sus_mod_stats.p_sign];
            temp_fr_mod = [popul_data.fr_sus_mod_stats.mod];
            %estim_nr_t(num2str(nr_num), 'fr Sustained') = {[num2str(temp_fr(nr_i)) ' (' num2str(temp_fr_mod(nr_i)) ')' ]};
            estim_nr_t{num2str(nr_num), 'fr Sustained Stat'} = temp_fr(nr_i);
            estim_nr_t{num2str(nr_num), 'fr Sustained Sign'} = temp_fr_mod(nr_i);
            
            % Indicate if neurons were modulated or non-modulated at all
            modulated = sum(popul_data.mod_matrix(nr_i, :));
            estim_nr_t{num2str(nr_num), 'Modulated'} = modulated;

            % Increment neuron number
            nr_num = nr_num + 1;
        end
    end
end

estim_nr_t
writetable(estim_nr_t, [figure_path 'Neuronwise' f 'Neurons_Estim_Only_Data.csv'], 'WriteRowNames', true);

%% Make table while grouped with each mouse, separated by condition
% For the stim only experiments

estim_mouse_t = table();
row_num = 1;
for f_region = fieldnames(region_data)' % TODO change back % {'r_V1'}
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequenciesd
    for f_stim = stims' % {'f_40'} %
        f_stim = f_stim{1};
        
        % Find all the neurons for current condition
        cond_idxs = (estim_nr_t.("Brain Region") == f_region(3:end)) & ...
            (estim_nr_t.("Stim Frequency") == str2num(f_stim(3:end)) );
        
        % Loop through each mouse for each condition
        cond_mice = unique(estim_nr_t.Mouse_ID(cond_idxs));
        for mouse = cond_mice'
            mouse
            mice_idxs = cond_idxs & (estim_nr_t.Mouse_ID == mouse);

            % Loop through each current for this mouse
            cur_amp = unique(estim_nr_t.Current(mice_idxs));
            for amp = cur_amp'
                amp_idxs = mice_idxs & (estim_nr_t.Current == amp);
                
                % Save these group of neuron's identifiable
                estim_mouse_t{row_num, 'Mouse'} = mouse;
                estim_mouse_t{row_num, 'Region'} = f_region(3:end);
                estim_mouse_t{row_num, 'Stim'} = str2num(f_stim(3:end));
                estim_mouse_t{row_num, 'Current'} = amp;
                
                % Count number of significant neurons for PLV
                sig_idx = (estim_nr_t.("PLV Sign")(amp_idxs) > 0);
                estim_mouse_t{row_num, 'PLV Sig'} = sum(sig_idx);
                sig_idx = (estim_nr_t.("PLV Sign")(amp_idxs) < 0);
                estim_mouse_t{row_num, 'Not PLV Sig'} = sum(sig_idx);

                % Count number of significant neurons for Vm transient
                sig_idx = (abs(estim_nr_t.("Vm Transient Sign")(amp_idxs)) > 0);
                estim_mouse_t{row_num, 'Vm Transient Sig'} = sum(sig_idx);
                sig_idx = (estim_nr_t.("Vm Transient Sign")(amp_idxs) == 0);
                estim_mouse_t{row_num, 'Not Vm Transient Sig'} = sum(sig_idx);

                % Count number of significant neurons for Vm sustained
                sig_idx = (abs(estim_nr_t.("Vm Sustained Sign")(amp_idxs)) > 0);
                estim_mouse_t{row_num, 'Vm Sustained Sig'} = sum(sig_idx);
                sig_idx = (estim_nr_t.("Vm Sustained Sign")(amp_idxs) == 0);
                estim_mouse_t{row_num, 'Not Vm Sustained Sig'} = sum(sig_idx);

                % Count number of significant neurons for fr transient
                sig_idx = (abs(estim_nr_t.("fr Transient Sign")(amp_idxs)) > 0);
                estim_mouse_t{row_num, 'fr Transient Sig'} = sum(sig_idx);
                sig_idx = (estim_nr_t.("fr Transient Sign")(amp_idxs) == 0);
                estim_mouse_t{row_num, 'Not fr Transient Sig'} = sum(sig_idx);

                % Count number of significant neurons for Vm sustained
                sig_idx = (abs(estim_nr_t.("fr Sustained Sign")(amp_idxs)) > 0);
                estim_mouse_t{row_num, 'fr Sustained Sig'} = sum(sig_idx);
                sig_idx = (estim_nr_t.("fr Sustained Sign")(amp_idxs) == 0);
                estim_mouse_t{row_num, 'Not fr Sustained Sig'} = sum(sig_idx);

                row_num = row_num + 1;
            end
        end
    end
end
estim_mouse_t = sortrows(estim_mouse_t, {'Region', 'Stim', 'Mouse', 'Current'})
writetable(estim_mouse_t, [figure_path 'Neuronwise' f 'Mouse_Estim_Only_Data.csv'], 'WriteRowNames', true);


%% Print out unique mice
stim_mice = unique(estim_nr_t.Mouse_ID);

% Print out the corresponding brain region
for m=stim_mice'
    m = m{1};
    idx = ismember(estim_nr_t.Mouse_ID, m);
    idx = find(idx == 1);
    idx = idx(1);
    disp([m estim_nr_t.("Brain Region")(idx)]);
end

%% Count average+spread, median+spread for all neurons imaged with TICO
for f_region = fieldnames(region_data)' % TODO change back % {'r_V1'}
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequenciesd
    for f_stim = stims' % {'f_40'} %
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Same mouse neurons
        looped_nrs = [];

        % Count number of neurons in same ROI
        nr_per_roi_count = [];

        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            if ismember(nr, looped_nrs)
                continue;
            end
                
            % Grab the Roi label
            tokens = split(popul_data.neuron_name{nr}, "_");
            left_str = join(tokens(1:end - 1), "_");
            right_str = tokens(end);
            
            % Check that the neuron comes from a TICO recording, greater
            % than 1
            if str2num(right_str{1}) > 1
                % Grab neurons from the same ROI
                same_roi_nrs = find(contains(popul_data.neuron_name, left_str) == 1);
                looped_nrs = [looped_nrs, same_roi_nrs];
                nr_per_roi_count(end + 1) = length(same_roi_nrs);

            end
        end

        % Print out the TICO recordings number of neurons info
        disp([f_region ' ' f_stim]);
        disp(nr_per_roi_count);
        disp([num2str(mean(nr_per_roi_count)) ' ' num2str(std(nr_per_roi_count))]);
        
    end
end

%% -- (Deprecated) -- This is an old method
% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

ses = dir([pv_data_path '*.mat']);
all_matfiles = {ses.name};

% Loop through each brain region
[region_matfiles] = find_region(all_matfiles);
for f_region = fieldnames(region_matfiles)'
    f_region = f_region{1};

    % Create brain region struct
    region_stats.(f_region) = struct();
 
    % Grab all of the mice for this brain region
    [mouse_matfiles] = find_mouse(region_matfiles.(f_region).names);
    for f_mouse = fieldnames(mouse_matfiles)'
        f_mouse = f_mouse{1};
        
        % Initialize mouse stats
        mouse_stats.(f_mouse) = struct();
        
        % Select matfiles by stim condition by this iteration's brain region
        [matfile_stim] = find_stim_cond(mouse_matfiles.(f_mouse).names);
        
        % Loop through each field of the struct and concatenate everything together
        % Keep track the stats of the number of trials and neurons per mouse
        cond_stats = struct();
        % Loop through each stimulation condition
        for field = fieldnames(matfile_stim)'
            field = field{1};
            matfiles = matfile_stim.(field).names;    
         
            % Initialize field to keep track of neurons and trials
            cond_stats.(field) = struct();
            cond_stats.(field).trial_nums = [];
            cond_stats.(field).num_fovs = [];
            cond_stats.(field).fovs_rej = [];
       
            % Loop through each matfile of the current stimulation condition
            for matfile =matfiles
                % Read in the mat file of the current condition
                data = load([pv_data_path matfile{1}]);
                temp_idx = find(~cellfun(@isempty, data.align.trial));
                
                for roi_idx=1:size(data.align.trial{temp_idx(1)}.detrend_traces, 2)
                    
                    % OLD CODE
                    % Grab the number of trials to ignore from this neuron
                    %ri = strsplit(matfile{1}, '_');
                    %if isfield(ignore_trial_dict.(['mouse_' ri{1}]), (['rec_' erase(ri{3}, 'rec')])) && ...
                    %    isfield(ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]), (ri{4})) && ...
                    %    isfield(ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}), (['f_' ri{5}]))
        
                    %    num_ignore = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).ROI1;
                    %else
                    %    num_ignore = 0;
                    %end
        
                    %TODO I should just encapsulate this into a try/catch block and then catch with num_ignore = 0
                    ri = strsplit(matfile{1}, '_');
                    try 
                        num_ignore = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
                        
                        num_ignore = length(num_ignore);
                    catch
                        num_ignore = 0;
                    end
        
                    % Add the number of trials for this specific neuron and check if at least one trial is included in calculation
                    cond_stats.(field).trial_nums(end + 1) = length(data.align.trial) - num_ignore;
                    
                    % Count the neuron if it has greater than or equal to 3 trials
                    if length(data.align.trial) - num_ignore >= 3
                        cond_stats.(field).num_fovs(end + 1) = str2num(erase(ri{4}, 'FOV')); %[cond_stats.(field).num_fovs, str2num(erase(ri{4}, 'FOV')) ];
                        %erase(ri{4}, 'FOV')
                    else
                        cond_stats.(field).fovs_rej(end + 1) = str2num(erase(ri{4}, 'FOV')); %[cond_stats.(field).num_fovs, str2num(erase(ri{4}, 'FOV')) ];

                    end
                end % End looping through ROIs
            end % End looping through FOVs of a condition
        end % End looping through all of the mice 
        region_stats.(f_region).(f_mouse).cond_stats = cond_stats;
    end
end

% Loop through and show each regions conditions
regions = fieldnames(region_stats);
totals = struct();
for f_region = regions'
    f_region = f_region{1};
    disp(f_region);

    total_40_acp = 0;
    total_40_rej = 0;
    total_140_acp = 0;
    total_140_rej = 0;

    % Loop through each mouse
    mice = fieldnames(region_stats.(f_region));
    for f_mouse = mice'
        f_mouse = f_mouse{1};
        disp(f_mouse);

        stims = fieldnames(region_stats.(f_region).(f_mouse).cond_stats);
        cond_stats = region_stats.(f_region).(f_mouse).cond_stats;
        for stim=stims'
            stim = stim{1};
            disp(stim);
            % disp(['Average trial per neuron: ' num2str(nanmean(cond_stats.(stim).trial_nums)) ]);
            disp(['Num kept neurons: ' num2str(numel(cond_stats.(stim).num_fovs)) ]);
            fprintf(['FOVs: ' num2str(cond_stats.(stim).num_fovs) '\n\n']);
            disp(['Rejected neurons: ' num2str(numel(cond_stats.(stim).fovs_rej)) ]);
            fprintf(['FOVs: ' num2str(cond_stats.(stim).fovs_rej) '\n\n']);
        
            % Aggregate the number of neurons between 40Hz and 140Hz
            if strcmp(stim, 'f_40') == 1
                total_40_acp = total_40_acp + numel(cond_stats.(stim).num_fovs);
                total_40_rej = total_40_rej + numel(cond_stats.(stim).fovs_rej);
            elseif strcmp(stim, 'f_140') == 1
                total_140_acp = total_140_acp + numel(cond_stats.(stim).num_fovs);
                total_140_rej = total_140_rej + numel(cond_stats.(stim).fovs_rej);
            end
        end
        fprintf('\n\n');
    end

    totals.(f_region).acp_40 =  total_40_acp;
    totals.(f_region).rej_40 =  total_40_rej;
    totals.(f_region).acp_140 = total_140_acp;
    totals.(f_region).rej_140 = total_140_rej;
end

% Print all of the total neurons across mice
v1_total = totals.r_V1

m1_total = totals.r_M1

% Return matfiles by stimulation condition
function [cond_struct] = find_stim_cond(matfile_names)
    cond_struct = struct();
    
    % Loop through each matfilename and group by stimulation conditions
    for i=1:length(matfile_names)
            file_parts = split(matfile_names{i}, '_');
            stim = file_parts{5};
            
            % Create stimulation field if it does not exist
            if ~isfield(cond_struct, ['f_' stim])
                cond_struct.(['f_' stim]).names = {};
            end

            cond_struct.(['f_' stim]).names{end+1} = matfile_names{i};
    end
end


%% Specific functions for determining which FOVs to look at
function [region_struct] = find_region(matfile_names)
    region_struct = struct();

    for i=1:length(matfile_names)
        file_parts = split(matfile_names{i}, '_');
        reg = file_parts{2};
        
        % Create region field if it does not exist
        if ~isfield(region_struct, ['r_' reg])
            region_struct.(['r_' reg]).names = {};
        end

        region_struct.(['r_' reg]).names{end+1} = matfile_names{i};
    end
end

% Find matfiles by mouse cage number
function [mouse_struct] = find_mouse(matfile_names)
    mouse_struct = struct();

    for i=1:length(matfile_names)
        file_parts = split(matfile_names{i}, '_');
        mouse = file_parts{1};
        
        % Create mouse field if it does not exist
        if ~isfield(mouse_struct, ['m_' mouse])
            mouse_struct.(['m_' mouse]).names = {};
        end

        mouse_struct.(['m_' mouse]).names{end+1} = matfile_names{i};
    end
end
