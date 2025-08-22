clc;
close all;
clear all;

%% Read in the matfiles
f_sep = filesep;
addpath('..');

interm_data_path = [f_sep 'home' f_sep 'pierfier' f_sep 'Projects' f_sep ...
      'Pierre Fabris' f_sep 'PV DBS neocortex' f_sep 'Interm_Data' f_sep];

server_root = '/home/pierfier/handata_server/';

data_root_path = [server_root 'eng_research_handata3' f_sep 'Yangyang_Wang' f_sep 'PV_V1_LED_SomArchon' f_sep];
matfiles = dir(fullfile(data_root_path, '**', '*flicker*DBS*fov*.mat'));

% Remove matfiles in TODO
indices = ~contains({matfiles.folder}, 'TODO');
matfiles = matfiles(indices);

% Openephys dual scope path
addpath([server_root 'Cara_Ravasio' f_sep 'Code' f_sep 'OpenEphys_Scripts']);

% Grab the ignore file csv
ignore_dict = Multi_func.csv_to_struct('flicker_ignore.csv');

%% Parameters for recordings
total_frames = 2500;
total_time = 5.0258; % in Sec
soft_Fs = total_frames/total_time; % soft here refers to software
theo_Fs = 500; % theo refers to theoretical

front_frame_drop = 14;

%% Test the alignment for a single matfile
start_chan = 37;
frame_chan = 33;
flick_chan = 35;
stim_chan = 36;

% Find matfiles with 40 Hz
% originally for debugging mat_40s = find(contains({matfiles.name}, 'DBS40hz'));

% Store all of the data here
data = struct;

% Store neuron recording that had more than 10 trials in, for debugging purposes
more_than_10 = {};

for i= 1:length(matfiles) % 34:36 %
    i
    matfile = fullfile(matfiles(i).folder, matfiles(i).name)
    
    % Grab all of the parameters for experiment
    try
        cur_flicker = Multi_func.get_param(matfiles(i).name, 'flicker');
        cur_freq = Multi_func.get_param(matfiles(i).name, 'DBS');
        cur_fov = Multi_func.get_param(matfiles(i).name, 'fov');
    catch ME
        continue;
    end

    % Create frequency condition field
    f_stim = ['f_' erase(cur_freq, 'DBS')];
    f_stim = f_stim(1:end - 2);
    
    % Create neuron number field
    f_nr = ['f_' num2str(i)];

    % Example neuron name: 109558_Vb_male_20240308_2
    mat_parts = split(matfiles(i).name, '_');
    name_parts = {mat_parts{2}, mat_parts{3}, cur_fov};
    nr_name = strjoin(name_parts, '_');    

    % Find ephys folder
    ephys = dir([matfiles(i).folder f_sep 'ephys' f_sep '*' cur_flicker '*' cur_freq '*' cur_fov '*']);
    
    % Only use new ephys file if it exists
    if length(ephys) > 0
        ephys_path = fullfile(ephys.folder, ephys.name);
        ephys_path = dir([ephys_path f_sep '**/*structure*']);
        ephys_path = fullfile(ephys_path.folder, ephys_path.name);
    
    % Read in the model ephys folder
    else
        ephys = dir(['.' f_sep '*' cur_flicker '*' cur_freq]);
        ephys_path = fullfile(ephys.folder, ephys.name);
        ephys_path = dir([ephys_path f_sep '**/*structure*']);
        ephys_path = fullfile(ephys_path.folder, ephys_path.name);
    end

    % Read ephys data
    D = load_open_ephys_binary(ephys_path, 'continuous', 1);
    
    % Load traces file
    trace_data = load(matfile);
    
    try
        % Grab the trial start times
        all_start_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(start_chan, :)');
        all_stim_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(stim_chan, :)');
        all_frame_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(frame_chan, :)');
        all_flick_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(flick_chan, :)');
        all_flickoff_times = flip(Multi_func.get_ephys_rise_times(flip(D.Timestamps), flip(D.Data(flick_chan, :)') ) );
        
    catch ME
        % Use the original ADC channel ports
        % Grab the trial start times
        all_start_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(5, :)');
        all_stim_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(4, :)');
        all_frame_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(1, :)');
        all_flick_times = Multi_func.get_ephys_rise_times(D.Timestamps, D.Data(3, :)');
        all_flickoff_times = flip(Multi_func.get_ephys_rise_times(flip(D.Timestamps), flip(D.Data(3, :)') ) );
    end
    
        % debugging purposes, found ephys file with more than 10 trials!!
        if length(all_start_times) > 10 % & ~contains(ephys_path, 'mod_')
            disp('foujnd');
            more_than_10{end + 1} = {ephys_path, length(all_start_times)};
        end

    % Grab the trials that need to be ignored
    try
        mouse_fn = fieldnames(ignore_dict);
        mouse_idx = find(contains(mouse_fn, name_parts{1})); % Use the name split
        mouse_f = mouse_fn(mouse_idx);
        fov_num = erase(cur_fov, 'fov');
        ign_trials = ignore_dict.(mouse_f{1}).(['rec_' num2str(name_parts{2})]).(['FOV' fov_num]).(f_stim).ROI1;
    catch ME
        ign_trials = [];
    end

    % Loop through each trial to align
    raw_vm = [];
    sp_amp_raster = [];
    camera_frame_times = [];
    stim_times = [];
    flicker_times = [];
    flicker_off_times = [];

    %TODO I am thinking I should add a condition if there are more trials than the ephys has. Then I can just concatenate extra trials to the ephys data
    
    % Repeat the last trial start time for neurons that have more trials than the ephys amount
    diff_trials = max(unique(trace_data.roi_list.trial_vec)) - length(all_start_times);
    all_start_times(end + 1:end + diff_trials) = all_start_times(end);
    % Quick check
    %if diff_trials > 1
    %error('Pause here');
    %end

    for tr_i = unique(trace_data.roi_list.trial_vec)
        
        % Skip trial if it is in the trials array
        if ismember(tr_i, ign_trials)
            continue;
        end

        trace = trace_data.roi_list.traces(tr_i == trace_data.roi_list.trial_vec);
        trace = trace(front_frame_drop + 1:end);
        
        % Detrend the trace
        trace_mov = movmean(trace, theo_Fs);
        trace = (trace - trace_mov)./trace_mov;
        
        % Grab the time frames
        tr_all_frame_times = all_frame_times(all_frame_times > all_start_times(tr_i));
        tr_all_frame_times = tr_all_frame_times(front_frame_drop + 1:end);
        tr_all_frame_times = tr_all_frame_times(1:length(trace));

        % Grab the flicker times
        tr_flicker_times = all_flick_times(all_flick_times > tr_all_frame_times(1) & all_flick_times < tr_all_frame_times(end) );

        tr_flickeroffset_times = all_flickoff_times(all_flickoff_times > tr_all_frame_times(1) & all_flickoff_times < tr_all_frame_times(end) );
 
        % Grab stim times
        tr_all_stim_times = all_stim_times(all_stim_times > tr_all_frame_times(1) & all_stim_times < tr_all_frame_times(end) );
   
        % Determine the spike raster
        sp_raster = nan(size(trace_data.spike_detect_SA_v4_info{tr_i}.roaster));
        sp_raster(trace_data.spike_detect_SA_v4_info{tr_i}.spike_idx{1}) = trace_data.spike_detect_SA_v4_info{tr_i}.spike_amplitude{1};
        sp_raster = sp_raster(front_frame_drop:end);
    
        % Store the spike amplitude information
        sp_amp_raster = [sp_amp_raster, sp_raster(:)];
        
        % Store the Vm
        raw_vm = [raw_vm, trace(:)];
    
        % Subtract the times by the flicker onset time
        tr_all_frame_times = tr_all_frame_times - tr_flicker_times(1);
        tr_all_stim_times = tr_all_stim_times - tr_flicker_times(1);
        tr_flickeroffset_times = tr_flickeroffset_times - tr_flicker_times(1);
        tr_flicker_times = tr_flicker_times - tr_flicker_times(1);
            
        % Store the times
        camera_frame_times = [camera_frame_times, tr_all_frame_times(:)];
        stim_times = [stim_times, tr_all_stim_times(:)];
        flicker_times = [flicker_times, tr_flicker_times(:)];
        flicker_off_times = [flicker_off_times, tr_flickeroffset_times(:)];

        % DEBUG
        %figure;
        %plot(tr_all_frame_times, ones(size(tr_all_frame_times)), '|');
        %hold on;
        %plot(tr_all_stim_times, 1 + ones(size(tr_all_stim_times)), '|');
        %hold on;
        %plot(tr_flicker_times, 2 + ones(size(tr_flicker_times)), '|');
        %hold on;
        %plot(tr_flickeroffset_times, 3 + ones(size(tr_flickeroffset_times)), '|');
        %hold on;

        % DEBUG
        %figure;
        %plot(trace);
        %hold on;
        %plot(find(~isnan(sp_raster)), trace(find(~isnan(sp_raster))), 'or');
        
        % Print the difference
        %DEBUG
        %tr_flicker_times(1) - tr_all_frame_times(1)
    end

    % If all trials were ignored, do not save the current neuron to the population data
    if isempty(raw_vm)
        continue
    end

    %-- Save all of the properties to data struct
    
    % Check if neuron data exists
    if ~isfield(data, f_stim)
        data.(f_stim) = struct;
        data.(f_stim).raw_vm = struct;
        data.(f_stim).sp_amp_raster = struct;
        data.(f_stim).camera_frame_times = struct;
        data.(f_stim).stim_times = struct;
        data.(f_stim).flicker_times = struct;
        data.(f_stim).flicker_off_times = struct;
        data.(f_stim).nr_name = struct;

    end

    % Save the neuron data
    data.(f_stim).raw_vm.(f_nr) = raw_vm;
    data.(f_stim).sp_amp_raster.(f_nr) = sp_amp_raster;
    data.(f_stim).camera_frame_times.(f_nr) = camera_frame_times;
    data.(f_stim).stim_times.(f_nr) = stim_times;
    data.(f_stim).flicker_off_times.(f_nr) = flicker_off_times;
    data.(f_stim).flicker_times.(f_nr) = flicker_times;
    data.(f_stim).nr_name.(f_nr) = nr_name;
end

%% Save to matfile
% Convert struct to individual fields
fields = fieldnames(data);
for i = 1:numel(fields)
    assignin('caller', fields{i}, data.(fields{i}));
end

% Save the data into a matfile
save([interm_data_path 'v1_flicker_aligned_dat.mat'], fields{:});

%% Check all of the continuous data
ephys_data = Multi_func.norm_signals(D.Data(33:37, :)');

figure;
plot(D.Timestamps, ephys_data + [1:size(ephys_data, 2)]);
title(num2str(i));


%% 

ignore_trials = {
    '109558_Vb_male':{
        '20240311':{
            1:{40:[3, 4, 5, 8, 10],
               140:[1, 2, 8, 9]},
            2:{40:[4, 5, 7, 8, 9, 10],
               140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
            3:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
            4:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
         },
         '20240308':{
             1:{140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
             2:{140:[2, 5, 6, 8, 9, 10, 11, 14, 16, 17, 18, 19]},
             3:{140:[6, 7, 14, 15, 16]}
         }
    },
    '109567_Vb_male':{
        '20240311':{
            1:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11],
               140:[4, 8, 9]},
            2:{40:[3, 5, 6, 7, 8],
               140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]},
            3:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               140:[]},
            4:{40:[2, 3, 6, 7, 8, 9, 10],
               140:[1, 2, 3, 4, 9, 10]},
           5:{40:[7, 8, 9],
               140:[1, 2, 3, 10]},
            6:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               140:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
        },
        '20240411':{
            1:{40:[6, 7, 8, 10],
               140:[2, 3, 5, 8, 9, 10]},
           2:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
               140:[2, 3, 4, 5, 8, 9, 10]},
            3:{40:[1, 6, 7, 8, 9],
               140:[1, 3, 5]},
            4:{40:[4, 5, 7, 9],
               140:[7, 9, 10]},
            5:{40:[1, 6, 7],
               140:[1, 3, 4, 5, 6, 7]}
        },
        '20240424':{
            2:{40:[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]}
        }
    }
}
