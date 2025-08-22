data_path = ['~/handata_server' f 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Get all FOV matfiles
matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);
    
    stim = strsplit(matfile_names{i}, '_');
    stim = stim{5};
        
    for tr = 1:length(data.raw.trial)
        rawdata_trial = data.raw.trial{tr};
        aligndata_trial = data.align.trial{tr};

        % Skip if empty
        if isempty(rawdata_trial)
            continue;
        end


        % Check if there are extra pulses
        if str2num(stim) ~= length(rawdata_trial.raw_stimulation_time) %|| rawdata_trial.raw_stimulation_time(end) > aligndata_trial.camera_frame_time(end)
            %figure('visible', 'off');
            %plot(rawdata_trial.raw_stimulation_time, 2, '|r');
            %hold on;
            %plot(aligndata_trial.camera_frame_time, 1, '|b');
            %ylim([0, 3]);
            %
            %% Show the incorrect files
            if abs(length(rawdata_trial.raw_stimulation_time) - str2num(stim) ) > 1
                num_stim_pulses = length(rawdata_trial.raw_stimulation_time);
                fprintf(['\n\n' matfile_names{i} ' ' num2str(tr) '\n']);
                fprintf(['Num pulses ' num2str(num_stim_pulses) '\n']);
                pause;
            end

            fprintf( num2str(length(rawdata_trial.raw_stimulation_time) - str2num(stim) ) );

        end
    end
end
