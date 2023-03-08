f = filesep;

%USER modify everything below here
%% Specific to the computer being used
%DMD
%server_rootpath = 'Z:\';

%Dual scope
%server_rootpath = '\\engnas.bu.edu\research\eng_research_handata\';

% Maingear office computer
server_rootpath = '~/Projects/';

% Folder location of all the saved and aligned data
data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

%USER END modification

matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each FOV matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);
    
    % Loop through each trial
    for  tr=1:length(data.align.trial)
        if ~isempty(data.align.trial{tr})
            % Use DMD spike detection script
            [spike_info] = spike_detect_SNR_v3b(data.raw.trial{tr}.raw_traces, 3.75);
            data.align.trial{tr}.spike_info375 = spike_info;
        end
    end
    
    % Update the matfile
    align = data.align;
    save([data_path matfile_names{i}], 'align', '-append');
    
    % Set the title of the Figure
    %sgtitle(matfile_names{i}, 'Interpreter', 'none');
end
