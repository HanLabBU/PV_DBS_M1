f = filesep;

%USER modify everything below here
%% Specific to the computer being used
%DMD
%server_rootpath = 'Z:\';

%Dual scope
%server_rootpath = '\\engnas.bu.edu\research\eng_research_handata\';

% Maingear office computer
server_rootpath = '/home/pierfier/handata_server/';

% Folder location of all the saved and aligned data
data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

frame_exp = 1.2;
frame_freq = 1/(1.2*10^-3);

%USER END modification

% Scripts folder for performing data alignment and saving relevant data to mat files
addpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'LFP_Trace_Alignment' ]);
addpath(genpath([server_rootpath 'Pierre Fabris' f 'DMD Project' f 'DMD Scripts' f 'In Vitro Analysis Scripts']));
addpath(genpath([server_rootpath 'EricLowet' f 'Scripts']));

matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each FOV matfile
for i=1:length(matfile_names)
    load([data_path matfile_names{i}]);
    
    %Initialize struct to save all of the spike detection data
    result.resultS = [];

    % Loop through each trial
    for  ne=unique(result.trial_vec)
        v= result.traces(result.trial_vec==ne, :); %eb edited 20211111
        %Dropping the first 50 frames, there appears to be artefacts
        %v(1:50, :) = NaN;
        
        % Use DMD spike detection script
        result.resultS{end+1} = spike_detect_SNR_v3b(v);
    end

    % Update the matfile
    save([data_path matfile_names{i}], 'result', '-append');
    
    % Set the title of the Figure
    %sgtitle(matfile_names{i}, 'Interpreter', 'none');
end
