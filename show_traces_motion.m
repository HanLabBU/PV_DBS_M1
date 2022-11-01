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
data_path = [server_rootpath 'EricLowet' f 'DBS' f 'PV' f];

frame_exp = 1.2;
frame_freq = 1/(1.2*10^-3);

%USER END modification

% Scripts folder for performing data alignment and saving relevant data to mat files
addpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'LFP_Trace_Alignment' ]);
addpath(genpath([server_rootpath 'EricLowet' f 'Scripts']));

matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each FOV matfile
for i=1:length(matfile_names)
    load([data_path matfile_names{i}]);
    
    figure('Position', [0, 0, 1000, 5000]);
    tiledlayout(length(unique(result.trial_vec)), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each trial
    for  ne=unique(result.trial_vec)
        
        %DEBUG
        traces_sice = size(result.traces);
    

        v= result.traces(result.trial_vec==ne); %eb edited 20211111
        %   v= result.traces(result.trial_vec==ne,3);%./result.tracesB(result.trial_vec==ne)   ;
        v=v-fastsmooth(v,1400,1,1); 
        %Dropping the first 50 frames, there appears to be artefacts
        v(1:50) = NaN;
        nexttile;
        plot([1:length(v)]/frame_freq, v);
    end
    
    % Set the title of the Figure
    sgtitle(matfile_names{i}, 'Interpreter', 'none');
end
