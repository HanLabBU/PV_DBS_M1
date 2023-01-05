f = filesep;

%USER modify everything below here
%% Specific to the computer being used
%DMD
server_rootpath = 'Z:\';

%Dual scope
%server_rootpath = '\\engnas.bu.edu\research\eng_research_handata\';

% Maingear office computer
%server_rootpath = '/home/pierfier/handata_server/';

% Folder location of all the saved and aligned data
data_path = [server_rootpath 'EricLowet' f 'DBS' f 'PV' f];

%USER END modification

frame_exp = 1.2;
Fs = 1/(1.2*10^-3);

front_frame_drop = 15;
back_frame_drop = 2500;

% Scripts folder for performing data alignment and saving relevant data to mat files
addpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'LFP_Trace_Alignment' ]);
addpath(genpath([server_rootpath 'EricLowet' f 'Scripts']));

matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each FOV matfile
for i=1:length(matfile_names)
    %TODO should change this to a variable such as 'data' to prevent variable remenance into subsequent loops
    load([data_path matfile_names{i}]);
    
    figure('Position', [0, 0, 1000, 5000]);
    tiledlayout(length(unique(result.trial_vec)), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    % Loop through each trial
    for  ne=unique(result.trial_vec)
        
        v= result.traces(result.trial_vec==ne); %eb edited 20211111
        
        % Including just the frames without any artefacts
        v = v(front_frame_drop:back_frame_drop);
        [x, y] = exp_fit(v(:), Fs);
        detrend_trace = v - y';
        %   v= result.traces(result.trial_vec==ne,3);%./result.tracesB(result.trial_vec==ne)   ;
        v=v-fastsmooth(v,1400,1,1); 
        nexttile;
        plot([1:length(v)]/Fs, detrend_trace, '-b');
        
        
    end
    
    % Set the title of the Figure
    sgtitle(matfile_names{i}, 'Interpreter', 'none');
end

% Perform exponential fit
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
