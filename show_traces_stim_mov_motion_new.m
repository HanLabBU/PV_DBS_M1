f = filesep;

%USER modification
% Server root path
server_rootpath = 'Z:\';

% Maingear office computer
%server_rootpath = '/home/pierfier/handata_server/';

% Folder location of the saved and aligned data
data_path = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];

%END modification

% Get all FOV matfiles
matfile_names = dir([data_path '*.mat']);
matfile_names = {matfile_names.name};

% Loop through each matfile
for i=1:length(matfile_names)
    data = load([data_path matfile_names{i}]);

    % Setup figure to show alignment data for all trials
    figure('Position', [0, 0, 1000, 5000]);
    tiledlayout(length(unique(result.trial_vec)), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
end

% Perform exponential fit
function [x, y]  = exp_fit(trace, Fs)
    t = 1:length(trace);
    f2 = fit(t', trace, 'exp2');
    y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
    x = t;
end
