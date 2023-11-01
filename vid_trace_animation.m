% Housekeeping clear all stuff
clc;
clear all;
close all;
f = filesep;

% Maingear office computer
local_root_path = '/home/pierfier/Projects/';
server_root_path = '~/handata_server/';

% Path to PV data matfiles
pv_data_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

% Grab matfile with trace
example_matfile = [pv_data_path '23072_V1_rec20220217_FOV3_40_220_.mat'];
tr_idx = 3; 

%Video file
mc_vid_file = ['/home/pierfier/Projects/Pierre Fabris/PV DBS neocortex/Stim Recordings/611284_V1/20210827/motion_corrected/' 'motcorrected_fov1_140_60_00005_denoised.tif'];

% Check the file format 
if strcmp(mc_vid_file(end - 3:end), '.tif') == 1
    frames = loadtiff(mc_vid_file);
elseif strcmp(mc_vid_file(end - 3:end), '.avi') == 1
    % Read in motion corrected video
    v = VideoReader(mc_vid_file);
    frames = read(v);
    frames = frames(:, :, 1, :);
    frames = squeeze(frames);
end


% Read in mat trace file and grab trial data
data = load(example_matfile);
trial_data = data.align.trial{tr_idx};
trace = trial_data.detrend_traces;

figure();
tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:size(frames, 3)
    nexttile(1);
    %imad = imadjust(frames(:, :, i)); %, [0.5 0.8]);
    imagesc(frames(:, :, i));
    colormap(gray);

    nexttile(2);
    plot(trace, 'k')
    xline(i, 'b');
    pause(.010);
end
