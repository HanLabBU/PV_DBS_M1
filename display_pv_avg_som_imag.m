
% Eric's room camera
camera_func = @read_dcimg_kk_scope_new;
% Eric's room other camera
camera_func = @read_dcimg_kk_scope;
% Treadmill room camera
%camera_func = @read_dcimg_treadmill_cam;


%% Average just single trial
close all;
% Path to dcimg recording file
%trial_file = '~/Projects/Pierre Fabris/PV DBS neocortex/Stim Recordings/31556noeartag_M1/20221206/';
trial_file = '~/Projects/Pierre Fabris/PV DBS neocortex/Stim Recordings/50464_M1/20230711/';

fov = 'FOV2'; 
stim = '_140_';

ses = dir([trial_file fov '*' stim '*.dcimg']);
all_trial_stack = [];
for file={ses.name}
    file = file{1};
    
    % Loop through all trials into one stack file
    stack = camera_func([trial_file file]);
    all_trial_stack(:, :, :, end + 1) = stack;

    display_cells_in_frame(mean(stack, 3), {}, [10, 10, 28, 4]);
end
all_trial_stack(:, :, :, 1) = [];


%% Average all trials for given FOV with stimulation frequency

avg_tr_frames = squeeze(mean(all_trial_stack, 3));

%%
% Show the average frame with a scale bar

frame = mean(avg_tr_frames, 3);

display_cells_in_frame(frame, {}, [40, 160, 28, 4]);
