%% Read in COMSOL model data
clear all;
close all;
clc;

f = filesep;

% Dropbox electrode data path
dropbox_path =  ['~/Dropbox/Pierre PV DBS Project Dropbox/Materials/'];


% Read in electric field with varying depth COMSOL data
comsol_model_t = readtable([dropbox_path 'Electrode COMSOL Data' f 'electrodeDepths.csv']);
%% Plot electric field profile with different electrode depths

figure('Renderer','painters');

% Loop through each variable name pairing
vars = comsol_model_t.Properties.VariableNames;

for i = 1:2:numel(vars)
    depthName = vars{i}
    
    % Parse out depth
    depth = extractBetween(depthName, "x", "_");
    depth = str2num(depth{1})
    
    plot(comsol_model_t{:, i}, comsol_model_t{:, i+1}, 'DisplayName', [num2str(depth) ' µm']);
    hold on;
end
legend();
xlabel('Distance from one electrode (m)');
ylabel('Electric field strength (V/m)');

Multi_func.set_default_axis(gca);
ax = gca;
ax.Units = 'centimeters';
ax.InnerPosition = [2 4 4.9 3.2];

fontsize(gcf, 10, 'points');

saveas(gcf, [dropbox_path 'Plots'  f 'Comsol' f 'Electrode_Depth_Varied.png']);
saveas(gcf, [dropbox_path 'Plots'  f 'Comsol' f 'Electrode_Depth_Varied.pdf']);