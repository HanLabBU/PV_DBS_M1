clc;
clear all;
close all;

%%
f= filesep;

addpath('..');

% This structure is taken from the python export file
interm_data_path = [f 'home' f 'pierfier' f 'Projects' f ...
      'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f];

savefig_path = [Multi_func.save_plot  'Flicker' f];


%% Read in both matfiles

save_data_file = [interm_data_path 'v1_flicker.mat'];
data_og = load(save_data_file);

save_data_file = [interm_data_path 'v1_flicker_aligned_dat.mat'];
data_align = load(save_data_file);

%% Compare files between each aligned matfile
% Loop through each stim frequency
for f_stim = fieldnames(data_og)'
    f_stim = f_stim{1};
    
    % Note: popu_data is from data_og
    popul_data = data_og.(f_stim);
    
    % Loop through the neurons of the original data
    for f_nr_og = fieldnames(popul_data.raw_vm)'
        f_nr_og = f_nr_og{1};
        

        og_nr_name = popul_data.nr_name.(f_nr_og);
        og_parts = split(og_nr_name, '_');

        % Ensure that neuron was found
        nr_found = 0;

        % Loop through the neurons of the aligned data
        for f_nr_al = fieldnames(data_align.(f_stim).raw_vm)'
            f_nr_al = f_nr_al{1};
        

            al_nr_name = data_align.(f_stim).nr_name.(f_nr_al);

            %TODO need to compare the individual points on 
            al_parts = split(al_nr_name, '_');

            % Check if the neurons match
            if strcmp(og_parts{1}, al_parts{1}) & strcmp(og_parts{4}, al_parts{2}) & ...
                strcmp(og_parts{5}, erase(al_parts{3}, 'fov'))
                
                nr_found = 1;

                if ~isequal(size(popul_data.raw_vm.(f_nr_og)), ...
                        size(data_align.(f_stim).raw_vm.(f_nr_al) ) )

                    disp('OG size: ');
                    size(popul_data.raw_vm.(f_nr_og))
                    
                    disp('Al size: ');
                    size(data_align.(f_stim).raw_vm.(f_nr_al) )
    

                    error('Neurons do not have the same trace size');
                
                else
                    
                    figure('Position', [5, 5, 2000, 1000]);
                    %tiledlayout(2, 1);
                    %nexttile;
                    data = Multi_func.norm_signals(popul_data.raw_vm.(f_nr_og));
                    plot(data + [1:size(popul_data.raw_vm.(f_nr_og), 2)], 'r', 'LineWidth', 4);
                    %nexttile;
                    hold on;
                    data = Multi_func.norm_signals(data_align.(f_stim).raw_vm.(f_nr_al));
                    plot(data + [1:size(data_align.(f_stim).raw_vm.(f_nr_al), 2)], 'b');

                    % Move to next neuron
                    break;
                end
            end
        end

        if ~nr_found
            disp('Not found');
        end
    end
end
