% Plot the currents used in the experiments
clear all;
close all;
clc;

f = filesep;

% Linux server
local_root_path = '~/Projects/';
% Handata Server on Linux
server_root_path = '~/handata_server/eng_research_handata3/';
% Windows server
%local_root_path = 'Z:\';

% List path where all of the matfiles are stored
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on local linux machine
%pv_data_path = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'PV_Data' f];
% Data on handata3 folder
pv_data_path = [server_root_path 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
savepath = Multi_func.save_plot;

exclude_200ms = 1;

% CSV file to determine which trials to ignore
ignore_trial_dict = Multi_func.csv_to_struct([local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f ...
                                       'Stim Recordings' f 'Data_Config' f 'byvis_ignore.csv']);

% Check that the server path exists
if ~isfolder(local_root_path)
    disp('Server rootpath does not exist!!!');
    return;
end

%% Load all region data

% Use congregated data to save all of the currents
% Read in the saved pv data and perform analysis
if ~exclude_200ms
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data.mat'];
else
    save_all_data_file = [local_root_path 'Pierre Fabris' f 'PV DBS neocortex' f 'Interm_Data' f 'pv_data_ex200.mat'];
end
tic;
%Load the data
load(save_all_data_file);
toc
disp('Finished Loading Data');

%% TODO <--- determine if the file that has the current info for the flicker experiments is older than the interm pickle file


% TODO Also read in that data here somehow
%% Loop through and save the current amplitude for each neuron
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);
    
    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Store the current amplitude for each neuron
        pop_amps = [];        
        
        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            tokens = regexp(popul_data.neuron_name{nr}, '_', 'split');
            
            % Parse out the current amperage
            fov_i = find(contains(tokens, 'FOV') == 1);
            amp_i = fov_i + 2;
            if contains(tokens{amp_i}, '.') == 1
                amp_str = regexp(tokens{amp_i}, '\.', 'split');
                amp_str = amp_str{1};
            else
                amp_str = tokens{amp_i};
            end
            
            % Store the amperage into array
            amp = str2num(amp_str);
            pop_amps(end + 1) = amp;

        end

        % Save the currents
        region_data.(f_region).(f_stim).currents = pop_amps;
    end
end

%% Read in the flicker amplitude data into our currents structure
% Note: Need to only do this once, otherwise neuron names will keep getting appended
% to the end of region_data

% Check if there are more than XX number of neurons for V1, if there are,
% it means that the flicker experiment data was already added
if ~(length(region_data.r_V1.f_40.neuron_name) > 27)
    T = readtable([savepath 'Current' f 'flicker_current.csv']);
    % Convert the recording session date to string
    T.session = string(T.session);
    
    for i=1:height(T)
        freq = T{i, 'stim_freq'};
        params = T{i, 'nr_params'}{:};
    
        amp_str = regexp(params, '\d+', 'match');
        amp = str2num(amp_str{:});
    
        % Add current data back into structure
        region_data.r_V1.(['f_' num2str(freq)]).currents ... 
            = [region_data.r_V1.(['f_' num2str(freq)]).currents, amp];
    
        % Parse out the mouse name and recording session
        mouse_name = T{i, 'mouse'};
        mouse_name = mouse_name{:};
        mouse_name = strsplit(mouse_name, '_');
        mouse_name = mouse_name{1};
       
        session_id = T{i, 'session'};
        
        % Construct full neuron name similar to the regular 3-second stimulation trials
        nr_name_array = string({mouse_name, "V1", session_id, "FOV", freq, amp});
        full_nr_name = join(nr_name_array, '_');
        
        % Add neuron name back into region_data structure
        region_data.r_V1.(['f_' num2str(freq)]).neuron_name{end + 1} = char(full_nr_name);
    end
end

%% Print the amplitude ranges for each brain region and stimulation frequency
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    stims = fieldnames(region_data.(f_region));

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);

        % Display the conditions
        disp([f_region ' ' f_stim]);
        disp(['Current Range: ' num2str(min(popul_data.currents)) ...
                               '-' num2str(max(popul_data.currents)) ]);        
        disp([num2str(mean(popul_data.currents)) '+-' num2str(std(popul_data.currents)) ]);        
        fprintf('\n');
    end
end

%% Plot the current amplitudes for each brain region and stimulation frequency
% Also do the current statistical test for current amplitude

% Store the data and labels for the violinplots
data = [];
labels = [];
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Add the current amplitudes to violin arrays
        data = [data, popul_data.currents];
        labels = [labels, repmat({[f_region ' ' f_stim]}, 1, length(popul_data.currents))];
    end
end

% Make the violinplots
figure('Renderer', 'Painters');
violinOpts = Multi_func.get_default_violin();
violinplot(data(:)', labels, 'GroupOrder', ...
    {'r_M1 f_140', 'r_M1 f_40', 'r_V1 f_140', 'r_V1 f_40'}, ...
    violinOpts);
ylabel('Current (uamp)');
Multi_func.set_default_axis(gca);

saveas(gcf, [savepath 'Current/' 'violin_current_by_condtion.png']);
saveas(gcf, [savepath 'Current/' 'violin_current_by_condtion.pdf']);

% Perform statistical analysis on the current amplitudes and log the stats
stats_log = [savepath 'Current' f 'region_stim_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

diary on;
[p, h, stats] = kruskalwallis(data, labels)
disp('Group Columns');
disp(stats.gnames');
c = multcompare(stats, 'CriticalValueType', 'dunn-sidak')
diary off;

%% Loop through and calculate charge density for each neuron
% Perform statistical analysis on the current amplitudes and log the stats
stats_log = [savepath 'Current' f 'charge_density_calcs'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off;


electrode_diam = 127*10^-4; % 127 um converted to cm
electrode_shave_len = 500 * 10^-4; % 500 um converted to cm

total_electrode_area = (pi * ((electrode_diam/2)^2)) + ...
    (pi * electrode_shave_len * electrode_diam);

phase_width = 200*10^-6; % S

data = [];
labels = [];
for f_region = fieldnames(region_data)'
    f_region = f_region{1};

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Charge densities
        charg_dens = popul_data.currents * (phase_width)/(total_electrode_area)

        data = [data, charg_dens];
        labels = [labels, repmat({[f_region ' ' f_stim]}, 1, length(charg_dens))];

        region_data.(f_region).(f_stim).charge_density = charg_dens;

        % Plot the average and standard deviation of the charge density used
        diary on;
        disp([f_region ' ' f_stim]);
        disp([num2str(mean(charg_dens)) '+- ' num2str(std(charg_dens)) ' (uamp/cm^2)']);
        fprintf('\n\n');
        diary off;
    end
end

% Make the violinplots
figure('Renderer', 'Painters');
violinOpts = Multi_func.get_default_violin();
violinplot(data(:)', labels, 'GroupOrder', ...
    {'r_M1 f_140', 'r_M1 f_40', 'r_V1 f_140', 'r_V1 f_40'}, ...
    violinOpts);
ylabel('Current (uamp)');
Multi_func.set_default_axis(gca);

saveas(gcf, [savepath 'Current/' 'violin_charge_density_by_condtion.png']);
saveas(gcf, [savepath 'Current/' 'violin_charge_density_by_condtion.pdf']);

%% Plotting current from days to surgery for each neuron
% Red: 140Hz
% Blue: 40Hz
% Circle: M1
% Plus: V1

%A map that keeps track of how many repetitions there are
%of number of days since surgery and amperage
days_amp_map = [];
mouse_color = struct();

figure('Position', [0 0 1000 800]);

% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Keep track of neurons plotted
        neurons_plotted = [];

        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            if ismember(nr, neurons_plotted)
                continue;
            end

            tokens = regexp(popul_data.neuron_name{nr}, '_', 'split');

            % Find all of the neurons that come from the same mouse
            same_m_nrs = find(contains(popul_data.neuron_name, tokens{1}) == 1);
            same_m_nrs_pts = [];
            neurons_plotted = cat(2, neurons_plotted, same_m_nrs);
            for i=same_m_nrs

                tokens = regexp(popul_data.neuron_name{i}, '_', 'split');

                % Calculate the days between imaging and surgery date
                surg_date = popul_data.surgery_date.(['m_' tokens{1}]);
                if isempty(surg_date)
                    continue;
                end
                tokens{3} = erase(tokens{3}, 'rec');

                surg_date = datetime(surg_date, 'InputFormat', 'yyyyMMdd');
                img_date = datetime(tokens{3}, 'InputFormat', 'yyyyMMdd');

                days_from_surg = days(img_date - surg_date);
            
                %DEBUG
                if days_from_surg > 400
                    disp(popul_data.neuron_name{i});
                end
                    
                % Parse out the current amperage
                fov_i = find(contains(tokens, 'FOV') == 1);
                amp_i = fov_i + 2;
                if contains(tokens{amp_i}, '.') == 1
                    amp_str = regexp(tokens{amp_i}, '\.', 'split');
                    amp_str = amp_str{1};
                else
                    amp_str = tokens{amp_i};
                end

                % Set labels based on brain region and stim frequency
                plot_str = '';
                if strcmp(f_region, 'r_M1') == 1
                    plot_str(end+1) = 'o';
                elseif strcmp(f_region, 'r_V1') == 1
                    plot_str(end+1) = '+';
                end

                % Set the frequency stuff
                if strcmp(f_stim, 'f_40') == 1
                    plot_str(end+1) = 'b';
                elseif strcmp(f_stim, 'f_140') == 1
                    plot_str(end+1) = 'r';
                end

                % Increment the marker size for repeated points
                if size(days_amp_map, 1) < days_from_surg | size(days_amp_map, 2) < str2num(amp_str)
                    days_amp_map(days_from_surg, str2num(amp_str)) = 1;
                else
                    days_amp_map(days_from_surg, str2num(amp_str)) = ...
                        1 + days_amp_map(days_from_surg, str2num(amp_str));
                end

                % Plot data point
                hold on;
                plot(days_from_surg, str2num(amp_str), plot_str, 'MarkerSize', 8 + days_amp_map(days_from_surg, str2num(amp_str)));
                hold on;

                % Keep track of points
                same_m_nrs_pts(:, end + 1) = [days_from_surg; str2num(amp_str)];
            end

            % Plot lines across all FOVs of the same mouse
            %TODO sort the x values from decreasing to increasing, makes figure easier to read
            if ~isempty(same_m_nrs_pts)
                
                [~, i] = sort(same_m_nrs_pts(1, :));
                if isfield(mouse_color, ['m_' tokens{1}])
                    color = mouse_color.(['m_' tokens{1}]);
                    plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i), '-', 'Color', color);
                else
                    plot_h = plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i),  '-');
                    mouse_color.(['m_' tokens{1}]) = get(plot_h, 'Color');
                end
                hold on;
            end
        end
    end
end


% Add legend label for each item
text(400, 600, '+: V1 region');
hold on;
text(400, 580, 'o: M1 region');
hold on;
text(400, 560, '140Hz', 'Color', 'r');
hold on;
text(400, 540, '40Hz', 'Color', 'b');
hold on;

y_val = 600;
for mouse_f=fieldnames(mouse_color)'
    mouse_f = mouse_f{1};
    text(450, y_val, '--', 'Color', mouse_color.(mouse_f), 'FontWeight', 'bold');
    y_val = y_val - 20;
end
% Labels
xlabel('Days since surgery (days)');
ylabel('Amperage (uamp)');

% Set the axis
Multi_func.set_default_axis(gca);

saveas(gcf, [savepath 'Current/' '2d_currents_per_mouse_by_region_frequency.png']);
saveas(gcf, [savepath 'Current/' '2d_currents_per_mouse_by_region_frequency.pdf']);

%% Plotting current from days to surgery for each neuron separate plots for each brain region
% Red: 140Hz
% Blue: 40Hz
% Circle: M1
% Plus: V1

%A map that keeps track of how many repetitions there are
%of number of days since surgery and amperage
days_amp_map = [];
mouse_color = struct();

% Loop through regions
for f_region = fieldnames(region_data)'

    figure('Position', [0 0 1000 800]);
    
    f_region = f_region{1};

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        % Keep track of neurons plotted
        neurons_plotted = [];

        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            if ismember(nr, neurons_plotted)
                continue;
            end

            tokens = regexp(popul_data.neuron_name{nr}, '_', 'split');

            % Find all of the neurons that come from the same mouse
            same_m_nrs = find(contains(popul_data.neuron_name, tokens{1}) == 1);
            same_m_nrs_pts = [];
            neurons_plotted = cat(2, neurons_plotted, same_m_nrs);
            for i=same_m_nrs

                tokens = regexp(popul_data.neuron_name{i}, '_', 'split');

                % Calculate the days between imaging and surgery date
                surg_date = popul_data.surgery_date.(['m_' tokens{1}]);
                if isempty(surg_date)
                    continue;
                end
                tokens{3} = erase(tokens{3}, 'rec');

                surg_date = datetime(surg_date, 'InputFormat', 'yyyyMMdd');
                img_date = datetime(tokens{3}, 'InputFormat', 'yyyyMMdd');

                days_from_surg = days(img_date - surg_date);
            
                %DEBUG
                if days_from_surg > 400
                    disp(popul_data.neuron_name{i});
                end
                    
                % Parse out the current amperage
                fov_i = find(contains(tokens, 'FOV') == 1);
                amp_i = fov_i + 2;
                if contains(tokens{amp_i}, '.') == 1
                    amp_str = regexp(tokens{amp_i}, '\.', 'split');
                    amp_str = amp_str{1};
                else
                    amp_str = tokens{amp_i};
                end

                % Set labels based on brain region and stim frequency
                plot_str = '';
                if strcmp(f_region, 'r_M1') == 1
                    plot_str(end+1) = 'o';
                elseif strcmp(f_region, 'r_V1') == 1
                    plot_str(end+1) = '+';
                end

                % Set the frequency stuff
                if strcmp(f_stim, 'f_40') == 1
                    plot_str(end+1) = 'b';
                elseif strcmp(f_stim, 'f_140') == 1
                    plot_str(end+1) = 'r';
                end

                % Increment the marker size for repeated points
                if size(days_amp_map, 1) < days_from_surg | size(days_amp_map, 2) < str2num(amp_str)
                    days_amp_map(days_from_surg, str2num(amp_str)) = 1;
                else
                    days_amp_map(days_from_surg, str2num(amp_str)) = ...
                        1 + days_amp_map(days_from_surg, str2num(amp_str));
                end

                % Plot data point
                hold on;
                plot(days_from_surg, str2num(amp_str), plot_str, 'MarkerSize', 8 + days_amp_map(days_from_surg, str2num(amp_str)));
                hold on;

                % Keep track of points
                same_m_nrs_pts(:, end + 1) = [days_from_surg; str2num(amp_str)];
            end

            % Plot lines across all FOVs of the same mouse
            %TODO sort the x values from decreasing to increasing, makes figure easier to read
            if ~isempty(same_m_nrs_pts)
                
                [~, i] = sort(same_m_nrs_pts(1, :));
                if isfield(mouse_color, ['m_' tokens{1}])
                    color = mouse_color.(['m_' tokens{1}]);
                    plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i), '-', 'Color', color);
                else
                    plot_h = plot(same_m_nrs_pts(1, i), same_m_nrs_pts(2, i),  '-');
                    mouse_color.(['m_' tokens{1}]) = get(plot_h, 'Color');
                end
                hold on;
            end
        end
    end
    
    % Add legend label for each item
    text(400, 600, '+: V1 region');
    hold on;
    text(400, 580, 'o: M1 region');
    hold on;
    text(400, 560, '140Hz', 'Color', 'r');
    hold on;
    text(400, 540, '40Hz', 'Color', 'b');
    hold on;
    
    y_val = 600;
    for mouse_f=fieldnames(mouse_color)'
        mouse_f = mouse_f{1};
        text(450, y_val, '--', 'Color', mouse_color.(mouse_f), 'FontWeight', 'bold');
        y_val = y_val - 20;
    end
    % Labels
    xlabel('Days since surgery (days)');
    ylabel('Amperage (uamp)');
    
    % Set the axis
    Multi_func.set_default_axis(gca);
    
    % Title
    title(f_region);
        
    saveas(gcf, [savepath 'Current/' f_region '_2d_currents_per_mouse_by_frequency.png']);
    saveas(gcf, [savepath 'Current/' f_region '_2d_currents_per_mouse_by_frequency.pdf']);

end


%TODO make a 3D plot so that each mouse's data is in a different axis
% Try to reuse as much as possible from the top figure plotting
%% Generates a 3D plot of the mouse current stuff



%% Plot a similar plot as above, but a linear regression for all points in each condition
% I will have 4 plots where I overlay the points for each group together
% For example, top row will be 140Hz one panel and 40Hz on the other panel.
% Both brain regions are plotted in each panel
% Bottom row is the same idea, but with each panel being a unique brain region
% Mouse information will not be preserved here

%Specify a color based on region and stimulation frequency
color_dict.r_M1 = 'g';
color_dict.r_V1 = 'm';
color_dict.f_40 = 'b';
color_dict.f_140 = 'r';

figure('Position', [0 0 1000 800])

tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
% Loop through regions
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);
    stims = fieldnames(data_bystim);

    % Skip CA1 neurons from this plot
    if strcmp(f_region, 'r_CA1') == 1
        continue;
    end

    % Loop through stim frequencies
    for f_stim = stims'
        f_stim = f_stim{1};
        popul_data = data_bystim.(f_stim);
        
        nrs_pts = [];
        days_amp_map = [];
        % Loop through each neuron
        for nr=1:length(popul_data.neuron_name)
            tokens = regexp(popul_data.neuron_name{nr}, '_', 'split');
            
            % Grab surgery date
            surg_date = popul_data.surgery_date.(['m_' tokens{1}]);
            if isempty(surg_date)
                continue;
            end

            % Grab imaging date
            tokens{3} = erase(tokens{3}, 'rec');

            surg_date = datetime(surg_date, 'InputFormat', 'yyyyMMdd');
            img_date = datetime(tokens{3}, 'InputFormat', 'yyyyMMdd');

            days_from_surg = days(img_date - surg_date);

            % Parse out the current amperage
            fov_i = find(contains(tokens, 'FOV') == 1);
            amp_i = fov_i + 2;
            if contains(tokens{amp_i}, '.') == 1
                amp_str = regexp(tokens{amp_i}, '\.', 'split');
                amp_str = amp_str{1};
            else
                amp_str = tokens{amp_i};
            end
            
            amp = str2num(amp_str);

            % Add brain region data
            nrs_pts(:, end + 1) = [days_from_surg; amp];

            % Increment the marker size for repeated points
            if size(days_amp_map, 1) < days_from_surg | size(days_amp_map, 2) < amp
                days_amp_map(days_from_surg, amp) = 1;
            else
                days_amp_map(days_from_surg, amp) = ...
                    1 + days_amp_map(days_from_surg, amp);
            end

            % plot point in respective panel
        
            % Add points to respective plots
            if strcmp(f_region, 'r_M1') == 1
                nexttile(1);
            elseif strcmp(f_region, 'r_V1') == 1
                nexttile(2);
            end

            % Perform plotting for the brain region
            hold on;
            plot(days_from_surg, amp, [color_dict.(f_stim) '.'], 'MarkerSize',...
                10 + 1.5*days_amp_map(days_from_surg, amp));
            
            if strcmp(f_stim, 'f_40') == 1
                nexttile(4);
            elseif strcmp(f_stim, 'f_140') == 1
                nexttile(3);
            end

            % Perform plotting for the stimulation frequency
            hold on;
            plot(days_from_surg, amp, [color_dict.(f_region) '.'], 'MarkerSize',...
                10 + 1.5*days_amp_map(days_from_surg, amp));

        end

        % Perform regression on each group of plots
        fit_results = polyfit(nrs_pts(1, :), nrs_pts(2, :), 1);
        fit_y = polyval(fit_results, nrs_pts(1, :));
        pearson_coeff = corrcoef(nrs_pts(1, :), nrs_pts(2, :));


        % Add points to respective plots
        if strcmp(f_region, 'r_M1') == 1
            nexttile(1);
        elseif strcmp(f_region, 'r_V1') == 1
            nexttile(2);
        end

        % Perform plotting for the brain region
        hold on;
        plot(nrs_pts(1, :), fit_y, [color_dict.(f_stim) '-']);
        hold on;
        text(450, max(nrs_pts(2, :))/2 - 10,...
            [f_stim ' ' num2str(pearson_coeff(1, 2))], ...
            'Interpreter', 'none', 'Color', color_dict.(f_stim));
        xlim([0 500]);
        ylim([0 350]);
        Multi_func.set_default_axis(gca);

        if strcmp(f_stim, 'f_40') == 1
            nexttile(4);
        elseif strcmp(f_stim, 'f_140') == 1
            nexttile(3);
        end

        % Perform plotting for the stimulation frequency
        hold on;
        plot(nrs_pts(1, :), fit_y, [color_dict.(f_region) '-']);
        hold on;
        text(450, max(nrs_pts(2, :))/2 - 10,...
            [f_region ' ' num2str(pearson_coeff(1, 2))], ...
            'Interpreter', 'none', 'Color', color_dict.(f_region));
        xlim([0 500]);
        ylim([0 350]);
        Multi_func.set_default_axis(gca);

    end
end

nexttile(1);
title('M1 Neurons');
xlabel('Time imaged from surgery (days)');
ylabel('Current amplitude (uamp)');
nexttile(2);
title('V1 Neurons');
xlabel('Time imaged from surgery (days)');
ylabel('Current amplitude (uamp)');
nexttile(3);
title('140Hz neurons');
xlabel('Time imaged from surgery (days)');
ylabel('Current amplitude (uamp)');
nexttile(4);
title('40Hz neurons');
xlabel('Time imaged from surgery (days)');
ylabel('Current amplitude (uamp)');


saveas(gcf, [savepath 'Current/' 'by_condition_currents.png']);
saveas(gcf, [savepath 'Current/' 'by_condition_currents.pdf']);




%% -- (Deprecated) --

all_currents = [];
all_Fs = [];

% Loop through and grab each neuron's current amplitude
for f_region = fieldnames(region_data)'
    f_region = f_region{1};
    data_bystim = region_data.(f_region);

   % % Skip CA1 neurons from this plot
   % if strcmp(f_region, 'r_CA1') == 1
   %     continue;
   % end

    for f_stim = fieldnames(data_bystim)'
        f_stim = f_stim{1};

        % Loop through each neuron
        for nr = 1:length(data_bystim.(f_stim).neuron_name)
            cur_nr_name = strsplit(data_bystim.(f_stim).neuron_name{nr}, '_');
            
            % Get the current from the matfile
            f_mouse = ['m_' cur_nr_name{1}];
            f_rec = ['r_' cur_nr_name{3}];               
            f_neuron= ['n_' cur_nr_name{4}];
            f_stim = ['f_' cur_nr_name{5}];

            neuron_num = str2num(erase(cur_nr_name{4}, 'FOV'));
            
            [region_data] = Multi_func.create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
            
            % Saving the current at neuron number index
            region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(cur_nr_name{6});
            
        end
    end
end

%for f_region = fieldnames(region_matfiles)'
%    f_region = f_region{1};
%
%
%    % Select matfiles by stim specific conditions for all regions
%    %[matfile_stim] = stim_cond(all_matfiles); 
%    % Select matfiles by stim condition for given region
%    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);
%
%    % Loop through each stimulation condition
%    for f_stim = fieldnames(matfile_stim)'
%        f_stim = f_stim{1};
%        matfiles = matfile_stim.(f_stim).names;    
%    
%        % Initialize field subthreshold array
%        data_bystim.(f_stim) = struct();
%
%        % Loop through each matfile of the current stimulation condition
%        for matfile = matfiles
%            matfile = matfile{1};
%            % Read in the mat file of the current condition
%            data = load([pv_data_path matfile]);
%            trial_idxs = find(~cellfun(@isempty, data.align.trial));
%            trial_data = data.align.trial{trial_idxs(1)};    
%         
%            
%            % Loop through each ROI
%            for roi_idx=1:size(trial_data.detrend_traces, 2)
%                %Determine whether this roi is to be ignored for this particular trial
%                ri = strsplit(matfile, '_');
%                try
%                    trial_ignr_list = ignore_trial_dict.(['mouse_' ri{1}]).(['rec_' erase(ri{3}, 'rec')]).(ri{4}).(['f_' ri{5}]).(['ROI' num2str(roi_idx)]);
%                catch
%                    trial_ignr_list = [];
%                end
%                
%                % Remove trials from trial idx list
%                trial_idxs = setdiff(trial_idxs, trial_ignr_list);
%
%                % Skip if there is at most 2 trials
%                if length(trial_idxs) <= 2
%                    continue;
%                end
%                
%                % Get the current from the matfile
%                f_mouse = ['m_' ri{1}];
%                f_rec = ['r_' ri{3}];               
%                f_neuron= ['n_' ri{4}];
%                f_stim = ['s_' ri{5}];
%                
%                % convert FOV number to single number
%                neuron_num = str2num(erase(ri{4}, 'FOV'));
%
%                [region_data] = create_struct(region_data, f_region, f_mouse, f_rec, f_neuron, f_stim);
%                
%                % Saving the current at neuron number index
%                region_data.(f_region).(f_mouse).(f_rec).(f_stim).currents(neuron_num) = str2num(ri{6});
%
%                % Save all of the currents
%                all_currents(end + 1) = str2num(ri{6});
%
%            end
%        end
%    end
%end


%% Plot the current by mouse name
% Clear all zeros in the current array
all_currents(find(all_currents == 0)) = [];

disp(['Average currents ' num2str(mean(all_currents))]);
disp(['Std currents ' num2str(std(all_currents))]);

% Plot all of the currents as a single plot
figure;
tick_labels = {};
tick_spots = [];
i = 1;
for f_region = {'r_M1'} %fieldnames(region_data)'
    f_region = f_region{1};

    for f_mouse = fieldnames(region_data.(f_region))'
        f_mouse = f_mouse{1};
        
        % Fields with leading m indicate a mouse entry
        if strcmp(f_mouse(1), 'm') ~= 1
            continue;
        end
    
        for f_rec = fieldnames(region_data.(f_region).(f_mouse))'
            f_rec = f_rec{1};
 
            % Check if each stimulation condition exists for the neuron
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_40')
                current = region_data.(f_region).(f_mouse).(f_rec).f_40.currents;
                current(find(current == 0)) = [];
                
                plot(i + rand(1, length(current)), current, '.r');
                hold on;
            end
            i = i + 1;
            
            if isfield(region_data.(f_region).(f_mouse).(f_rec), 'f_140')
                current = region_data.(f_region).(f_mouse).(f_rec).f_140.currents;
                current(find(current == 0)) = [];
                plot(i + rand(1, length(current)), current, '.b');
                hold on;
            end

            % Add the mouse and recording date for the plot
            tick_labels{end + 1} = [erase(f_mouse, 'm_')]; % ' ' erase(f_rec, 'r_')
            tick_spots(end + 1) = i;

            i = i + 2;
        end 

        i = i + 3;
    end
end
xticks(tick_spots);
xticklabels(tick_labels);
ylabel('uamp');

saveas(gcf, [savepath 'Current/' 'all_mice_currents.png']);
saveas(gcf, [savepath 'Current/' 'all_mice_currents.pdf']);

% Plot the Normalized Vm Change vs Current amplitude
