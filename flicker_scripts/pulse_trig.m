clc;
clear all;
close all;

%%
f = filesep;

extra_trace = 3;

addpath('..');

save_data_file = [interm_data_path 'v1_flicker.mat'];

data = load(save_data_file);

%% Store and calculate the visual pulse-triggered 
stim_time = [1, 2];

% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;

    % Find all of the pulse onset times
    diff_arr = diff(popul_data.flicker_raster);
    diff_arr = [0, diff_arr];
    [~, pulse_idxs] = findpeaks(diff_arr);

    % Calculate idxs that correspond to the width of a full flicker wave
    pulse_width = floor(mean(diff(pulse_idxs)));
    data.(f_stim).pulse_width = pulse_width;

    % Plot to ensure that flicker onset idxs are in the right place
    %figure;
    %plot(popul_data.flicker_raster);
    %hold on;
    %plot(pulse_idxs, popul_data.flicker_raster(pulse_idxs),  '|')
    %hold on;
    %plot([pulse_idxs(1) + [0  pulse_width]], [1 1], 'LineWidth', 10);

    % Store each period pulses
    all_pre_pulse_vm = [];
    all_stim_pulse_vm = [];
    all_post_pulse_vm = [];

    % Loop through each neuron
    for f_nr = fieldnames(popul_data.raw_vm)'
        f_nr = f_nr{1};
        
        % Loop through each trial
        for i = 1:size(popul_data.raw_vm.(f_nr), 2)
            tr_vm = popul_data.raw_vm.(f_nr)(:, i);
                       
            % Spike normalize the trace
            sp_amp = popul_data.spike_amp_raster.(f_nr)(:, i);
            avg_sp_amp = mean(sp_amp, 'omitnan');

            % Replace if NaN
            if isnan(avg_sp_amp)
                avg_sp_amp = 1;
            end
            
            tr_norm_vm = tr_vm/avg_sp_amp;
            
            % Loop through each pulse
            for p_idx = pulse_idxs
                cur_pulse = tr_norm_vm(p_idx - extra_trace:p_idx+pulse_width+extra_trace);
                
                % Pulse subtract the current pulse
                cur_pulse = cur_pulse - cur_pulse(extra_trace + 1);

                % Check what period the flicker idx is
                
                % Pre-stim
                if timeline(p_idx) < stim_time(1)
                    all_pre_pulse_vm = cat(2, all_pre_pulse_vm, cur_pulse(:));
                
                % During stim
                elseif timeline(p_idx) >= stim_time(1) && timeline(p_idx) < stim_time(2)
                    all_stim_pulse_vm = cat(2, all_stim_pulse_vm, cur_pulse(:));
                % Post stim
                else
                    all_post_pulse_vm = cat(2, all_post_pulse_vm, cur_pulse(:));
                end
            end
        end
    end

    % Save pulses to struct
    data.(f_stim).pre_pulse_vm = all_pre_pulse_vm;
    data.(f_stim).stim_pulse_vm = all_stim_pulse_vm;
    data.(f_stim).post_pulse_vm = all_post_pulse_vm;
    
end

%% Compare point-by-point with t-test
% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
 
    popul_data = data.(f_stim);

    % Store the main window pulse points
    cur_pre_pulses = popul_data.pre_pulse_vm(extra_trace + 1:end - extra_trace, :);
    cur_stim_pulses = popul_data.stim_pulse_vm(extra_trace + 1:end - extra_trace, :);
    cur_post_pulses = popul_data.post_pulse_vm(extra_trace + 1:end - extra_trace, :);
    
    % Adjust the alpha value
    p_val_adj = 0.05/size(cur_pre_pulses, 1);
    data.(f_stim).adj_pval = p_val_adj;

    % Pre to stim comparison
    [~, p_vals] = ttest(cur_pre_pulses', cur_stim_pulses');
    % Save the p_vals that were less than the adjusted p_value
    data.(f_stim).ttest_pre_stim_sig_idx = find(p_vals < p_val_adj);
    
    % Pre to post comparison
    [~, p_vals] = ttest(cur_pre_pulses', cur_post_pulses');
    % Save the p_vals that were less than the adjusted p_value
    data.(f_stim).ttest_pre_post_sig_idx = find(p_vals < p_val_adj);

    % Stim to post comparison
    [~, p_vals] = ttest(cur_stim_pulses', cur_post_pulses');
    % Save the p_vals that were less than the adjusted p_value
    data.(f_stim).ttest_stim_post_sig_idx = find(p_vals < p_val_adj);
end

%% Plot pulse-triggered curves together

% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time(1:popul_data.pulse_width + 2*extra_trace + 1)';
    timeline = timeline - timeline(1) - extra_trace*(mean(diff(timeline)));

    % Plot the pre, stim, post at once
    avg_pre_vm = mean(popul_data.pre_pulse_vm, 2);
    avg_stim_vm = mean(popul_data.stim_pulse_vm, 2);
    avg_post_vm = mean(popul_data.post_pulse_vm, 2);

    % Calculate the SEM for each pulse vm
    sem_pre_vm = std(popul_data.pre_pulse_vm, 0, 2)./sqrt(size(popul_data.pre_pulse_vm, 2));
    sem_stim_vm = std(popul_data.stim_pulse_vm, 0, 2)./sqrt(size(popul_data.stim_pulse_vm, 2));
    sem_post_vm = std(popul_data.post_pulse_vm, 0, 2)./sqrt(size(popul_data.post_pulse_vm, 2));
    
    % Plot all together
    figure;
    %tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
    %-- Pre to Stim
    nexttile;
    fill_h = fill([timeline; flip(timeline)], [avg_pre_vm + sem_pre_vm; flipud(avg_pre_vm - sem_pre_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_pre_vm, '-b', 'DisplayName', 'Pre', 'color', Multi_func.base_color);
    hold on;

    fill_h = fill([timeline; flip(timeline)], [avg_stim_vm + sem_stim_vm; flipud(avg_stim_vm - sem_stim_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_stim_vm, '-g', 'DisplayName', 'Stim', 'color', Multi_func.stim_color );
    hold on;
    
    % Plot the significant between pre to stim
    sig_idx = popul_data.ttest_pre_stim_sig_idx + extra_trace;
    plot(timeline(sig_idx), avg_stim_vm(sig_idx) , '.b', 'MarkerSize', 10, 'HandleVisibility', 'off');
    hold on;

    % Plot the dotted midline
    xline(timeline(popul_data.pulse_width/2), '--');

    Multi_func.set_default_axis(gca);
    legend('Location', 'bestoutside');
    title('Pre to Stim');


    %-- Stim to Post
    nexttile;
    fill_h = fill([timeline; flip(timeline)], [avg_stim_vm + sem_stim_vm; flipud(avg_stim_vm - sem_stim_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_stim_vm, '-b', 'DisplayName', 'Stim', 'color', Multi_func.stim_color);
    hold on;

    fill_h = fill([timeline; flip(timeline)], [avg_post_vm + sem_post_vm; flipud(avg_post_vm - sem_post_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_post_vm, '-g', 'DisplayName', 'Post', 'color', Multi_func.post_color);
    hold on;
    
    % Plot the significant between Stim to Post
    sig_idx = popul_data.ttest_stim_post_sig_idx + extra_trace;
    plot(timeline(sig_idx), avg_post_vm(sig_idx) , '.b', 'MarkerSize', 10, 'HandleVisibility', 'off');
    hold on;

    % Plot the dotted midline
    xline(timeline(popul_data.pulse_width/2), '--');

    Multi_func.set_default_axis(gca);
    legend('Location', 'bestoutside');
    title('Stim to Post');


    %-- Pre to Post
    nexttile;
    fill_h = fill([timeline; flip(timeline)], [avg_pre_vm + sem_pre_vm; flipud(avg_pre_vm - sem_pre_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_pre_vm, '-b', 'DisplayName', 'Pre', 'color', Multi_func.base_color);
    hold on;

    fill_h = fill([timeline; flip(timeline)], [avg_post_vm + sem_post_vm; flipud(avg_post_vm - sem_post_vm)], [0.5 0.5 0.5], ...
        'HandleVisibility', 'off');
    Multi_func.set_fill_properties(fill_h);
    hold on;
    plot(timeline, avg_post_vm, '-g', 'DisplayName', 'Post', 'color', Multi_func.post_color);
    hold on;
    
    % Plot the significant between Stim to Post
    sig_idx = popul_data.ttest_pre_post_sig_idx + extra_trace;
    plot(timeline(sig_idx), avg_post_vm(sig_idx) , '.b', 'MarkerSize', 10, 'HandleVisibility', 'off');
    hold on;

    % Plot the dotted midline
    xline(timeline(popul_data.pulse_width/2), '--');

    Multi_func.set_default_axis(gca);
    legend('Location', 'bestoutside');
    title('Pre to Post');

    % Plot the flicker window
    %xline([], 'Color', [170, 176, 97]./255, 'LineWidth', 0.5);

    %fill_h = fill([timeline; flip(timeline)], [avg_post_vm + sem_post_vm; flipud(avg_post_vm - sem_post_vm)], [0.5 0.5 0.5]);
    %Multi_func.set_fill_properties(fill_h);
    %hold on;
    %plot(timeline, avg_post_vm, '-m');
    
    sgtitle(f_stim, 'Interpreter', 'none');
end
