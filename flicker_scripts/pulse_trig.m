clc;
clear all;
close all;

%%
f_sep = filesep;

extra_trace = 3;

addpath('..');

% This structure is taken from the python export file
interm_data_path = [f_sep 'home' f_sep 'pierfier' f_sep 'Projects' f_sep ...
      'Pierre Fabris' f_sep 'PV DBS neocortex' f_sep 'Interm_Data' f_sep];

%save_data_file = [interm_data_path 'v1_flicker.mat'];
save_data_file = [interm_data_path 'v1_flicker_aligned_dat.mat'];

data = load(save_data_file);

%% Calculate the idealized frame time
total_frames = 2500;
total_time = 5.0258; % in Sec
soft_Fs = total_frames/total_time;
exp_time = 1/soft_Fs;

front_frame_drop = 14;

% Start time taken as difference from the ephys time
true_time = 1000*exp_time*[front_frame_drop+1:total_frames] - 1; % Convert to ms

%% Store and calculate the visual pulse-triggered 
stim_time = [1, 2];

% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);
    timeline = true_time; %timeline = popul_data.interp_time;

    % Find all of the pulse onset times
    diff_arr = diff(popul_data.flicker_raster); %TODO may need to change a theoretical flicker raster
    diff_arr = [0, diff_arr];
    [~, pulse_idxs] = findpeaks(diff_arr);

    % Calculate idxs that correspond to the width of a full flicker wave
    pulse_width = floor(mean(diff(pulse_idxs), 'omitnan'));
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
    
    % Store each period pulses
    avg_pre_pulse_vm = [];
    avg_stim_pulse_vm = [];
    avg_post_pulse_vm = [];

    % Loop through each neuron
    for f_nr = fieldnames(popul_data.raw_vm)'
        f_nr = f_nr{1};
        
        % Store each period pulses
        nr_pre_pulse_vm = [];
        nr_stim_pulse_vm = [];
        nr_post_pulse_vm = [];

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
                    nr_pre_pulse_vm = cat(2, nr_pre_pulse_vm, cur_pulse(:));
                
                % During stim
                elseif timeline(p_idx) >= stim_time(1) && timeline(p_idx) < stim_time(2)
                    all_stim_pulse_vm = cat(2, all_stim_pulse_vm, cur_pulse(:));
                    nr_stim_pulse_vm = cat(2, nr_stim_pulse_vm, cur_pulse(:));
                
                % Post stim
                else
                    all_post_pulse_vm = cat(2, all_post_pulse_vm, cur_pulse(:));
                    nr_post_pulse_vm = cat(2, nr_post_pulse_vm, cur_pulse(:));
                end
            end
        end

        avg_pre_pulse_vm = cat(2, avg_pre_pulse_vm, mean(nr_pre_pulse_vm, 2, 'omitnan'));
        avg_stim_pulse_vm = cat(2, avg_stim_pulse_vm, mean(nr_stim_pulse_vm, 2, 'omitnan'));
        avg_post_pulse_vm = cat(2, avg_post_pulse_vm, mean(nr_post_pulse_vm, 2, 'omitnan'));
    end

    % Save pulses to struct
    data.(f_stim).pre_pulse_vm = all_pre_pulse_vm;
    data.(f_stim).stim_pulse_vm = all_stim_pulse_vm;
    data.(f_stim).post_pulse_vm = all_post_pulse_vm;
    
    data.(f_stim).avg_pre_pulse_vm =  avg_pre_pulse_vm;
    data.(f_stim).avg_stim_pulse_vm = avg_stim_pulse_vm;
    data.(f_stim).avg_post_pulse_vm = avg_post_pulse_vm;
end

%% Calculate pulse-triggered average using the individual frame times and flicker onset and offset time
extra_trace = 0;
stim_time = [1, 2];
% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);

    % Store each period pulses
    avg_pre_pulse_vm = [];
    avg_stim_pulse_vm = [];
    avg_post_pulse_vm = [];

    avg_pre_offset_t = [];
    avg_stim_offset_t = [];
    avg_post_offset_t = [];

    all_pulse_frame_time = [];

    % Loop through each neuron
    for f_nr = fieldnames(popul_data.raw_vm)'
        f_nr = f_nr{1};
 
        % Store each period pulses
        nr_pre_pulse_vm = [];
        nr_stim_pulse_vm = [];
        nr_post_pulse_vm = [];
        
        % Store the pulse period swtich off
        nr_pre_offset_t = [];
        nr_stim_offset_t = [];
        nr_post_offset_t = [];

        % DEBUG the trace timepoint between the camera frames and the theoretical time
        %figure;
        
        % TODO plot the flicker timepoints and maybe have it so that all of the figures are saved somewhere
        % for quick checking
        pulse_width = diff(popul_data.flicker_times.(f_nr) );

        % Loop through each trial
        for i = 1:size(popul_data.raw_vm.(f_nr), 2)
            tr_vm = popul_data.raw_vm.(f_nr)(:, i);
            
            % Spike normalize the trace
            sp_amp = popul_data.sp_amp_raster.(f_nr)(:, i);
            avg_sp_amp = mean(sp_amp, 'omitnan');

            % Replace if NaN
            if isnan(avg_sp_amp)
                avg_sp_amp = 1;
            end
            
            tr_norm_vm = tr_vm/avg_sp_amp;

            % DEBUG plot the trial with the true time and the frame time
            %plot(true_time, tr_norm_vm + i, '-r'); % , 'LineWidth', 3
            %hold on;
            %f_time = popul_data.camera_frame_times.(f_nr);
            %plot(f_time(:, i), tr_norm_vm + i, '-b');
            %hold on;
            %flick_time = popul_data.flicker_times.(f_nr);
            %plot(flick_time(:, i), repmat(max(tr_norm_vm + i), length(flick_time(:, i)), 1), '|');

            % Loop through each pulse onset times and get the timepoints within each flicker on and off time
            pulse_on_time = popul_data.flicker_times.(f_nr);
            pulse_off_time = popul_data.flicker_off_times.(f_nr);
            trace_time = popul_data.camera_frame_times.(f_nr);
            for p_idx = 1:length(pulse_on_time) - 1
                pulse_on = pulse_on_time(p_idx, i);
                pulse_next = pulse_on_time(p_idx + 1, i);
                
                % Store the pulse switch within the duty cycle
                pulse_offset = pulse_off_time(p_idx, i);

                % Store the current pulse
                cur_pulse = tr_norm_vm(trace_time(:, i) > pulse_on &  trace_time(:, i) < pulse_next);
                cur_pulse = cur_pulse - cur_pulse(1);
                
                % Store the frame time for the current pulse window
                cur_pulse_frame_window = trace_time(trace_time(:, i) > pulse_on & trace_time(:, i) <= pulse_next, i);
                cur_pulse_frame_window = cur_pulse_frame_window - pulse_on; % cur_pulse_frame_window(1) 
            
                frame_diff = mean(diff(cur_pulse_frame_window));
                num_frames = floor(mean(pulse_width, 'all')./frame_diff);

                % Chop end piece if the number of frames are too large for the window
                cur_pulse_frame_window = cur_pulse_frame_window(1:num_frames);
                cur_pulse = cur_pulse(1:num_frames);

                % Consolidate all pulse windows into one array
                all_pulse_frame_time = horzcat_pad(all_pulse_frame_time, cur_pulse_frame_window(:));

                %DEBUG check if the array is in sorted order
                %if length(cur_pulse) == 63 %  ~isnan(all_pulse_frame_time(end, end))
                %    figure, plot(cur_pulse);
                %    input('En');
                %    %error('Testing single pulses');
                %end

                % Pre-stim
                if pulse_on < stim_time(1)
                    nr_pre_pulse_vm = horzcat_pad(nr_pre_pulse_vm, cur_pulse(:));
                    nr_pre_offset_t = cat(2, nr_pre_offset_t, pulse_offset - pulse_on);
                % During stim
                elseif pulse_on > stim_time(1) & pulse_on < stim_time(2)
                    nr_stim_pulse_vm = horzcat_pad(nr_stim_pulse_vm, cur_pulse(:));
                    nr_stim_offset_t = cat(2, nr_stim_offset_t, pulse_offset - pulse_on);
                % Post stim
                else
                    nr_post_pulse_vm = horzcat_pad(nr_post_pulse_vm, cur_pulse(:));
                    nr_post_offset_t = cat(2, nr_post_offset_t, pulse_offset - pulse_on);
                end
            end
        end
        
        % Store the average Vm
        avg_pre_pulse_vm = horzcat_pad(avg_pre_pulse_vm, mean(nr_pre_pulse_vm, 2, 'omitnan'));
        avg_stim_pulse_vm = horzcat_pad(avg_stim_pulse_vm, mean(nr_stim_pulse_vm, 2, 'omitnan'));
        avg_post_pulse_vm = horzcat_pad(avg_post_pulse_vm, mean(nr_post_pulse_vm, 2, 'omitnan'));

        avg_pre_offset_t = cat(2, avg_pre_offset_t, mean(nr_pre_offset_t, 'omitnan'));
        avg_stim_offset_t = cat(2, avg_stim_offset_t, mean(nr_stim_offset_t, 'omitnan'));
        avg_post_offset_t = cat(2, avg_post_offset_t, mean(nr_post_offset_t, 'omitnan'));
    end

    % Save all of the data for the respective frequency
    data.(f_stim).all_pulse_frame_time = all_pulse_frame_time;

    data.(f_stim).avg_pre_pulse_vm = avg_pre_pulse_vm;
    data.(f_stim).avg_stim_pulse_vm = avg_stim_pulse_vm;
    data.(f_stim).avg_post_pulse_vm = avg_post_pulse_vm;

    data.(f_stim).pre_offset_t = avg_pre_offset_t;
    data.(f_stim).stim_offset_t = avg_stim_offset_t;
    data.(f_stim).post_offset_t = avg_post_offset_t;
end

%% Compare point-by-point with t-test
% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
 
    popul_data = data.(f_stim);

    % Store the main window pulse points
    cur_pre_pulses = popul_data.avg_pre_pulse_vm(extra_trace + 1:end - extra_trace, :);
    cur_stim_pulses = popul_data.avg_stim_pulse_vm(extra_trace + 1:end - extra_trace, :);
    cur_post_pulses = popul_data.avg_post_pulse_vm(extra_trace + 1:end - extra_trace, :);
    
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

%% Plot pulse-triggered curves together using the true, "theoretical" time

% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);
    %timeline = popul_data.interp_time(1:popul_data.pulse_width + 2*extra_trace + 1)';
    timeline = true_time(1:popul_data.pulse_width + 2*extra_trace + 1)';
    timeline = timeline - timeline(1) - extra_trace*(mean(diff(timeline)));

    % Plot the pre, stim, post at once
    avg_pre_vm = mean(popul_data.avg_pre_pulse_vm, 2, 'omitnan');
    avg_stim_vm = mean(popul_data.avg_stim_pulse_vm, 2, 'omitnan');
    avg_post_vm = mean(popul_data.avg_post_pulse_vm, 2, 'omitnan');

    % Calculate the SEM for each pulse vm
    sem_pre_vm = std(popul_data.avg_pre_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_pre_pulse_vm, 2));
    sem_stim_vm = std(popul_data.avg_stim_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_stim_pulse_vm, 2));
    sem_post_vm = std(popul_data.avg_post_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_post_pulse_vm, 2));
    
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

%% Plot pulse-triggered curves from the camera frame triggers
%TODO need to change all of the variables here for 

% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    
    popul_data = data.(f_stim);
    %timeline = popul_data.interp_time(1:popul_data.pulse_width + 2*extra_trace + 1)';
    %timeline = true_time(1:popul_data.pulse_width + 2*extra_trace + 1)';
    %timeline = timeline - timeline(1) - extra_trace*(mean(diff(timeline)));

    timeline = mean(popul_data.all_pulse_frame_time, 2, 'omitnan') + 2/1000; %; %  

    % Plot the pre, stim, post at once
    avg_pre_vm = mean(popul_data.avg_pre_pulse_vm, 2, 'omitnan');
    avg_stim_vm = mean(popul_data.avg_stim_pulse_vm, 2, 'omitnan');
    avg_post_vm = mean(popul_data.avg_post_pulse_vm, 2, 'omitnan');

    % Calculate the SEM for each pulse vm
    sem_pre_vm = std(popul_data.avg_pre_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_pre_pulse_vm, 2));
    sem_stim_vm = std(popul_data.avg_stim_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_stim_pulse_vm, 2));
    sem_post_vm = std(popul_data.avg_post_pulse_vm, 0, 2, 'omitnan')./sqrt(size(popul_data.avg_post_pulse_vm, 2));
    
    % Plot all together
    figure;
    %tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 20, 3.5, 5]);
    %tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    tiledlayout(2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    
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

    % Plot the dotted midline from the offset_t
    xline(popul_data.pre_offset_t, 'r--', 'HandleVisibility', 'off');
    hold on;
    xline(popul_data.stim_offset_t, 'b--', 'HandleVisibility', 'off');

    Multi_func.set_default_axis(gca);
    legend('Location', 'bestoutside');
    title('Pre to Stim');


   % %-- Stim to Post
   % nexttile;
   % fill_h = fill([timeline; flip(timeline)], [avg_stim_vm + sem_stim_vm; flipud(avg_stim_vm - sem_stim_vm)], [0.5 0.5 0.5], ...
   %     'HandleVisibility', 'off');
   % Multi_func.set_fill_properties(fill_h);
   % hold on;
   % plot(timeline, avg_stim_vm, '-b', 'DisplayName', 'Stim', 'color', Multi_func.stim_color);
   % hold on;

   % fill_h = fill([timeline; flip(timeline)], [avg_post_vm + sem_post_vm; flipud(avg_post_vm - sem_post_vm)], [0.5 0.5 0.5], ...
   %     'HandleVisibility', 'off');
   % Multi_func.set_fill_properties(fill_h);
   % hold on;
   % plot(timeline, avg_post_vm, '-g', 'DisplayName', 'Post', 'color', Multi_func.post_color);
   % hold on;
   % 
   % % Plot the significant between Stim to Post
   % sig_idx = popul_data.ttest_stim_post_sig_idx + extra_trace;
   % plot(timeline(sig_idx), avg_post_vm(sig_idx) , '.b', 'MarkerSize', 10, 'HandleVisibility', 'off');
   % hold on;

   % % Plot the dotted midline
   % xline(timeline(popul_data.pulse_width/2), '--');

   % Multi_func.set_default_axis(gca);
   % legend('Location', 'bestoutside');
   % title('Stim to Post');


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
    xline(popul_data.pre_offset_t, 'r--', 'HandleVisibility', 'off');
    hold on;
    xline(popul_data.post_offset_t, 'b--', 'HandleVisibility', 'off');

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

%% Plot a handful of the pulse-triggered examples to try to observe general trends
test_pulses = 100;
% Loop through each stim frequency
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);

    % Grab some pulses
    pulses = popul_data.stim_pulse_vm(:, 1:test_pulses);

    % Plot with separation
    figure;
    plot(pulses);
    hold on;
    xline(popul_data.pulse_width/2, '--');
    
    title(f_stim, 'Interpreter', 'none');
end
