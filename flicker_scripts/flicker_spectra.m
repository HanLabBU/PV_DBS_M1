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

% Index for 'nr_mod' array the description of each  is storing
mod_str = {'inc', 'dec', 'unc'};

%res =  24; % Change to 24 for higher res, 10 for low res
voices_per_octaves = {'vpo_10', 'vpo_17', 'vpo_24'};

%% Calculate parameters for recording

total_frames = 2500;
total_time = 5.0258; % in Sec
soft_Fs = total_frames/total_time;

front_frame_drop = 14;

%%
save_data_file = [interm_data_path 'v1_flicker.mat'];
%save_data_file = [interm_data_path 'v1_flicker_aligned_dat.mat'];

data = load(save_data_file);

%% Calculate the power spectra for each trial
freqLimits = [0 200];
Fs = 500;
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    
    % Initialize power spectra fields
    %data.(f_stim).spec_pow = {};
    %data.(f_stim).spec_f = {};

    % Loop through each neuron
    for f_nr = fieldnames(popul_data.raw_vm)'
        f_nr = f_nr{1};
        
        % Store the power spectra for the current neuron
        nr_spec_pow = struct;
        nr_spec_f = struct;

        for vpo = voices_per_octaves
            vpo = vpo{1};
            nr_spec_pow.(vpo) = [];
            nr_spec_f.(vpo) = [];
        end

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
            
            % Calculate frequencies with different frequency step resolutions
            for vpo = voices_per_octaves
                vpo = vpo{1};
                
                fb = cwtfilterbank(SignalLength=length(tr_norm_vm), ...
                            SamplingFrequency=Fs, ...
                            FrequencyLimits=freqLimits, ...
                            VoicesPerOctave=str2num(vpo(5:end)));
                [wt, frq] = cwt(tr_vm, FilterBank=fb);
                
                % Concatenate all of the power data
                nr_spec_pow.(vpo)(:, :, end + 1) = wt;
                nr_spec_f.(vpo)(:, :, end + 1) = frq;
            end
        end

        % Remove empty element
        for vpo = voices_per_octaves
            vpo = vpo{1};
            nr_spec_pow.(vpo)(:, :, 1) = [];
            nr_spec_f.(vpo)(:, :, 1) = [];
        end

        % Add to full structure
        data.(f_stim).spec_pow.(f_nr) = nr_spec_pow;
        data.(f_stim).spec_f.(f_nr) = nr_spec_f;
    end
end

%% Average and Plot the power spectra of the Vm
flick_base_tmps = [0 1];
stim_timestamps = [1 2];

vpo = 'vpo_24';

for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;

    flicker_raster = popul_data.flicker_raster;

    get_base_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp < stim_tmstp(1));
    get_stim_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp >= stim_tmstp(1) & ...
                                               tr_tmstmp <= stim_tmstp(end));

    calc_time_pow = @(trial_spec, time_idxs) mean(trial_spec(:, time_idxs, :), 2);
    
    calc_trial_spec = @(trial_spec, base_pow, stim_pow) (trial_spec - base_pow)./(base_pow + stim_pow);
    
    calc_nr_spec = @(trial_spec, tr_tmstmp, stim_tmstp) mean(calc_trial_spec(trial_spec, ...
        calc_time_pow(trial_spec, get_base_idxs(tr_tmstmp, stim_tmstp)),  ...
        calc_time_pow(trial_spec, get_stim_idxs(tr_tmstmp, stim_tmstp))), 3, 'omitnan');

    cur_spec_pow = cellfun(@(f_nr) calc_nr_spec(abs(popul_data.spec_pow.(f_nr).(vpo)), ...
                            interp_time, ...
                            flick_base_tmps), ...
                            fieldnames(popul_data.spec_pow)', 'UniformOutput', false);
    
    cur_spec_pow = cat(3, cur_spec_pow{:});
    
    freqs_p = mean(popul_data.spec_f.n_0.(vpo), 3, 'omitnan');
    
    figure; % 3.3672 cm height, width= 5 cm
    surface(interp_time, ... 
            freqs_p, ...
            mean(cur_spec_pow, 3, 'omitnan'), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
    %colormap(jet*.8);
    colormap(Multi_func.red_purple_blue_color);
    hold on;

    % Plot the stim protocol
    plot(interp_time, 1.*flicker_raster + max(freqs_p) + 1, '-b');
    hold on;
    plot(stim_timestamps, 0.5.*[1 1] + max(freqs_p) + 2, 'color', Multi_func.dbs_color);

    % Color bar and labels
    a = colorbar;
    a.Label.String = "Relative Power";
    title([f_stim ' Hz'], 'Interpreter', 'none');
    xlabel('Time from onset (ms)');
    ylabel('Frequency (Hz)');
    
    % Set limits from timeline and the frequency range
    xlim([min(interp_time), max(interp_time)]);
    ylim([min(freqs_p), max(freqs_p)]); %  8 + 4 TODO change here for the range

    % Set the default properties for the plot
    ax = gca;
    Multi_func.set_default_axis(ax);
    ylimits = ylim(ax);
    xlimits = xlim(ax);
    caxis_lts = caxis;

    % Specify dimensions
    ax.Units = 'centimeters';
    ax.InnerPosition = [2, 2, 5.5, 3.367];

    % First save as a rasterized object
    exportgraphics(gcf, [Multi_func.save_plot 'Flicker' f 'Spectra' f 'Power_spec' f_stim(3:end)  'Hz_' num2str(ylimits(2)) '_' vpo '.png'], ...
                    'Resolution', 600);
    
    % Delete the plot and save the vector graphics axis and heatmap
    surf_h = findobj(ax, "Type", "Surface");
    delete(surf_h);

    % Convert figure to vector graphics and reset the axis
    xlim(ax, xlimits);
    ylim(ax, ylimits);
    caxis(caxis_lts);
    
    set(gcf, 'Renderer', 'painters');
    exportgraphics(gcf, [Multi_func.save_plot 'Flicker' f 'Spectra' f 'Power_spec' f_stim(3:end) 'Hz_' num2str(ylimits(2)) '_' vpo '.pdf'], ...
        'ContentType','vector');
    
end

%% Population calculation for period power spectra for flicker onset, stim, to flicker offset
% Note: use vpo from the population average plots
stats_filename =[Multi_func.save_plot 'Flicker' f 'Spectra' f 'pop_stats_file.txt']; 
stats_file = fopen(stats_filename, 'w');
fclose(stats_file);

stats_t = table();

for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;

    freqs = mean(popul_data.spec_f.n_0.(vpo), 3, 'omitnan');
    
    pop_pow = cellfun(@(f_nr) mean(abs(popul_data.spec_pow.(f_nr).(vpo)), 3, 'omitnan'), ...
                fieldnames(popul_data.spec_pow)', 'UniformOutput', false);

    pop_pow = cat(3, pop_pow{:});
    
    % Calculate the power between each period
    flicker_onset_idxs = find(0.250 < interp_time & 0.750 > interp_time );
    stim_idxs = find(1.250 < interp_time & 1.750 > interp_time );
    flicker_offset_idxs = find(2.250 < interp_time & 2.750 > interp_time );

    % Find the two points that surround 8 Hz
    flick_freq_idxs = find(freqs < 8, 1, 'first');
    flick_freq_idxs = [flick_freq_idxs, flick_freq_idxs - 1]; % Note: frequencies are in decreaseing order

    %TODO switch to narrow band
    %flick_freq_idxs = find(7.5 < freqs & 8.5 > freqs);

    flicker_onset_pow = mean(mean(pop_pow(flick_freq_idxs, flicker_onset_idxs, :), 1), 2);
    stim_pow = mean(mean(pop_pow(flick_freq_idxs, stim_idxs, :), 1), 2);
    flicker_offset_pow = mean(mean(pop_pow(flick_freq_idxs, flicker_offset_idxs, :), 1), 2);

    flicker_onset_pow = squeeze(flicker_onset_pow);
    stim_pow = squeeze(stim_pow);
    flicker_offset_pow = squeeze(flicker_offset_pow);
        
    % Test significnnce of the different periods
    [kw_pval, ~, kw_stats] = kruskalwallis([flicker_onset_pow, stim_pow, flicker_offset_pow], {'pre-stim', 'stim', 'post-stim'});
    close(gcf);
    [kw_c, kw_2, kw_3, mul_h] = multcompare(kw_stats, 'CriticalValueType', 'dunn-sidak');
    close(gcf);
    [sr_pval, ~, sr_stats] = signrank(flicker_onset_pow, stim_pow);
    stats_file = fopen(stats_filename, 'a');
    
    fprintf(stats_file, '-- %s Hz --\n', num2str(f_stim(3:end)));
    fprintf(stats_file, 'Krusk-Wallis test \n');
    fprintf(stats_file, 'KW stat: %s\n', num2str(kw_stats.meanranks));
    fprintf(stats_file, 'KW P-Value: %s\n\n', num2str(kw_pval));
    
    fprintf(stats_file, 'Dunn-sidak test \n');
    fprintf(stats_file, strjoin(mul_h'));
    fprintf(stats_file, '\n');

    for i=1:size(kw_c, 1)
        fprintf(stats_file, '%d\t%d\t%d\n', kw_c(i, 1), kw_c(i, 2), kw_c(i, end));
    end
    fprintf(stats_file, '\n');

    fprintf(stats_file, 'sign-rank test \n');
    fprintf(stats_file, 'SR stat: %s\n', num2str(sr_stats.signedrank));
    fprintf(stats_file, 'SR P-Value: %s\n', num2str(sr_pval));

    fprintf(stats_file, '\n\n');
    fclose(stats_file);

    % Prepare data for violins
    data_pow = [flicker_onset_pow, stim_pow, flicker_offset_pow];
    labels = [repmat({'Flic Onset'}, 1, length(flicker_onset_pow)), ...
            repmat({'Stim'}, 1, length(stim_pow)), ...
            repmat({'Flic Offset'}, 1, length(flicker_offset_pow))  ];
    
    % Create violin plot
    figure('Renderer', 'Painters');
    
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data_pow, labels, 'GroupOrder', {'Flic Onset', 'Stim', 'Flic Offset'}, ...
                ViolinOpts);
    Multi_func.set_default_axis(gca);
    violins(1).ViolinColor = {Multi_func.base_color};
    violins(2).ViolinColor = {Multi_func.stim_color};
    violins(3).ViolinColor = {Multi_func.post_color};
    hold on;
    
    % Plot the individual value points
    plot([1, 2, 3], [flicker_onset_pow, stim_pow, flicker_offset_pow], '-k')
    f_stim
    title([f_stim ' Hz'], 'Interpreter', 'none');

    % Save figure
    saveas(gcf, [Multi_func.save_plot 'Flicker' f 'Spectra' f 'Power_violin_' f_stim(3:end) 'Hz.png']);
    saveas(gcf, [Multi_func.save_plot 'Flicker' f 'Spectra' f 'Power_violin_' f_stim(3:end) 'Hz.pdf']);
    close(gcf);
end

%% Stats for single-Cell level power spectra
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);

    interp_time = popul_data.interp_time;

    % Store the modulation of power spectra for each neuron
    pow_mod = struct; % 2D cell array where the first one is improved, second is worsened, and third is unchanged
    pow_vals = struct;

    for vpo = voices_per_octaves
        vpo = vpo{1};
        pow_mod.(vpo) = {};
        pow_mod.(vpo){1} = {};
        pow_mod.(vpo){2} = {};
        pow_mod.(vpo){3} = {};

        pow_vals.(vpo).onset = struct;
        pow_vals.(vpo).stim = struct;
        pow_vals.(vpo).stats = struct;
    end

    % Loop through each neuron and calculate significance
    for f_nr = fieldnames(popul_data.spec_pow)'
        f_nr = f_nr{1};
        
        for vpo = voices_per_octaves
            vpo = vpo{1};
            
            nr_pow = abs(popul_data.spec_pow.(f_nr).(vpo));
            freqs = mean(popul_data.spec_f.n_0.(vpo), 3, 'omitnan');

            % TODO swtich to narrow band
            %flick_freq_idxs = find(7.5 < freqs & 8.5 > freqs);
            flick_freq_idxs = find(freqs < 8, 1, 'first');
            flick_freq_idxs = [flick_freq_idxs, flick_freq_idxs - 1]; % Note: frequencies are in decreaseing order

            % Calculate each trial's power for each period
            flick_onset_idxs = find(interp_time > 0.250 & interp_time < 0.750);
            stim_idxs = find(interp_time > 1.250 & interp_time < 1.750);
            flick_offset_idxs = find(interp_time > 2.250 & interp_time < 2.750);

            flicker_onset_pow = mean(mean(nr_pow(flick_freq_idxs, flicker_onset_idxs, :), 2), 1);
            stim_pow = mean(mean(nr_pow(flick_freq_idxs, stim_idxs, :), 2), 1);
            flicker_offset_pow = mean(mean(nr_pow(flick_freq_idxs, flicker_offset_idxs, :), 2), 1);
            
            flicker_onset_pow = squeeze(flicker_onset_pow);
            stim_pow = squeeze(stim_pow);
            flicker_offset_pow = squeeze(flicker_offset_pow);

            % Test for significance
            [sr_p, ~, sr_stats] = signrank(flicker_onset_pow, stim_pow);

            % Store the power values
            pow_vals.(vpo).onset.(f_nr) = flicker_onset_pow(:);
            pow_vals.(vpo).stim.(f_nr) = stim_pow(:);

            pow_vals.(vpo).stats.(f_nr) = sr_p(:);

            if sr_p < 0.05
                if mean(stim_pow -  flicker_onset_pow) > 0
                    pow_mod.(vpo){1}{end + 1} = f_nr;
                else
                    pow_mod.(vpo){2}{end + 1} = f_nr;
                end
            else
                pow_mod.(vpo){3}{end + 1} = f_nr;
            end
        end
    end

    % Add to the data structure
    data.(f_stim).pow_mod = pow_mod;
    data.(f_stim).pow_vals = pow_vals;
end

%% Create table indicating the power modulation for each stimulation frequency

stats_t = table();
vpo = 'vpo_24';
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
   
    stats_t([f_stim], 'Pow Imprv') = {length(data.(f_stim).pow_mod.(vpo){1})};
    stats_t([f_stim], 'Pow Wors') = {length(data.(f_stim).pow_mod.(vpo){2})};
    stats_t([f_stim], 'Pow Unch') = {length(data.(f_stim).pow_mod.(vpo){3})};
end
disp(stats_t);
writetable(stats_t, [Multi_func.save_plot 'Flicker' f 'Spectra' f 'power_modulation_' vpo '.csv'], 'WriteRowNames', true);

%% Calculate the wavelet coherence
freqLimits = [0 200];
Fs = 500;
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    
    flicker_raster = popul_data.flicker_raster;

    % Loop through each neuron
    for f_nr = fieldnames(popul_data.raw_vm)'
        f_nr = f_nr{1};

        % Store the power spectra for the current neuron
        nr_spec_wcoh = struct;
        nr_spec_f = struct;
        for vpo = voices_per_octaves
            vpo = vpo{1};
            nr_spec_wcoh.(vpo) = [];
            nr_spec_f.(vpo) = [];
        end
        
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
            
            for vpo = voices_per_octaves
                vpo = vpo{1};
                fb = cwtfilterbank(SignalLength=length(tr_norm_vm), ...
                            SamplingFrequency=Fs, ...
                            FrequencyLimits=freqLimits);
                
                [wcoh, ~, frq] = wcoherence(flicker_raster, tr_norm_vm, Fs, ...
                                'VoicesPerOctave', str2num(vpo(5:end))); % TODO change to 24 for higher resolution
                
                % Concatenate all of the power data
                nr_spec_wcoh.(vpo)(:, :, end + 1) = wcoh;
                nr_spec_f.(vpo)(:, :, end + 1) = frq;
            end
        end

        % Remove empty element
        for vpo = voices_per_octaves
            vpo = vpo{1};
            nr_spec_wcoh.(vpo)(:, :, 1) = [];
            nr_spec_f.(vpo)(:, :, 1) = [];
        end

        % Add to full structure
        data.(f_stim).spec_wcoh.(f_nr) = nr_spec_wcoh;
        data.(f_stim).spec_coh_f.(f_nr) = nr_spec_f;
    end
end

%% Average and plot coherence spectra
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;

    for vpo = voices_per_octaves
        vpo = vpo{1};
        cur_coh = cellfun(@(f_nr) mean(abs(popul_data.spec_wcoh.(f_nr).(vpo)), 3, 'omitnan'), ...
                    fieldnames(popul_data.spec_wcoh)', 'UniformOutput', false);
        cur_coh = cat(3, cur_coh{:});

        figure;
        surface(interp_time, ... 
                mean(popul_data.spec_coh_f.n_0.(vpo), 3, 'omitnan'), ...
                mean(cur_coh, 3, 'omitnan'), 'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
        colormap(jet*.8);
        colorbar;
        ylim([0 20]);
        title([f_stim ' Hz ' vpo], 'Interpreter', 'none');
    end
end

%% Population calculate the period coherences for flicker onset, stim, to flicker offset
stats_filename =[Multi_func.save_plot 'Flicker' f 'Coherence' f 'pop_stats_file.txt']; 
stats_file = fopen(stats_filename, 'w');
fclose(stats_file);
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;

    freqs = mean(popul_data.spec_coh_f.n_0, 3, 'omitnan');
    
    % calculate the population coherence 
    pop_coh = cellfun(@(f_nr) mean(abs(popul_data.spec_wcoh.(f_nr)), 3, 'omitnan'), ...
                fieldnames(popul_data.spec_wcoh)', 'UniformOutput', false);

    pop_coh = cat(3, pop_coh{:});
    
    % Calculate the power between each period
    flicker_onset_idxs = find(0.250 < interp_time & 0.750 > interp_time );
    stim_idxs = find(1.250 < interp_time & 1.750 > interp_time );
    flicker_offset_idxs = find(2.250 < interp_time & 2.750 > interp_time );

    flick_freq_idxs = find(7.5 < freqs & 8.5 > freqs);

    flicker_onset_coh = mean(mean(pop_coh(flick_freq_idxs, flicker_onset_idxs, :), 1), 2);
    stim_coh = mean(mean(pop_coh(flick_freq_idxs, stim_idxs, :), 1), 2);
    flicker_offset_coh = mean(mean(pop_coh(flick_freq_idxs, flicker_offset_idxs, :), 1), 2);

    flicker_onset_coh = squeeze(flicker_onset_coh);
    stim_coh = squeeze(stim_coh);
    flicker_offset_coh = squeeze(flicker_offset_coh);
    

    data_coh = [flicker_onset_coh, stim_coh, flicker_offset_coh];
    labels = [repmat({'Flic Onset'}, 1, length(flicker_onset_coh)), ...
            repmat({'Stim'}, 1, length(stim_coh)), ...
            repmat({'Flic Offset'}, 1, length(flicker_offset_coh))  ];
    
    % Create violin plot
    figure;
    
    ViolinOpts = Multi_func.get_default_violin();
    violins = violinplot(data_coh, labels, 'GroupOrder', {'Flic Onset', 'Stim', 'Flic Offset'}, ...
                ViolinOpts);
    
    violins(1).ViolinColor = {Multi_func.base_color};
    violins(2).ViolinColor = {Multi_func.stim_color};
    violins(3).ViolinColor = {Multi_func.post_color};
    hold on;
    
    % Plot the individual value points
    plot([1, 2, 3], [flicker_onset_coh, stim_coh, flicker_offset_coh], '-k')

    title([f_stim ' Hz'], 'Interpreter', 'none');
    
    % Test significnnce of the different periods
    [kw_pval, ~, kw_stats] = kruskalwallis([flicker_onset_coh, stim_coh, flicker_offset_coh]);
    kw_c = multcompare(kw_stats, 'CriticalValueType', 'dunn-sidak');
    [sr_pval, ~, sr_stats] = signrank(flicker_onset_coh, stim_coh);
    stats_file = fopen(stats_filename, 'a');
    
    fprintf(stats_file, '-- %s Hz --\n', num2str(f_stim(3:end)));
    fprintf(stats_file, 'Krusk-Wallis test \n');
    fprintf(stats_file, 'KW stat: %s\n', num2str(kw_stats.meanranks));
    fprintf(stats_file, 'KW P-Value: %s\n\n', num2str(kw_pval));
    
    fprintf(stats_file, 'Dunn-sidak test \n');
    for i=1:size(kw_c, 1)
        fprintf(stats_file, '%d\t%d\t%d\n', kw_c(i, 1),kw_c(i, 2), kw_c(i, end));
    end
    fprintf(stats_file, '\n');

    fprintf(stats_file, 'sign-rank test \n');
    fprintf(stats_file, 'SR stat: %s\n', num2str(sr_stats.signedrank));
    fprintf(stats_file, 'SR P-Value: %s\n', num2str(sr_pval));

    fprintf(stats_file, '\n\n');
    fclose(stats_file);
end

%% Stats for single-cell coherence
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;

    % Store the modulation of power spectra for each neuron
    coh_mod = struct; % 2D cell array where the first one is improved, second is worsened, and third is unchanged
    coh_vals = struct;
    for vpo = voices_per_octaves
        vpo = vpo{1};
        coh_mod.(vpo) = {};
        coh_mod.(vpo){1} = {};
        coh_mod.(vpo){2} = {};
        coh_mod.(vpo){3} = {};

        coh_vals.(vpo).onset = struct;
        coh_vals.(vpo).stim = struct;
        coh_vals.(vpo).stats = struct;
    end
    
    % Loop through each neuron and calculate significance
    for f_nr = fieldnames(popul_data.spec_wcoh)'
        f_nr = f_nr{1};
        
        % Loop through each voices
        for vpo = voices_per_octaves
            vpo = vpo{1};
            freqs = mean(popul_data.spec_coh_f.(f_nr).(vpo), 3, 'omitnan');
        
            nr_coh = abs(popul_data.spec_wcoh.(f_nr).(vpo));
            
            % Find the two points that surround 8 Hz
            flick_freq_idxs = find(freqs < 8, 1, 'first');
            flick_freq_idxs = [flick_freq_idxs, flick_freq_idxs - 1]; % Note: frequencies are in decreaseing order
                
            % OLD way and do not want to use
            %flick_freq_idxs = find(7.5 < freqs < 8.5);


            % Calculate each trial's coher for each period
            flicker_onset_idxs = find(interp_time > 0.250 & interp_time < 0.750);
            stim_idxs = find(interp_time > 1.250 & interp_time < 1.750);
            flicker_offset_idxs = find(interp_time > 2.250 & interp_time < 2.750);

            flicker_onset_coh = mean(mean(nr_coh(flick_freq_idxs, flicker_onset_idxs, :), 2), 1);
            stim_coh = mean(mean(nr_coh(flick_freq_idxs, stim_idxs, :), 2), 1);
            flicker_offset_coh = mean(mean(nr_coh(flick_freq_idxs, flicker_offset_idxs, :), 2), 1);
            
            flicker_onset_coh = squeeze(flicker_onset_coh);
            stim_coh = squeeze(stim_coh);
            flicker_offset_coh = squeeze(flicker_offset_coh);

            % Test for significance
            [sr_p, ~, sr_stats] = signrank(flicker_onset_coh, stim_coh);
            
            % Store the onset and stim coh values
            coh_vals.(vpo).onset.(f_nr) = flicker_onset_coh(:);
            coh_vals.(vpo).stim.(f_nr) = stim_coh(:);
            
            coh_vals.(vpo).stats.(f_nr) = sr_p(:);

            % DEBUG to show where the 8 Hz band is
            %figure;
            %surface(interp_time, freqs, mean(nr_coh, 3, 'omitnan'), ...
            %    'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            %ylim([7 9]);
            %hold on;
            %yline(freqs(flick_freq_idxs)); % Show the 8 Hz frequency band that is taken
            %title([vpo], 'Interpreter', 'none');

            if sr_p < 0.05
                if mean(stim_coh -  flicker_onset_coh) > 0
                    coh_mod.(vpo){1}{end + 1} = f_nr;
                else
                    coh_mod.(vpo){2}{end + 1} = f_nr;
                end
            else
                coh_mod.(vpo){3}{end + 1} = f_nr;
            end
        end
    end

    % Add to the data structure
    data.(f_stim).coh_mod = coh_mod;
    data.(f_stim).coh_vals = coh_vals;
end

%% Plotting to check the difference each neuron's coherence values with p values
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    interp_time = popul_data.interp_time;
    
    % Loop through each neuron's coherence values
    for f_nr = fieldnames(popul_data.spec_wcoh)'
        f_nr = f_nr{1};
 
        % Store individual neuron's coherence values across different voicesPerOctave
        cur_onset_vals = [];
        cur_stim_vals = [];

        % Loop through each voices
        for vpo = voices_per_octaves
            vpo = vpo{1};
            
            cur_onset_vals(:, end + 1) = popul_data.coh_vals.(vpo).onset.(f_nr);
            cur_stim_vals(:, end + 1) = popul_data.coh_vals.(vpo).stim.(f_nr);
        end

        % Plot the voices per octave data
        figure;
        %tiledlayout(1, 3);
        tiledlayout(1, 2);
        nexttile;
        plot(cur_onset_vals');
        title('Onset Values');

        nexttile;
        plot(cur_stim_vals');
        title('Stim values');

        % Plot with change from onset to stim
        %for i=1:2
        %    nexttile;
        %    plot([1, 2], [cur_onset_vals(:, i), cur_stim_vals(:, i) ]);
        %    delta = mean(cur_stim_vals(:, i) - cur_onset_vals(:, i));
        %    legend([voices_per_octaves{i} ' ' num2str(delta)], 'Interpreter', 'none');
        %    title(['p=' num2str(popul_data.coh_vals.(voices_per_octaves{i}).stats.(f_nr) )]);
        %end
        
        sgtitle([f_stim ' ' f_nr], 'Interpreter', 'none');
    end   
end

%% Create table indicating the coherence modulation for each stimulation frequency

stats_t = table();
for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
   
    stats_t([f_stim], 'Coh Imprv') = {length(popul_data.coh_mod{1})};
    stats_t([f_stim], 'Coh Wors') = {length(popul_data.coh_mod{2})};
    stats_t([f_stim], 'Coh Unch') = {length(popul_data.coh_mod{3})};
end
disp(stats_t);
writetable(stats_t, [Multi_func.save_plot 'Flicker' f 'Coherence' f 'coh_modulation.csv'], 'WriteRowNames', true);

%% Show individual trials with filtered 8 Hz for neuron that improved their coherence during stim
flick_base_tmps = [0 1];
stim_timestamps = [1 2];
gray = [0.5 0.5 0.5];

% Set the voices_per_octave to use
vpo = 'vpo_24';

for f_stim = fieldnames(data)'
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;

    flicker_raster = popul_data.flicker_raster;

    % loop only through the increased and decreased neurons
    for f_mod = 1:2
        % Loop through neuron per category coherence
               
        % For debugging purposes
        nrs = popul_data.coh_mod.(vpo){f_mod};
        if length(nrs) > 1
            nrs = nrs(1);
        end

        for f_nr = popul_data.coh_mod.(vpo){f_mod} %nrs % TODO change to popul_data.coh_mod{f_mod} again
            f_nr = f_nr{1};

            % Normalize trials and filter 8 Hz component
            norm_trs = Multi_func.norm_signals(popul_data.raw_vm.(f_nr));
            trs_8hz = Multi_func.raw_filt(norm_trs, 8, 500);
            avg_8hz_tr = mean(trs_8hz, 2, 'omitnan');
            
            %Spike normalized and baseline subtracted
            sp_amp = popul_data.spike_amp_raster.(f_nr);
            
            sp_amp_trs = mean(sp_amp, 1, 'omitnan');
            sp_amp_trs(isnan(sp_amp_trs)) == 1;
            amp_norm_trs = popul_data.raw_vm.(f_nr)./sp_amp_trs;
            
            avg_vm = mean(amp_norm_trs, 2, 'omitnan');
            sem_vm = std(amp_norm_trs, 0, 2, 'omitnan')./sqrt(size(amp_norm_trs, 2));

            % Calculate the neurons average coherence
            avg_coh = mean(abs(popul_data.spec_wcoh.(f_nr).(vpo)), 3, 'omitnan');
            avg_coh_f = mean(popul_data.spec_coh_f.(f_nr).(vpo), 3, 'omitnan');


            % Calculate the neuron's average power spectra
            get_base_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp < stim_tmstp(1));
            get_stim_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp >= stim_tmstp(1) & ...
                                                       tr_tmstmp <= stim_tmstp(end));

            calc_time_pow = @(trial_spec, time_idxs) mean(trial_spec(:, time_idxs, :), 2);
            
            calc_trial_spec = @(trial_spec, base_pow, stim_pow) (trial_spec - base_pow)./(base_pow + stim_pow);
            
            calc_nr_spec = @(trial_spec, tr_tmstmp, stim_tmstp) mean(calc_trial_spec(trial_spec, ...
                calc_time_pow(trial_spec, get_base_idxs(tr_tmstmp, stim_tmstp)),  ...
                calc_time_pow(trial_spec, get_stim_idxs(tr_tmstmp, stim_tmstp))), 3, 'omitnan');
            avg_pow = calc_nr_spec(abs(popul_data.spec_pow.(f_nr).(vpo)), timeline, flick_base_tmps);
            avg_pow_f = mean(popul_data.spec_f.(f_nr).(vpo), 3, 'omitnan');

            % Plot figure for individual neuron
            figure('Renderer', 'Painters', 'Units', 'centimeters', 'Position', [4 20 21.59 27.94]);
            %figure;
            %tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact', 'Units', 'centimeters', 'InnerPosition', [4, 4, 7.2, 9]);
            tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
            
            % Plot all of the trials
            %nexttile;
            
            nexttile([2, 1]); % Make the trials span two rows

            plot(timeline, norm_trs - 0.5 + repmat(1:size(norm_trs, 2), size(norm_trs, 1), 1), '-k');
            hold on;
            plot(timeline, trs_8hz + repmat(1:size(trs_8hz, 2), size(trs_8hz, 1), 1), '-b');
            hold on;
            xline([0, 3], '--', 'color', gray);
            hold on;
            xline([1, 2], '--', 'color', gray);
            hold on;
            plot(timeline, 1.*flicker_raster + size(norm_trs, 2) + 1, '-b');
            hold on;
            plot(stim_timestamps, 0.5.*[1 1] + size(norm_trs, 2) + 2, 'color', Multi_func.dbs_color);
            Multi_func.set_default_axis(gca);
            title('All Trials');

            % Plot trial average Vm
            nexttile;
            Multi_func.set_default_axis(gca);
            tr_fill = fill([timeline, flip(timeline)], ...
                [avg_vm + sem_vm; flip(avg_vm - sem_vm)]', ...
                [0.5 0.5 0.5]);
            Multi_func.set_fill_properties(tr_fill);
            hold on;
            plot(timeline, avg_vm, '-k');
            hold on;
            plot(timeline, avg_8hz_tr, '-b');
            Multi_func.set_default_axis(gca);

            % Plot the coherence spectrum
            nexttile;
            surface(timeline, ...
                avg_coh_f, avg_coh, ...
                'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
            
            %imagesc('XData', timeline, 'YData', avg_coh_f, 'CData', avg_coh);

            %colormap(jet*.8);
            colormap(Multi_func.red_purple_blue_color);
            colorbar;
            hold on;
            plot(timeline, 10.*flicker_raster + max(avg_coh_f) + 1);
            ax = gca;
            Multi_func.set_default_axis(ax);
            
            ylim([0 20]);
            title('Coherence');
            
            % Save the X and Y axes limits
            xlimits = xlim(ax);
            ylimits = ylim(ax);
            caxis_lts = caxis;

           % % Plot the power spectra
           % nexttile;
           % surface(timeline, ...
           %     avg_pow_f, avg_pow, ...
           %     'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'edgecolor', 'none');
           % 
           % colormap(jet*.8);
           % colorbar;
           % hold on;
           % plot(timeline, 10.*flicker_raster + max(avg_coh_f) + 1);
           % Multi_func.set_default_axis(gca);
           % 
           % ylim([0 20]);
           % title('Power');

            sgtitle([popul_data.nr_name.(f_nr) ' mod: ' mod_str{f_mod} ' stim: ' f_stim], ...
                'Interpreter', 'none');
            
            exportgraphics(gcf, [savefig_path 'Coherence' f 'Mod' f ...
                'Mod_' mod_str{f_mod} '_' f_stim '_' f_nr '_' vpo '.png'], 'Resolution', 600);
            
            %TODO need to take this back
            % Delete the surface plot and reset the axis
            surf_h = findobj(ax, "Type", "Surface");
            delete(surf_h);

            xlim(ax, xlimits);
            ylim(ax, ylimits);
            caxis(caxis_lts);

            exportgraphics(gcf, [savefig_path 'Coherence' f 'Mod' f ...
                'Mod_' mod_str{f_mod} '_' f_stim '_' f_nr '_' vpo '.pdf'], 'ContentType', 'vector');
        end
    end
end


%% Plotting the signal to demonstrate the concept of wavelet coherence
figure;

% Plot the flicker raster
plot(timeline, 5+Multi_func.norm_signals(flicker_raster'));
hold on;

% Plot the trial average
plot(timeline, 4+Multi_func.norm_signals(avg_vm));

% Filter the flicker raster
rast_8Hz = Multi_func.raw_filt(flicker_raster', 8, 500);
plot(timeline, 3+Multi_func.norm_signals(rast_8Hz));

% Plot the 8 Hz trial average
plot(timeline, 2+Multi_func.norm_signals(avg_8hz_tr));
hold on;

% Plot the absolute value of the raster
plot(timeline, 1+Multi_func.norm_signals(abs(rast_8Hz)));
hold on;

% Plot the absolute value of the average trial 8 Hz
plot(timeline, Multi_func.norm_signals(abs(avg_8hz_tr)));
hold on;
