%% This file requires flicker_plv.m and flicker_spectra.m to show the single-cell stats
savefig_path = [Multi_func.save_plot]; 


%% Plot each neuron's PLV distribution, power, and coherence with 8 Hz
% violin plots for all the properties and 
flick_base_tmps = [0 1];
stim_timestamps = [1 2];
gray = [0.5 0.5 0.5];
freqs = Multi_func.entr_freqs;
set_xlim = [0, 10];
polar_edges = linspace(0, 2*pi, 24);

% Interpolate the frequency space        
freq_lin = linspace(Multi_func.entr_freqs(1), Multi_func.entr_freqs(end), 200);
for f_stim = fieldnames(data)' %{'f_140'} %
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;

    % Loop through all of the neurons
    for f_nr = fieldnames(popul_data.nr_name)' %{'n_10'} % {'n_12'}%
        f_nr = f_nr{1};
        
        fig = figure('Position', [400 400 1800 1200]);
        
        t = tiledlayout(4, 2, 'Padding', 'tight');

        % Plot the power
        nexttile();
        get_base_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp < stim_tmstp(1));
        get_stim_idxs = @(tr_tmstmp, stim_tmstp) find(tr_tmstmp >= stim_tmstp(1) & ...
                                                   tr_tmstmp <= stim_tmstp(end));
        calc_time_pow = @(trial_spec, time_idxs) mean(trial_spec(:, time_idxs, :), 2);
        
        calc_trial_spec = @(trial_spec, base_pow, stim_pow) (trial_spec - base_pow)./(base_pow + stim_pow);
        
        calc_nr_spec = @(trial_spec, tr_tmstmp, stim_tmstp) mean(calc_trial_spec(trial_spec, ...
            calc_time_pow(trial_spec, get_base_idxs(tr_tmstmp, stim_tmstp)),  ...
            calc_time_pow(trial_spec, get_stim_idxs(tr_tmstmp, stim_tmstp))), 3, 'omitnan');
        avg_pow = calc_nr_spec(abs(popul_data.spec_pow.(f_nr).vpo_24), timeline, flick_base_tmps);
        avg_f = mean(popul_data.spec_f.(f_nr).vpo_24, 3, 'omitnan');
        
        % Interpolate frequency points
        pow_interp = interp1(avg_f, avg_pow, freq_lin);

        %surface(timeline, ...
        %       freq_lin, pow_interp, ...
        %        'CDataMapping', 'scaled', 'FaceColor', 'texturemap', 'EdgeColor', 'none');
        imagesc(timeline, freq_lin, pow_interp);
        ylim([0, 20]);
        axis xy;
        xlabel('Time');
        ylabel('Frequency (Hz)');
        %ylim([0, 20]);
        %pcolor(timeline, avg_f, avg_pow);
        
        % Check if neuron's power spec is significant
        sig = [];
        if ismember(f_nr, popul_data.pow_mod.vpo_24{1})
            sig = 'Inc';
        elseif ismember(f_nr, popul_data.pow_mod.vpo_24{2})
            sig = 'Dec';
        else
            sig = '';
        end
        title(['Power Spectra ' sig ' ' num2str(popul_data.pow_vals.vpo_24.stats.(f_nr))]);

        % Power violin
        nexttile;
        pow_f = mean(popul_data.spec_f.n_0.(vpo), 3, 'omitnan');
        % TODO change the freqs in this block to pow_f

        flick_freq_idxs = find(pow_f < 8, 1, 'first');
        flick_freq_idxs = [flick_freq_idxs, flick_freq_idxs - 1]; % Note: frequencies are in decreaseing order

        flick_onset_idxs = find(timeline > 0.250 & timeline < 0.750);
        stim_idxs = find(timeline > 1.250 & timeline < 1.750);
        flick_offset_idxs = find(timeline > 2.250 & timeline < 2.750);

        nr_pow = abs(popul_data.spec_pow.(f_nr).(vpo));

        flicker_onset_pow = mean(mean(nr_pow(flick_freq_idxs, flick_onset_idxs, :), 2), 1);
        stim_pow = mean(mean(nr_pow(flick_freq_idxs, stim_idxs, :), 2), 1);
        flicker_offset_pow = mean(mean(nr_pow(flick_freq_idxs, flick_offset_idxs, :), 2), 1);

        flicker_onset_pow = squeeze(flicker_onset_pow);
        stim_pow = squeeze(stim_pow);
        flicker_offset_pow = squeeze(flicker_offset_pow);

        data_pow = [flicker_onset_pow, stim_pow, flicker_offset_pow];
        labels = [repmat({'Pre-stim'}, 1, length(flicker_onset_pow)), ...
            repmat({'Stim'}, 1, length(stim_pow)), ...
            repmat({'Post-stim'}, 1, length(flicker_offset_pow))  ];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data_pow, labels, 'GroupOrder', {'Pre-stim', 'Stim', 'Post-stim'}, ...
            ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        violins(3).ViolinColor = {Multi_func.post_color};
        hold on;

        % Plot the lines between correspond pow values
        plot([1, 2, 3], [flicker_onset_pow, stim_pow, flicker_offset_pow], '-k');
        hold on;
        Multi_func.set_default_axis(gca);

        title('Power');

        % Plot the coherence
        nexttile;
        avg_coh = mean(abs(popul_data.spec_wcoh.(f_nr).vpo_24), 3, 'omitnan');
        avg_coh_f = mean(popul_data.spec_coh_f.(f_nr).vpo_24, 3, 'omitnan');
        
        % Interpolate frequency points
        coh_interp = interp1(avg_coh_f, avg_coh, freq_lin);
            
        imagesc(timeline, freq_lin, coh_interp);
        ylim([0, 20]);
        axis xy;
        xlabel('Time');
        ylabel('Frequency (Hz)');

        % Check if neuron's coherence is significant
        sig = [];
        if ismember(f_nr, popul_data.coh_mod.vpo_24{1})
            sig = 'Inc';
        elseif ismember(f_nr, popul_data.coh_mod.vpo_24{2})
            sig = 'Dec';
        else
            sig = '';
        end
        title(['Coherence ' sig ' ' num2str(popul_data.coh_vals.vpo_24.stats.(f_nr))]);

       
        % Coherence violin
        nexttile;
        % TODO freqs here needs to be changed for the coherence calculated
        % values for frequencies as these are different from the straight
        % linear version
        coh_f = mean(popul_data.spec_coh_f.n_0.(vpo), 3, 'omitnan');
        
        flick_freq_idxs = find(coh_f < 8, 1, 'first');
        flick_freq_idxs = [flick_freq_idxs, flick_freq_idxs - 1]; % Note: frequencies are in decreaseing order

        flick_onset_idxs = find(timeline > 0.250 & timeline < 0.750);
        stim_idxs = find(timeline > 1.250 & timeline < 1.750);
        flick_offset_idxs = find(timeline > 2.250 & timeline < 2.750);

        nr_coh = abs(popul_data.spec_wcoh.(f_nr).(vpo));

        flicker_onset_coh = mean(mean(nr_coh(flick_freq_idxs, flick_onset_idxs, :), 2), 1);
        stim_coh = mean(mean(nr_coh(flick_freq_idxs, stim_idxs, :), 2), 1);
        flicker_offset_coh = mean(mean(nr_coh(flick_freq_idxs, flick_offset_idxs, :), 2), 1);

        flicker_onset_coh = squeeze(flicker_onset_coh);
        stim_coh = squeeze(stim_coh);
        flicker_offset_coh = squeeze(flicker_offset_coh);

        data_coh = [flicker_onset_coh, stim_coh, flicker_offset_coh];
        labels = [repmat({'Pre-stim'}, 1, length(flicker_onset_coh)), ...
            repmat({'Stim'}, 1, length(stim_coh)), ...
            repmat({'Post-stim'}, 1, length(flicker_offset_coh))  ];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data_coh, labels, 'GroupOrder', {'Pre-stim', 'Stim', 'Post-stim'}, ...
            ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        violins(3).ViolinColor = {Multi_func.post_color};
        hold on;

        % Plot the lines between correspond pow values
        plot([1, 2, 3], [flicker_onset_coh, stim_coh, flicker_offset_coh], '-k');
        hold on;
        Multi_func.set_default_axis(gca);

        title('Coherence');

        %-- Plot the PLV
        nexttile();
        % Plotting the flicker baseline
        flicker_base_mean = mean(popul_data.nr_plvs.(f_nr).base_PLV2s, 2, 'omitnan');
        flicker_base_std = std(popul_data.nr_plvs.(f_nr).base_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).base_PLV2s, 2);
        flicker_base_sem = flicker_base_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [flicker_base_mean' + flicker_base_sem', ...
            flip(flicker_base_mean' - flicker_base_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, flicker_base_mean, 'color', Multi_func.base_color, DisplayName='Base');
        hold on;

        % Plotting the estim
        estim_mean = mean(popul_data.nr_plvs.(f_nr).estim_PLV2s, 2, 'omitnan');
        estim_std = std(popul_data.nr_plvs.(f_nr).estim_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).estim_PLV2s, 2);
        estim_sem = estim_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [estim_mean' + estim_sem', ...
            flip(estim_mean' - estim_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, estim_mean, 'color', Multi_func.stim_color, DisplayName='Estim');
        hold on;

        % Plotting the flicker offset
        flicker_offset_mean = mean(popul_data.nr_plvs.(f_nr).offset_PLV2s, 2, 'omitnan');
        flicker_offset_std = std(popul_data.nr_plvs.(f_nr).offset_PLV2s, [], 2, 'omitnan');
        num_flick_base = size(popul_data.nr_plvs.(f_nr).offset_PLV2s, 2);
        flicker_offset_sem = flicker_offset_std./sqrt(num_flick_base);

        fill_h = fill([freqs, flip(freqs)], [flicker_offset_mean' + flicker_offset_sem', ...
            flip(flicker_offset_mean' - flicker_offset_sem')], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);

        Multi_func.set_default_axis(gca);
        hold on;
        plot(freqs, flicker_offset_mean, 'color', Multi_func.post_color, DisplayName='Offset');
        hold on;

        % Check if neuron's PLV is significant
        sig = [];
        if ismember(f_nr, popul_data.plv_mod{1})
            sig = 'Inc';
        elseif ismember(f_nr, popul_data.plv_mod{2})
            sig = 'Dec';
        else
            sig = '';
        end
        title(['PLV ' sig ' ' num2str(popul_data.plv_stats.(f_nr))]);
            
        %TODO need to doube check wat the stats are being present here,
        % Check how they are being saved

        %ax = gca;
        %set(ax,'Xscale','log');
        xlabel('Frequency (Hz)');
        ylabel('PLV');
        xlim(set_xlim);

        % Show only plots that have DisplayNames
        lines = findall(gca, 'Type', 'Line');
        lines = lines(~cellfun(@isempty, get(lines, 'DisplayName')));
        legend(lines, 'Location', 'east');
        %-- End of PLV plotting

        % PLV violin
        nexttile;

        % Get PLVs from each period
        flicker_onset_plvs = popul_data.nr_plvs.(f_nr).base_PLV2s(Multi_func.entr_freqs == 8, :);
        stim_plvs = popul_data.nr_plvs.(f_nr).estim_PLV2s(Multi_func.entr_freqs == 8, :);
        flicker_offset_plvs = popul_data.nr_plvs.(f_nr).offset_PLV2s(Multi_func.entr_freqs == 8, :);
        
        % Make the violin plot
        data_plvs = [flicker_onset_plvs, stim_plvs, flicker_offset_plvs];
        labels = [repmat({'Pre-stim'}, 1, length(flicker_onset_plvs)), ...
            repmat({'Stim'}, 1, length(stim_plvs)), ...
            repmat({'Post-stim'}, 1, length(flicker_offset_plvs))  ];

        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data_plvs, labels, 'GroupOrder', {'Pre-stim', 'Stim', 'Post-stim'}, ...
            ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        violins(3).ViolinColor = {Multi_func.post_color};
        hold on;

        % Plot the lines between correspond PLV values
        plot([1, 2, 3], [flicker_onset_plvs', stim_plvs', flicker_offset_plvs'], '-k');
        hold on;
        Multi_func.set_default_axis(gca);
        
        title('PLV');

        num_trials = length(flicker_onset_plvs);

        % Plot the flicker phases for each period
        nexttile([1 2]);
        delete(gca);
        t_nest = tiledlayout(t, 1, 3);
        t_nest.Layout.Tile = 7;
        t_nest.Layout.TileSpan = [1 2];

        nexttile(t_nest);
        base_vecs = popul_data.nr_plvs.(f_nr).base_norm_vecs;
        base_vecs = cellfun(@(ar) ar(find(Multi_func.entr_freqs == 8), :), base_vecs, 'UniformOutput', false);
        base_vecs = [base_vecs{:}];

        polarhistogram(angle(base_vecs), polar_edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.base_color);
        title('Flicker Base Phases');

        nexttile(t_nest);
        estim_vecs = popul_data.nr_plvs.(f_nr).estim_norm_vecs;
        estim_vecs = cellfun(@(ar) ar(find(Multi_func.entr_freqs == 8), :), estim_vecs, 'UniformOutput', false);
        estim_vecs = [estim_vecs{:}];

        % Circular statistical testing
        [p_vals, ~] = circ_wwtest(angle(base_vecs), angle(estim_vecs));

        polarhistogram(angle(estim_vecs), polar_edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.stim_color);
        title(['Flicker Estim Phases: ' num2str(p_vals)]);

        nexttile(t_nest);
        offset_vecs = popul_data.nr_plvs.(f_nr).offset_norm_vecs;
        offset_vecs = cellfun(@(ar) ar(find(Multi_func.entr_freqs == 8), :), offset_vecs, 'UniformOutput', false);
        offset_vecs = [offset_vecs{:}];

        polarhistogram(angle(offset_vecs), polar_edges, 'Normalization', 'probability', 'DisplayStyle', 'stairs', 'EdgeColor', Multi_func.post_color);
        title('Flicker Offset Phases');

        sgtitle([f_stim 'Hz ' f_nr ' ' popul_data.nr_name.(f_nr) ' Trial#: ' num2str(num_trials)], 'Interpreter', 'none');

        fontsize(fig, 10, "points");
        %set(findall(fig,'-property','FontSize'), 'FontSize', 20);

        exportgraphics(fig, [savefig_path 'Flicker' f 'Neuronwise' f '8Hz_plots' f ...
           f_stim(3:end) '_' popul_data.nr_name.(f_nr) '_.png'], 'Resolution', 600);
   
    end
end

%% Create a table displaying all of the info for single-cells
% Have the relabels ready
mouse_rename = struct();
mice_names = fieldnames(Multi_func.mouse_color)';
for f_i = 1:length(mice_names)
    mouse_rename.(mice_names{f_i}) = ['Mouse ' num2str(f_i)];
end

flicker_nr_t = table();

row_num = 1;
for f_stim = fieldnames(data)' %{'f_140'} %
    f_stim = f_stim{1};
    popul_data = data.(f_stim);

    for f_nr = fieldnames(popul_data.spec_pow)' %{'n_12'}%
        f_nr = f_nr{1};

        % Parse out recording info
        name_parts = strsplit(popul_data.nr_name.(f_nr), "_");
        mouse_name = name_parts{1};

        rec_date = datetime(name_parts{4}, 'InputFormat', 'yyyyMMdd');
        
        % Get neurons with the same name
        fn = fieldnames(popul_data.nr_name);
        vals = struct2cell(popul_data.nr_name);
        idx = contains(vals, mouse_name);
        same_m_names = vals(idx);
        
        % Find the earliest recording date for this neuron
        splitNames = cellfun(@(x) strsplit(x, '_'), same_m_names, 'UniformOutput', false);
        dateStrs = cellfun(@(x) x{4}, splitNames, 'UniformOutput', false);
        dates = datetime(dateStrs, 'InputFormat', 'yyyyMMdd');
        earliestDate = min(dates);

        % Parse out the current amplitude
        amp_str = popul_data.stim_amp.(f_nr);
        start_idx = find(amp_str == ':') + 1;
        end_idx = find(amp_str == 'u', 1, 'last') - 1;
        
        amp_str = amp_str(start_idx:end_idx);

        % TODO re-export the data from Python to MATLAB
        % Make sure that the save matfile has essentially the same data and
        % its just the current amplitude that is being added

        % Add parameters to table
        flicker_nr_t{row_num, 'Mouse_ID'} = string(mouse_rename.(['m' mouse_name]));
        flicker_nr_t{row_num, 'Stim Frequency'} = str2num(f_stim(3:end));
        flicker_nr_t{row_num, 'Current'} = str2num(amp_str);
        
        % Enter statistical significance of PLV
        sig = [];
        if ismember(f_nr, popul_data.plv_mod{1})
            sig = '1';
        elseif ismember(f_nr, popul_data.plv_mod{2})
            sig = '-1';
        else
            sig = '0';
        end
        flicker_nr_t{row_num, 'PLV Stats'} = popul_data.plv_stats.(f_nr);
        flicker_nr_t{row_num, 'PLV Sign'} = str2num(sig);

        % Enter coherrence significance
        sig = [];
        if ismember(f_nr, popul_data.coh_mod.vpo_24{1})
            sig = '1';
        elseif ismember(f_nr, popul_data.coh_mod.vpo_24{2})
            sig = '-1';
        else
            sig = '0';
        end
        flicker_nr_t{row_num, 'Coherence stats'} = popul_data.coh_vals.vpo_24.stats.(f_nr);
        flicker_nr_t{row_num, 'Coherence Sign'} = str2num(sig);

        % Enter power significance
        sig = [];
        if ismember(f_nr, popul_data.pow_mod.vpo_24{1})
            sig = '1';
        elseif ismember(f_nr, popul_data.pow_mod.vpo_24{2})
            sig = '-1';
        else
            sig = '0';
        end
        flicker_nr_t{row_num, 'Power stats'} = popul_data.pow_vals.vpo_24.stats.(f_nr);
        flicker_nr_t{row_num, 'Power Sign'} = str2num(sig);

        % Add the recording day
        flicker_nr_t{row_num, 'Recording Day'} = days(rec_date - earliestDate);

        row_num = row_num + 1;
    end
end

%DEBUG
flicker_nr_t
writetable(flicker_nr_t, [savefig_path f 'Flicker' f 'Neuronwise' f 'Neurons_Flicker_Data.csv'], 'WriteRowNames', true);

%% Create a table grouped my mouse to show significances across categories
flicker_mouse_t = table();
row_num = 1;
for f_stim = fieldnames(data)' %{'f_140'} %
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    
    cond_idxs = flicker_nr_t.("Stim Frequency") == str2num(f_stim(3:end));

    cond_mice = unique(flicker_nr_t.Mouse_ID(cond_idxs));
    for mouse = cond_mice'
        mice_idxs = cond_idxs & (flicker_nr_t.Mouse_ID == mouse);

        cur_amp = unique(flicker_nr_t.Current(mice_idxs));
        for amp = cur_amp'
            amp_idxs = mice_idxs & (flicker_nr_t.Current == amp);

            % Save these group of neuron's identifiables
            flicker_mouse_t{row_num, 'Mouse'} = mouse;
            flicker_mouse_t{row_num, 'Stim'} = str2num(f_stim(3:end));
            flicker_mouse_t{row_num, 'Current'} = amp;

            % Count number of significant neurons for PLV
            sig_idx = (abs(flicker_nr_t.("PLV Sign")(amp_idxs)) > 0);
            flicker_mouse_t{row_num, 'PLV Sig'} = sum(sig_idx);
            sig_idx = (flicker_nr_t.("PLV Sign")(amp_idxs) == 0);
            flicker_mouse_t{row_num, 'Not PLV Sig'} = sum(sig_idx);

            % Count number of significant neurons for Coherence
            sig_idx = (abs(flicker_nr_t.("Coherence Sign")(amp_idxs)) > 0);
            flicker_mouse_t{row_num, 'Coherence Sig'} = sum(sig_idx);
            sig_idx = (flicker_nr_t.("Coherence Sign")(amp_idxs) == 0);
            flicker_mouse_t{row_num, 'Not Coherence Sig'} = sum(sig_idx);

            % Count number of significant neurons for Power
            sig_idx = (abs(flicker_nr_t.("Power Sign")(amp_idxs)) > 0);
            flicker_mouse_t{row_num, 'Power Sig'} = sum(sig_idx);
            sig_idx = (flicker_nr_t.("Power Sign")(amp_idxs) == 0);
            flicker_mouse_t{row_num, 'Not Power Sig'} = sum(sig_idx);

            row_num = row_num + 1;
        end
    end
end

flicker_mouse_t = sortrows(flicker_mouse_t, {'Stim', 'Mouse', 'Current'} )
writetable(flicker_mouse_t, [savefig_path f 'Flicker' f 'Neuronwise' f 'Mouse_Flicker_Only_Data.csv'], 'WriteRowNames', true);

