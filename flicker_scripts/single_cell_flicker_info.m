%% This file requires flicker_plv.m and flicker_spectra.m to show the single-cell stats

%% Plot each neuron's PLV distribution, power, and coherence with 8 Hz
flick_base_tmps = [0 1];
stim_timestamps = [1 2];
gray = [0.5 0.5 0.5];
freqs = Multi_func.entr_freqs;
set_xlim = [0, 10];

% Interpolate the frequency space        
freq_lin = linspace(Multi_func.entr_freqs(1), Multi_func.entr_freqs(end), 200);
for f_stim = fieldnames(data)' %{'f_140'} %
    f_stim = f_stim{1};
    popul_data = data.(f_stim);
    timeline = popul_data.interp_time;

    % Loop through all of the neurons
    for f_nr = fieldnames(popul_data.nr_name)' %{'n_12'}%
        f_nr = f_nr{1};
        
        fig = figure('Position', [400 400 1800 1200]);
        
        tiledlayout(3, 1);

        % Plot the power
        nexttile;
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

        % Check if neuron's power spec is significant
        sig = [];
        if ismember(f_nr, popul_data.coh_mod.vpo_24{1})
            sig = 'Inc';
        elseif ismember(f_nr, popul_data.coh_mod.vpo_24{2})
            sig = 'Dec';
        else
            sig = '';
        end
        title(['Coherence ' sig ' ' num2str(popul_data.coh_vals.vpo_24.stats.(f_nr))]);

        % Plot the PLV
        nexttile;
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

        sgtitle([f_stim 'Hz ' f_nr ' ' popul_data.nr_name.(f_nr)], 'Interpreter', 'none');

        %fontsize(fig, 20, "points");
        set(findall(fig,'-property','FontSize'), 'FontSize', 20);

        exportgraphics(fig, [savefig_path 'Flicker' f 'Neuronwise' f '8Hz_plots' f ...
           f_stim(3:end) '_' popul_data.nr_name.(f_nr) '_.png'], 'Resolution', 600);
   
    end
end

%% Create a table displaying all of the info for single-cells
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

        % Grab the recording date
        recording_date = name_parts{4};

        % Parse out the current amplitude
        amp_str = popul_data.stim_amp.(f_nr);
        start_idx = find(amp_str == ':') + 1;
        end_idx = find(amp_str == 'u', 1, 'last') - 1;
        
        amp_str = amp_str(start_idx:end_idx);

        % TODO re-export the data from Python to MATLAB
        % Make sure that the save matfile has essentially the same data and
        % its just the current amplitude that is being added

        % Add parameters to table
        flicker_nr_t{row_num, 'Mouse_ID'} = mouse_name;
        flicker_nr_t{row_num, 'Stim Frequency'} = str2num(f_stim(3:end));
        flicker_nr_t{row_num, 'Current'} = amp_str;
        
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

        row_num = row_num + 1;
    end
end

%DEBUG
flicker_nr_t
writetable(flicker_nr_t, [savefig_path f 'Neuronwise' f 'Neurons_Flicker_Data.csv'], 'WriteRowNames', true);

