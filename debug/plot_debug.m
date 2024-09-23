load('detrend_debug.mat');

for f_region = fieldnames(region_data)'
    f_region = f_region{1};

    for f_stim = fieldnames(region_data.(f_region))'
        f_stim = f_stim{1};
        

        figure('Position', [1, 1, 800, 600]);
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_lin_vm, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_lin_vm, 0, 2, 'omitnan');
        num_vm = size(region_data.(f_region).(f_stim).all_base_lin_vm, 2);
        base_sem = base_std./sqrt(num_vm);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_lin_vm, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Linear fit baseline only');
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_offset_exp_vm, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_offset_exp_vm, 0, 2, 'omitnan');
        num_vm = size(region_data.(f_region).(f_stim).all_base_offset_exp_vm, 2);
        base_sem = base_std./sqrt(num_vm);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_offset_exp_vm, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Exponential baseline and offset period');
        
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_filt_vm, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_filt_vm, 0, 2, 'omitnan');
        num_vm = size(region_data.(f_region).(f_stim).all_filt_vm, 2);
        base_sem = base_std./sqrt(num_vm);
        timeline = [1:size(region_data.(f_region).(f_stim).all_filt_vm, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window across all points');

        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_off_filt_vm, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_off_filt_vm, 0, 2, 'omitnan');
        num_vm = size(region_data.(f_region).(f_stim).all_base_off_filt_vm, 2);
        base_sem = base_std./sqrt(num_vm);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_off_filt_vm, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window, base and offset points (average across stim period)');

        sgtitle([f_region '_' f_stim], 'Interpreter', 'none');

        % Plot the fits together TODO add all of this
        figure('Position', [1, 1, 800, 600]);
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_lin_fit, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_lin_fit, 0, 2, 'omitnan');
        num_fit = size(region_data.(f_region).(f_stim).all_base_lin_fit, 2);
        base_sem = base_std./sqrt(num_fit);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_lin_fit, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Linear baseline only fits');
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_offset_exp_fit, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_offset_exp_fit, 0, 2, 'omitnan');
        num_fit = size(region_data.(f_region).(f_stim).all_base_offset_exp_fit, 2);
        base_sem = base_std./sqrt(num_fit);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_offset_exp_fit, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Exponential baseline and offset period fits');
        
        
        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_filter_fit, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_filter_fit, 0, 2, 'omitnan');
        num_fit = size(region_data.(f_region).(f_stim).all_filter_fit, 2);
        base_sem = base_std./sqrt(num_fit);
        timeline = [1:size(region_data.(f_region).(f_stim).all_filter_fit, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window across all points fits');

        nexttile;
        base_avg = mean(region_data.(f_region).(f_stim).all_base_offset_filter_fit, 2, 'omitnan');
        base_std = std(region_data.(f_region).(f_stim).all_base_offset_filter_fit, 0, 2, 'omitnan');
        num_fit = size(region_data.(f_region).(f_stim).all_base_offset_filter_fit, 2);
        base_sem = base_std./sqrt(num_fit);
        timeline = [1:size(region_data.(f_region).(f_stim).all_base_offset_filter_fit, 1)]';
        fill_h = fill([timeline; flip(timeline)], [base_avg + base_sem; flipud(base_avg - base_sem)], [0.5 0.5 0.5]);
        Multi_func.set_fill_properties(fill_h);
        hold on;
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window, base and offset points (average across stim period) fits');

        sgtitle([f_region '_' f_stim], 'Interpreter', 'none');

        % Plot all of the fits
        timeline = repmat(timeline(:), 1, num_fit);

        figure('Position', [400, 1, 800, 600]);
        tiledlayout(4, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        nexttile;
        base_avg = region_data.(f_region).(f_stim).all_base_lin_fit;
        %num_fit = size(region_data.(f_region).(f_stim).all_base_only_exp_fit, 2);
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Linear baseline only fits');
        
        nexttile;
        base_avg = region_data.(f_region).(f_stim).all_base_offset_exp_fit;
        %num_fit = size(region_data.(f_region).(f_stim).all_base_offset_exp_fit, 2);
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Exponential baseline and offset period fits');
        
        nexttile;
        base_avg = region_data.(f_region).(f_stim).all_filter_fit;
        %num_fit = size(region_data.(f_region).(f_stim).all_filter_fit, 2);
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window across all points fits');

        nexttile;
        base_avg = region_data.(f_region).(f_stim).all_base_offset_filter_fit;
        %num_fit = size(region_data.(f_region).(f_stim).all_base_offset_filter_fit, 2);
        plot(timeline, base_avg, 'k', 'LineWidth', 0.3);
        title('Moving window, base and offset points (average across stim period) fits');

        sgtitle([f_region '_' f_stim '_all_fits'], 'Interpreter', 'none');
    end
end
