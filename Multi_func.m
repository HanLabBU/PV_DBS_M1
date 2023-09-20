%% This class is used to have the implementation of commonly called routines consolidated into a single file.
% The functions here are hopefully specific to the PV DBS project, and therefore may look similar to other previously made functions but with slight differences

classdef Multi_func
    properties (Constant)
        
        % Colors for violin plots
        trans_color = [153, 51, 51]/255;
        sus_color = [51, 51, 153]/255;
        %stim_color = [76, 149, 108]/255; nevermind lol

        base_color = [57, 77, 161]/255;
        stim_color = [131, 195, 65]/255;
        post_color = [128, 56, 149]/255;

        % Colors for heatmaps
        light_gray_color = [linspace(256, 178, 256)', linspace(256, 178, 256)', linspace(256, 178, 256)']./256;
        tang_blue_color = [[linspace(38, 233, 256/2)', linspace(70, 196, 256/2)', linspace(83, 106, 256/2)'];
                          [linspace(233, 231, 256/2)', linspace(196, 111, 256/2)', linspace(106, 81, 256/2)']]./256;                    

        warm_cold_color = [[linspace(58, 256, 256/2)', linspace(12, 256, 256/2)', linspace(163, 256, 256/2)'];
                         [linspace(256, 217, 256/2)', linspace(256, 4, 256/2)', linspace(256, 41, 256/2)']]./256;
        
        warm_cold_gray_color = [[linspace(58, 200, 256/2)', linspace(12, 200, 256/2)', linspace(163, 200, 256/2)'];
                         [linspace(200, 217, 256/2)', linspace(200, 4, 256/2)', linspace(200, 41, 256/2)']]./256;
               
        green_warm_cold_color = [[linspace(202, 160, 256/2)', linspace(255, 196, 256/2)', linspace(191, 255, 256/2)'];
                         [linspace(160, 217, 256/2)', linspace(196, 4, 256/2)', linspace(255, 41, 256/2)']]./256;
    end

    methods(Static)

        % Read in .csv file to a dictionary structure that stores which traces to ignore    
        function [result_dict] = csv_to_struct(csv_pathname)    
            global func_debug;    
    
            % Grab the Variables names at the top of the csv    
            fid = fopen(csv_pathname);    
            varNames = strsplit(fgetl(fid), ';');    
            fclose(fid);    
    
            % Read in csv file        
            opts = detectImportOptions(csv_pathname);    
            opts = setvartype(opts, {'mouse', 'rec_date', 'FOV', 'cond', 'roi_trial_dict'}, 'string');    
    
            T = readtable(csv_pathname, opts);    
    
            % Loop through each row and continusly add fields to the dictionary, while ignoring comments       
            for i=1:length(T.Var1)    
                mouse_key = ['mouse_' char(T.mouse(i))];    
                rec_key = ['rec_' char(T.rec_date(i))];    
                FOV_key = ['FOV' char(T.FOV(i))];
                cond_key = ['f_' char(T.cond(i))];
                result_dict.(mouse_key).(rec_key).(FOV_key).(cond_key) = Ignore_trials_class.str_to_dict(T.roi_trial_dict(i));    
            end    
        end

        % Select specific parameters from matfiles

        %% Specific functions for determining which FOVs to look at
        function [region_struct] = find_region(matfile_names)
            region_struct = struct();
        
            for i=1:length(matfile_names)
                file_parts = split(matfile_names{i}, '_');
                reg = file_parts{2};
                
                % Create region field if it does not exist
                if ~isfield(region_struct, ['r_' reg])
                    region_struct.(['r_' reg]).names = {};
                end
        
                region_struct.(['r_' reg]).names{end+1} = matfile_names{i};
            end
        end
        
        % Return matfiles by stimulation condition
        function [cond_struct] = stim_cond(matfile_names)
            cond_struct = struct();
            
            % Loop through each matfilename and group by stimulation conditions
            for i=1:length(matfile_names)
                    file_parts = split(matfile_names{i}, '_');
                    stim = file_parts{5};
                    
                    % Create stimulation field if it does not exist
                    if ~isfield(cond_struct, ['f_' stim])
                        cond_struct.(['f_' stim]).names = {};
                    end
        
                    cond_struct.(['f_' stim]).names{end+1} = matfile_names{i};
            end
        end

        % Return line of best exponential fit
        function [baseline, coeff]  = exp_fit(trace)
            t = 1:length(trace);
            f2 = fit(t', trace, 'exp2');
            y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
            baseline = y;
            coeff = f2;
        end

        % Return sophisticated exponential fit accounting for stimulation depolarization
        function [ fitbaseline, coeff]=exp_fit_Fx(v,FS)    
            %% fits an expoential to estimate the photobleaching    
            v1=v;    
            v1([FS-10:FS*2+20])= mean([ v([FS-40:FS-20 2*FS+20:2*FS+40 ])]);    
            v1(1)=v1(2);    
            F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3)); %+ x(3)*exp(- xdata./x(4))  ;    
            x0 = [mean(v1) 40 1.5] ;    
            OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');    
            t = (1:length(v))./FS;
            % Whole trace
            %tsel=[1:length(t)];    
            %Base and offset period consideration
            tsel=[1:FS 2*FS:length(t)];    
            [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);    
            fitbaseline=F(xunc, t);    
            %DEBUG
            %figure; tiledlayout(2, 1); nexttile; plot(v); hold on; plot(fitbaseline);
            %nexttile;
            %plot(v' - fitbaseline);
            coeff=xunc;    
        end

        % TODO change so only the baseline points are considerred
        % Exponential fit that only accounts for the exponential in the baseline
        function [ fitbaseline, coeff]=exp_fit_Fx_Base(v,FS)    
            %% fits an expoential to estimate the photobleaching    
            v1=v;    
            v1([FS-10:end])= mean([ v([FS-40:end ])]);    
            v1(1)=v1(2);    
            F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3));%+ x(3)*exp(- xdata./x(4))  ;    
            x0 = [mean(v1) 40 1.5   ] ;    
            OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');    
            t = (1:length(v))./FS;
            tsel=1:length(t);    
            [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);    
            fitbaseline=xunc(1)+xunc(2)*exp(-t./xunc(3));
            %DEBUG
            %figure; tiledlayout(2, 1); nexttile; plot(v); hold on; plot(fitbaseline);
            %nexttile;
            %plot(v' - fitbaseline);
           
            coeff=xunc;    
        end

        % Calculate cwt for input signal and 
        function [wt, f] = get_power_spec(signal, samp_freq)
            freqLimits = [0 150];
            fb = cwtfilterbank(SignalLength=length(signal),...
                               SamplingFrequency=samp_freq,...
                               FrequencyLimits=freqLimits);
            [wt, f] = cwt(signal, FilterBank=fb);
        end

        % Sets default values of the violin plots  
        function [opts] = get_default_violin()
            opts.ShowMedian = true;
            opts.ShowMean = false;
            opts.MedianColor = [1 1 1];
            opts.MarkerSize = 5;
            opts.MedianMarkerSize = 5;
            opts.BoxWidth = 0.1;
            opts.BoxColor = [0 0 0];
            opts.ViolinAlpha = {[0.3], [0.3]};
        end

        % Specify the fill property for all figures
        % fill_handle -> the fill to set all of these properties to make them uniform
        function [result] = set_fill_properties(fill_handle)
            fill_handle.EdgeAlpha = 1;
            fill_handle.FaceAlpha = 0.2;
            fill_handle.LineWidth = 0.2;
        end
        
        % Specify default axis on plots
        function [result] = set_default_axis(ax)
            set(ax, 'Color', 'none', 'Box', 'off', 'TickDir', 'out', 'linewidth', 0.2);
            ax.FontSize = 7;
        end

        % Space tick marks starting from 0 and going towards the limits
        function [result] = set_spacing_axis(ax, spacing, num_dec)
            upper_ticks = 0:spacing:ax.Limits(2);
            lower_ticks = 0:-spacing:ax.Limits(1);
            ticks = round(union(lower_ticks, upper_ticks), num_dec);    
            ax.TickValues = ticks;
            ax.TickLabels = num2cell(ticks);
        end

        % Plot DBS bar above specified value
        function [result] = plot_dbs_bar(x_pts, y, text_str)
            offset = 1;
            plot(x_pts, [y + offset, y + offset], '-', 'LineWidth', 2.5, 'Color', Fig_color_props.dbs_color);
            hold on;
            t1 = text(x_pts(1), y + 2*offset, text_str, 'FontSize', 7);
            txt_width = t1.Extent(3);
            x_offset = (diff(x_pts) - txt_width)/2;
            delete(t1);
            text(x_pts(1) + x_offset, y + 2*offset, text_str, 'FontSize', 7);
        end

        %TODO change the data_bystim part of the struct
        % Combine all regions into a single 'region' structure called 'f_combined'
        function [combine_struct] = combine_regions(region_data)
            % Initialize combined region
            combine_struct.r_combine = struct();
        
            % Grab all of the fields from the strutures
            f_regions = fieldnames(region_data)';
            f_stims = fieldnames(region_data.(f_regions{1}))';
            f_data = fieldnames(region_data.(f_regions{1}).(f_stims{1}))';
            
            % Initialize data strutures in the combined field
            data_bystim = struct();
            for f_stim = f_stims
                f_stim = f_stim{1};
                data_bystim.(f_stim) = cell2struct(cell(size(f_data)), f_data, 2);
            end
        
            % Loop through each region
            for f_region = f_regions
                f_region = f_region{1};
                % Loop through each stimulation condition
                for f_stim = f_stims
                    f_stim = f_stim{1};
                    % Loop through each data field
                    for f_datum = f_data
                        f_datum = f_datum{1};
                        dim = length( size( region_data.(f_region).(f_stim).(f_datum) ) );
                        if dim == 2
                            data_bystim.(f_stim).(f_datum) = horzcat_pad(data_bystim.(f_stim).(f_datum), region_data.(f_region).(f_stim).(f_datum));
                        else
                            data_bystim.(f_stim).(f_datum) = cat(dim, data_bystim.(f_stim).(f_datum), region_data.(f_region).(f_stim).(f_datum));       
                        end
                    end
                end
            end
            combine_struct.r_combine = data_bystim;
        end
        
        %-- Depracated ---
        % Combine all regions into a single 'region' structure called 'f_combined'
        function [combine_struct] = combine_regions_old(region_data)
            % Initialize combined region
            combine_struct.r_combine = struct();
            combine_struct.r_combine.data_bystim = struct();
        
            % Grab all of the fields from the strutures
            f_regions = fieldnames(region_data)';
            f_stims = fieldnames(region_data.(f_regions{1}).data_bystim)';
            f_data = fieldnames(region_data.(f_regions{1}).data_bystim.(f_stims{1}))';
            
            % Initialize data strutures in the combined field
            data_bystim = struct();
            for f_stim = f_stims
                f_stim = f_stim{1};
                data_bystim.(f_stim) = cell2struct(cell(size(f_data)), f_data, 2);
            end
        
            % Loop through each region
            for f_region = f_regions
                f_region = f_region{1};
                % Loop through each stimulation condition
                for f_stim = f_stims
                    f_stim = f_stim{1};
                    % Loop through each data field
                    for f_datum = f_data
                        f_datum = f_datum{1};
                        dim = length( size( region_data.(f_region).data_bystim.(f_stim).(f_datum) ) );
                        if dim == 2
                            data_bystim.(f_stim).(f_datum) = horzcat_pad(data_bystim.(f_stim).(f_datum), region_data.(f_region).data_bystim.(f_stim).(f_datum));
                        else
                            data_bystim.(f_stim).(f_datum) = cat(dim, data_bystim.(f_stim).(f_datum), region_data.(f_region).data_bystim.(f_stim).(f_datum));       
                        end
                    end
                end
            end
            combine_struct.r_combine.data_bystim = data_bystim;
        end

        % Perform the Wilcoxon test and save results into a file
        % The test label will indicate what this test is for
        %function [result] wilc_test(x1, x2, test_label, stats_table)
        %    [p, h, stats] = ranksum(x1, x2);
        %    
        %end

        % This will recursively perform regular sliding averages with decreasing point size
        function [result] = recurse_filter(arr)
            if isempty(arr)
                result = [];
            elseif length(arr) == 1
                result = arr;
            else
                result = [Multi_func.recurse_filter(arr(1:end - 1)), mean(arr, 'omitnan')];
            end
        end

        % Perform spike rate estimation
        function [result] = estimate_spikerate(spike_raster, Fs, wind)
            impulse_raster = spike_raster.*Fs;
            trans_coeff = ones(wind, 1)./wind;
            impulse_raster = [impulse_raster(wind:-1:2), impulse_raster];

            result = filter(trans_coeff, 1, impulse_raster);
            result(1:wind - 1) = [];
            % This method unforunately created weird edges
            %% Calculate the front edge with smaller window sizes
            %edge_rate = Multi_func.recurse_filter(impulse_raster(1:wind - 1));
            %result(1:wind - 1) = edge_rate;
        end
    end
end
