% This class is used to have the implementation of commonly called routines consolidated into a single file.
% The functions here are hopefully specific to the PV DBS project, and therefore may look similar to other previously made functions but with slight differences

classdef Multi_func
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
            F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3));%+ x(3)*exp(- xdata./x(4))  ;    
            x0 = [mean(v1) 40 1.5   ] ;    
            OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');    
            t = (1:length(v))./FS;    
            tsel=1:length(t);    
            [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);    
            fitbaseline=xunc(1)+xunc(2)*exp(-t./xunc(3));    
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

        % Re-space axis
        function [result] = set_spacing_axis(ax, spacing, num_dec)
            ticks = round(linspace(ax.Limits(1), ax.Limits(2), spacing), num_dec);    
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
    end
end
