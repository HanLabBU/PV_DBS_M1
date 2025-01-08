% The functions here are hopefully specific to the PV DBS project, and therefore may look similar to other previously made functions but with slight differences
classdef Multi_func
    properties (Constant)
        
        % Path to save plot files
        % Into Dropbox
        save_plot = ['~/Dropbox/RKC-HanLab/Pierre PV DBS Project Dropbox/Materials/Plots/'];
        % Into the server
        %save_plot = ['~/handata_server/eng_research_handata3/Pierre Fabris/PV Project/Plots/'];

        % Specify the theta range use for filtering stuff
        theta_range = [2 10];
        theta_frequencies = [2:10];
        entr_freqs = [1:200];

        % Specify the transient and sustained time period
        trans_ped = [0, 100];
        sus_ped = [100, 1000];
        base_ped = [-500 0];
        stim_ped = [0 1000];
        offset_trans_ped = [1000, 1150];

        % Specify points for zoom-ins of pulses
        onset_ped = [-50, 100];
        mid_stim_ped = [500, 650];


        % Colors for violin plots
        trans_color = [153, 51, 51]/255;
        sus_color = [51, 51, 153]/255;
        %stim_color = [76, 149, 108]/255; nevermind lol

        % Colors indicating different stimulation periods
        base_color = [57, 77, 161]/255;
        stim_color = [131, 195, 65]/255;
        post_color = [128, 56, 149]/255;
        flick_color = [189, 130, 57]/255;

        % Colors indicating different brain regions
        ca1_color = [0.8500 0.3250 0.0980];
        m1_color = [0 0.4470 0.7410];
        v1_color = [106, 189, 69]/255;

        % Colors for activated, non-modulated, suppressed
        sup_color = [0 119 182]/255;
        non_color = [0 0 0]/255;
        act_color = [214 40 40]/255;

        % Colors for shuffled data
        shuf_color = [92, 161, 255]/255;

        % Colors indicating pulses
        pulse_color = [170, 176, 97]./255; %TODO should use dbs_color instead of pulse_color
        dbs_color = [202, 141,25]./255;        

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

        % Filter and perform hilbert transform for range
        function [filt_sig] = filt_range(sig, range, FS)
            % If range is a gradient, just set to minimum and maximum of gradient
            if length(range) > 2
                range = [min(range) max(range)];
            end
            Fn = FS/2;
            FB = [0.95 1.05].*range;

            [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
            filt_sig = hilbert(filtfilt(B,A,sig));
        end

        % Filter each trial by given frequency
        function [filt_sigs] = raw_filt(sigs, fr, FS)
            Fn = FS/2;
            FB = [fr*0.95 fr*1.05];
            for i = 1:size(sigs, 2)
                [B, A] = butter(2, [min(FB) max(FB)]./Fn);
                filt_sigs(:, i) = filtfilt(B, A, sigs(:, i));
            end
        end

        % Using FIR filter
        function [filt_sig] = fir_filt(sig, fr, Fs)
            N = 499;
            Wn = [fr*0.95 fr*1.05] / (Fs/2);
        
            b = fir1(N, Wn, 'bandpass', blackman(N + 1));
            filt_sig = filtfilt(b, 1, sig);
        end

        % Filter and grab the hilbert transform for at each frequency in the specified range
        function [filt_sig]=filt_data(sig,frs, FS)
            Fn = FS/2;
            for steps=frs;
                FB=[ frs(steps)*0.95 frs(steps)*1.05];
                
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                filt_sig(steps,:)= hilbert(filtfilt(B,A,sig));
            end
        end

        %Calculates spike phase-locking value
        function [PLV_output, PLV_output2, mat] = spike_field_ppcDBS_Pierre(wavD, spikes, timshift)
            if size(wavD,4)>1    
                [tr chs frs t]= size(wavD);    
            else    
                [t frs tr]= size(wavD);    
            end    
            mat=[];    
            for ind=1:tr    
                if size(wavD,4)>1    
                     M=exp(1i.*((squeeze(wavD(ind,CH,:,:)))));    
                else    
                    M=exp(1i.*((squeeze(wavD(:,:,ind))))); 
                end    
                s=find(spikes{ind});
                s=s-timshift; 
                s(s<=0)=[];    
                mat=[mat, M(s,:)']; % The normal PLV formula   
            end    
                
            NT= sum(~isnan(mat(1,:)));    
            Z=abs(nanmean(mat,2));  % MVL    
            T= Z.^2;    
            PLV_output= (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)    
            PLV_output2= mat;    
            
            % Exclusion criteria for how many events were detected
            if NT<=10    
                PLV_output=PLV_output.*NaN;    
                %PLV_output2=PLV_output2.*NaN;
            end
        end

        % Calculates the spike-phase locking value from spike to signal
        function [PLV, PLV2, phase_vecs] = spike_field_PLV(nr_phases, spikes, timeshift, exc_criteria)
            % Loop through each trial
            phase_vecs = [];
            for i=1:size(nr_phases, 3)
                trial_spikes = spikes(:, i);
                spike_idx = find(trial_spikes == 1);
                spike_idx - timeshift;
                spike_idx(spike_idx <= 0) = [];
                indi_phase_vecs = exp(1i.*(nr_phases(:, spike_idx, i)));
                
                %DEBUG
                if sum(~isnan(indi_phase_vecs), 'all') == 0
                    %sum(spikes)
                    %size(spikes)
                end

                phase_vecs = [phase_vecs, indi_phase_vecs];
            end

            % Transpose so columns indicate frequencies and rows are the spike phases
            phase_vecs = phase_vecs';
            
            num_spikes = size(phase_vecs, 1);
            PLV = (1/num_spikes)*(abs(sum(phase_vecs, 1)));
            PLV2 = (1/((num_spikes - 1))).*((PLV.^2).*num_spikes - 1);

            if num_spikes <= exc_criteria
                PLV2 = PLV.*NaN;
            end
        end

        % Calculate cwt for input signal and 
        function [wt, f] = get_power_spec(signal, samp_freq)
            freqLimits = [0 200];
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

        % Specify colormap of red and blue
        function [cmap] = get_red_blue_cmap()
            [red_blue_color_cmap] = (cbrewer('div', 'RdBu',500));
            red_blue_color_cmap(red_blue_color_cmap > 1) = 1;
            red_blue_color_cmap(red_blue_color_cmap < 0) = 0;
            red_blue_color_cmap = flipud(red_blue_color_cmap);
            cmap = red_blue_color_cmap;
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
            offset = 0;
            plot(x_pts, [y + offset, y + offset], '-', 'LineWidth', 2.5, 'Color', Multi_func.dbs_color);
            hold on;
            t1 = text(x_pts(1), y + 2*offset, text_str, 'FontSize', 7);
            txt_width = t1.Extent(3);
            x_offset = (diff(x_pts) - txt_width)/2;
            delete(t1);
            text(x_pts(1) + x_offset, y + 2*offset, text_str, 'FontSize', 7);
        end

        % Determine if the population neuron data has a significant low frequency data
        % Will return an array that will indicate this
        function [has_low_freq] = get_low_freq_neurons(popul_data)
            has_low_freq = [];
            
            base_all_power_norm = [];
            power_norm_neuron_name = {};

            % Range of frequencies to zscore by
            low_freqs = [1:50];

            % Loop through each neuron
            for nr=1:size(popul_data.neuron_Vm, 2)
                base_idx = find(popul_data.trace_timestamps(:, nr) < popul_data.stim_timestamps(1, nr));
                nr_filt_data = popul_data.neuron_hilbfilt{nr};

                base_power = [];
                % Loop through each trial
                for tr=1:size(nr_filt_data, 3)
                    base_power(:, end + 1) = nanmean(abs(nr_filt_data(low_freqs, base_idx, tr)), 2);
                end
                
                % Calculate the average power
                avg_power = nanmean(base_power, 2);
                std_power = std(base_power, 0, 2, 'omitnan');
                num_trials = size(base_power, 2);
                sem_power = std_power./sqrt(num_trials);
                
                timeline = popul_data.trace_timestamps(:, nr);
                freqs = 1:size(base_power, 1);

                % z-score the power across frequency
                base_all_power_norm(:, end + 1) = zscore(avg_power(:, end), [], 1);
                power_norm_neuron_name{end + 1} = popul_data.neuron_name{nr};

                % Determine if neuron has low frequency oscillations
                has_low_freq(end + 1) = nanmean(base_all_power_norm(Multi_func.theta_frequencies, end), 1) > 1; % Currently set to higher than 1 std of all frequencies
            end
        end

        % Non-linearly normalize each trial in the matrix, this function also assumes that the trials
        % are columnwise
        function [result] = norm_signals(x)
            result = (x - min(x, [], 1))./(max(x, [], 1) - min(x, [], 1));
        end

        % Create fields in structure if they do not exist
        function [result] = create_struct(src, f_region, f_mouse, f_rec, f_neuron, f_stim)
            result = src;
            % Create region field
            if isfield(result, f_region) == 0
                result.(f_region) = struct();
            end
            
            % Create mouse field
            if isfield(result.(f_region), f_mouse) == 0
                result.(f_region).(f_mouse) = struct();
            end
        
            % Create recording field
            if isfield(result.(f_region).(f_mouse), f_rec) == 0
                result.(f_region).(f_mouse).(f_rec) = struct();
            end
        
            % Create stim field
            if isfield(result.(f_region).(f_mouse).(f_rec), f_stim) == 0
                result.(f_region).(f_mouse).(f_rec).(f_stim) = struct();
                result.(f_region).(f_mouse).(f_rec).(f_stim).currents = [];
            end
        
            % Create neuron field
            % Old code for saving everything into the neuron field's
            %if isfield(result.(f_region).(f_mouse).(f_rec).(f_stim), f_neuron) == 0
            %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron) = struct();
            %    result.(f_region).(f_mouse).(f_rec).(f_stim).(f_neuron).currents = [];
            %    
            %end
        end

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
            % Calculate the front edge with smaller window sizes
            %edge_rate = Multi_func.recurse_filter(impulse_raster(1:wind - 1));
            %result(1:wind - 1) = edge_rate;
        end
    end
end
