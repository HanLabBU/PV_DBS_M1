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
        function [x, y]  = exp_fit(trace)
            t = 1:length(trace);
            f2 = fit(t', trace, 'exp2');
            y = f2.a*exp(f2.b*t) + f2.c*exp(f2.d*t);
            x = t;
        end

        % Calculate cwt for input signal and 
        function [wt, f] = get_power_spec(signal, samp_freq)
            freqLimits = [0 150];
            fb = cwtfilterbank(SignalLength=length(signal),...
                               SamplingFrequency=samp_freq,...
                               FrequencyLimits=freqLimits);
            [wt, f] = cwt(signal, FilterBank=fb);
        end
    end
end
