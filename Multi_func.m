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
    end
end
