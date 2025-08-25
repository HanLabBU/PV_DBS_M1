 clear all
 
f = filesep;

rootpath = '/home/pierfier/Projects/';

cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\PV_Data\')


ses = dir(['*.mat']);
all_matfiles = {ses.name}
          
   [region_matfiles] = find_region(all_matfiles);        
            all_matfiles = {ses.name};
         
            region=0;
for f_region = fieldnames(region_matfiles)'
    region=region+1;
    f_region = f_region{1};

    % Create brain region struct
    region_stats.(f_region) = struct();
    
      %% Select matfiles by stim condition by this iteration's brain region
    [matfile_stim] = stim_cond(region_matfiles.(f_region).names);
  
      %% Loop through each field of the struct and concatenate everything together
    % Keep track the stats of the number of trials and neurons per mouse
    cond_stats = struct();
    % Loop through each stimulation condition
    cond=0
    for field = fieldnames(matfile_stim)'
        field = field{1};  
        cond=cond+1;
     
              matfiles = matfile_stim.(field).names;    
     
        % Initialize field to keep track of neurons and trials
        cond_stats.(field) = struct();
        cond_stats.(field).trial_nums = [];
        cond_stats.(field).num_fovs = 0;
    
    
        % Loop through each matfile of the current stimulation condition
        mr=0;
        for matfile =matfiles
              frs=1:2:150;FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([ matfile{1}]);
            clear alls wavD  allVm
            SNR=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
               if ~isempty(  data.align.trial{tr})
                    tr_iter= tr_iter+1;
                voltsig=data.raw.trial{tr}.raw_traces;
               resultS=spike_detect_SNR_sim3( voltsig,4,4,7);
              data.align.trial{tr}.spike_info2=  resultS;
               align= data.align; raw=data.raw;
              save([ matfile{1}], 'align', 'raw')
               end      
            end
            
              
        end
        
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


function  [filt_sig]=filt_data(sig,frs, FS)
Fn = FS/2;

for steps=1:length(frs);


FB=[ frs(steps)*0.9 frs(steps)*1.1];

 [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
 filt_sig(:,steps)= hilbert(filtfilt(B,A,sig));
            
 
end
end
