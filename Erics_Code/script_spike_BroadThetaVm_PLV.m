clear all
 
f = filesep;

% Linux computer paths
server_root_path = '~/handata_server/';
local_root_path = '~/Projects/';
pv_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
addpath([local_root_path f 'Pierre Fabris' f 'PV DBS neocortex' f 'Scripts' f]);

savepath = Multi_func.save_plot;

for id=1:2; for id2=1:3; allPLVb2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVs2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVp2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allVmDBS{id,id2}=[];end;end

ses = dir([pv_path '*.mat']);
all_matfiles = {ses.name}
          
[region_matfiles] = Multi_func.find_region(all_matfiles);        
all_matfiles = {ses.name};

%% Loop through each region matfile
region=0;
for f_region = fieldnames(region_matfiles)'
    region=region+1;
    f_region = f_region{1};

    % Create brain region struct
    region_stats.(f_region) = struct();
    
    % Select matfiles by stim condition by this iteration's brain region
    [matfile_stim] = Multi_func.stim_cond(region_matfiles.(f_region).names);
  
    % Loop through each field of the struct and concatenate everything together
    % Keep track the stats of the number of trials and neurons per mouse
    cond_stats = struct();
    % Loop through each stimulation condition
    %DBS_cond= {'f_40'; 'f_140' ; 'f_1000'};
    DBS_cond= {'f_40'; 'f_140'};
    DBS_reg = fieldnames(region_matfiles);
    cond=0
    for dbs =1:2
        field = DBS_cond{dbs};  
        cond=cond+1;
     
              matfiles = matfile_stim.(field).names;    
     
        % Initialize field to keep track of neurons and trials
        cond_stats.(field) = struct();
        cond_stats.(field).trial_nums = [];
        cond_stats.(field).num_fovs = 0;
    
    
        % Loop through each matfile of the current stimulation condition
        mr=0;
        for matfile =matfiles

              %frs=[1:1:10 11:1:150 152:2:180];
              frs = [2:10];
              FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([ matfile{1}]);
            clear alls wavD  allVm
            SNR=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
               if ~isempty(  data.align.trial{tr}) 
                   %actually filter by 4 SNR here
                   if mean(data.align.trial{tr}.spike_info375.spike_snr{1})>0 % Trial included if spike SNR above 4SNR
                    tr_iter= tr_iter+1;
                voltsig=data.align.trial{tr}.spike_info375.trace_ws;
                [ fitbaseline, coeff]=Multi_func.exp_fit_Fx(voltsig',FS); %remove photobleaching
                voltsig=(voltsig-fitbaseline);
                
                if  0
                    
                   Camtime=  data.raw.trial{tr}.raw_camera_frame_time- data.raw.trial{tr}.raw_camera_start_time(1);
DBStim=data.raw.trial{tr}.raw_stimulation_time- data.raw.trial{tr}.raw_camera_start_time(1);

traceDBS=zeros(1,length(data.raw.trial{tr}.raw_traces));
for id=1:length(DBStim)
z=find(Camtime > DBStim(id));try
traceDBS(z(1))=1;end
end
 
                    
                end
                allVm(:, tr_iter)=voltsig;
                
                [filt_sig]=Multi_func.filt_range(voltsig,frs, FS); % filtering+Hilbert
                wavD(:,:, tr_iter)=filt_sig';
              SNR  =[SNR;(data.align.trial{tr}.spike_info375.spike_snr{1})];
                alls(:, tr_iter)=data.align.trial{tr}.spike_info375.roaster;
                   end 
               end
            end
       if    tr_iter   >0 
            baseTimsel= [10:FS-10];
            StimTimsel= [FS:FS*2];
            PostTimsel= [FS*2+10:size(alls,1)];
            
   
if  mean(SNR) >0 % mean spike SNR for the neuron to be over 4SNR (can be changed)
 
    
            shift_time=0; % shift from spike backwards
clear s M
M=angle(wavD(baseTimsel,:,:));
  for tr1=1:size(alls,2); 
      s{tr1}=alls(baseTimsel,tr1)';
  end
[ PLV_b, PLV_b2]= Multi_func.spike_field_ppcDBS_Pierre(M,s,  shift_time);  % PLV estimatiom

    clear s M;M=angle(wavD( StimTimsel,:,:));
  for tr1=1:size(alls,2); 
      s{tr1}=alls( StimTimsel,tr1)';
  end
[ PLV_s, PLV_s2]= Multi_func.spike_field_ppcDBS_Pierre(M,s,  shift_time); 

    clear s M;M=angle(wavD(  PostTimsel,:,:));
  for tr1=1:size(alls,2); 
      s{tr1}=alls(  PostTimsel,tr1)';
  end
[ PLV_p, PLV_p2]= Multi_func.spike_field_ppcDBS_Pierre(M,s,  shift_time); 

mr=mr+1


  allVmDBS{region,cond}(1:length(allVm),mr) = nanmean(allVm,2);

 allPLVb{region,cond}(1:length(PLV_b),mr)=PLV_b;
  allPLVs{region,cond}(1:length(PLV_s),mr)=PLV_s;
   allPLVp{region,cond}(1:length(PLV_p),mr)=PLV_p;
   
    allPLVb2{region,cond}=[allPLVb2{region,cond},  PLV_b2];
      allPLVs2{region,cond}=[allPLVs2{region,cond},  PLV_s2];
      allPLVp2{region,cond}=[allPLVp2{region,cond},  PLV_p2];
end  
end
end
end
end

pheight = 350;

%% Violin plots of broad band low frequency spike-Vm PLV
for area = 1:2
    for cond = 1:2
        figure('COlor','w','Position', [ 300 400 380 pheight],'Renderer', 'painters')
        num_vals = length(allPLVb{area,cond});
        data = [allPLVb{area,cond}, allPLVs{area,cond}]; 
        labels = [repmat({'Base'}, num_vals, 1); repmat({'Stim'}, num_vals, 1)];

        % Plot the violins
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);

        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};

        % Statistical Data
        [p, h, stats] = signrank(allPLVb{area,cond}, allPLVs{area,cond});
        fprintf('\n\n');
        disp('Signrank for LowFreq spike-Vm PLV');
        disp([DBS_reg{area} ' ' DBS_cond{cond}]);
        disp(['p=' num2str(p)]);
        disp(['stats=']);
        disp(struct2table(stats));
        fprintf('\n\n');

        title([DBS_reg{area} ' ' DBS_cond{cond}], 'Interpreter', 'none');
        saveas(gcf, [savepath 'PLV' f 'Violin_PLVBroadTheta_' DBS_reg{area} '_' DBS_cond{cond} '.pdf']);
        saveas(gcf, [savepath 'PLV' f 'Violin_PLVBroadTheta_' DBS_reg{area} '_' DBS_cond{cond} '.png']);
    end
end

% NOTE: This does not make sense because we need the individual frequency points
%% Plotting Spike-Theta Vm PLV using broad band filtering
for area=1:2
    for cond=1:2
            M=  [allPLVb{area,cond}] ; 
            M1= [allPLVs{area,cond} ] ; 
            %  M2= [allPLVp{1,cond},allPLVp{2,cond}] ;
            figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters')
            fill_error_area2(frs, nanmean(M,2),nanstd(M,[],2)./sqrt(sum(~isnan((M(1,:))))  ) ,[ .5 .5 .5] )
            fill_error_area2(frs, nanmean(M1,2),nanstd(M1,[],2)./sqrt(sum(~isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
            %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] )
            plot(frs,nanmean(M,2),'b','Linewidth',1); hold on; plot(frs,nanmean(M1,2),'g','Linewidth',1)
            %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
            axis tight  ; ylim([ -0.02 0.07]); 
            xlim([ 2 10])
            % set(gca,'Xscale','log')
            title([DBS_reg{area}(3:end) ' ' DBS_cond{cond}(3:end)], 'Interpreter', 'none');
            print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'PLV_BroadTHETA' DBS_reg{area}(3:end) '_area_high_' DBS_cond{cond}(3:end) '.pdf'])
    end
end
