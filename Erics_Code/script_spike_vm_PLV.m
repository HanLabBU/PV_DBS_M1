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

              frs=[1:1:10 11:1:150 152:2:180];
              FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([ matfile{1}]);
            clear alls wavD  allVm
            SNR=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
               if ~isempty(  data.align.trial{tr}) 
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
                
                [filt_sig]=Multi_func.filt_data(voltsig,frs, FS); % filtering+Hilbert
                wavD(:,:, tr_iter)=filt_sig;
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
    
% Original savepath
%savepath='\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\Figures_Eric\PLV_VM\'
pheight=150;

%% Show PLV for each regiona and frequency with baseline in blue and stimulation in red
mm=0;  
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
for region=1:2;
    for cond=1:2
        mm=mm+1;
        subplot(2,3,mm)
        plot(frs,nanmean(allPLVb{region,cond},2),'b','Linewidth',1); hold on
        plot(frs,nanmean(allPLVs{region,cond},2),'r','Linewidth',1)
        %plot(frs,nanmean(allPLVp{region,cond},2),'m','Linewidth',1)
        axis tight;set(gca,'Xscale','log')
        ylim([ -0.05 0.33]);
        title([DBS_reg{region}(3:end) ' ' DBS_cond{cond}(3:end)]);
    end;
end;
    
%  for cond=1:2
%   M= [allPLVb{1,cond},allPLVb{1,cond}] ; M1= [allPLVs{1,cond},allPLVs{1,cond}] ;
%     M2= [allPLVp{1,cond},allPLVp{2,cond}] ;
%     figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
% fill_error_area2(frs, nanmean(M,2),nanstd(M,[],2)./sqrt(sum(isnan((M(1,:))))  ) ,[ .5 .5 .5] )
% fill_error_area2(frs, nanmean(M1,2),nanstd(M1,[],2)./sqrt(sum(isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
% %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] ) 
% plot(frs,nanmean(M,2),'b','Linewidth',1); hold on; plot(frs,nanmean(M1,2),'g','Linewidth',1)
%  %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
%  axis tight  ; ylim([ -0.05 0.2]); xlim([ 2 20])
%    % set(gca,'Xscale','log')
%    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_across_area_low_' num2str(cond) '.pdf'])
%  end
%

%% PLV for each region and stimulation frequency
for area=1:2
    for cond=1:2
        M= [allPLVb{area,cond}] ;
        M1= [allPLVs{area,cond}] ;
        %  M2= [allPLVp{1,cond},allPLVp{2,cond}] ;
        figure('COlor','w','Position', [ 300 400 400 pheight],'Renderer', 'painters')
        fill_error_area2(frs, nanmean(M,2),nanstd(M,[],2)./sqrt(sum(~isnan((M(1,:))))  ) ,[ .5 .5 .5] )
        fill_error_area2(frs, nanmean(M1,2),nanstd(M1,[],2)./sqrt(sum(~isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
        %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] )
        plot(frs,nanmean(M,2),'b','Linewidth',1); hold on; plot(frs,nanmean(M1,2),'g','Linewidth',1)
        %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
        axis tight; 
        ylim([ -0.05 0.3]); 
        xlim([ 2 160]);
        ax = gca;
        ax.Units = "centimeters";
        ax.InnerPosition = [1 1 3.5118 3.2964];
        % set(gca,'Xscale','log')
        title([DBS_reg{area}(3:end) ' ' DBS_cond{cond}(3:end) ' Spike-Vm PLV']);
        print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'PLV_' DBS_reg{area}(3:end) '_' DBS_cond{cond}(3:end) '_spike_vm.pdf']);
        saveas(gcf, [ savepath 'PLV' f 'PLV_' DBS_reg{area}(3:end) '_' DBS_cond{cond}(3:end) '_spike_vm.png']);
    end
end
    

%% Violin (stim freq) of spike-vm PLV at 140Hz stimulation for 140Hz neurons
%for area=1:2;
%    cond=2;
%    M= [allPLVb{area,cond}] ; M1= [allPLVs{area,cond}] ;M2= [allPLVp{area,cond}] ;
%    fsel=frs>138 & frs< 142;  
%     MM1= nanmean( M(fsel,:),1)'; MM= nanmean( M1(fsel,:),1)'; MM2= nanmean( M2(fsel,:),1)';
%    figure('COlor','w','Position', [ 300 400 120 140],'Renderer', 'painters')
%    violinplot2(MM1,[1.3 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
%    hold on,
%    violinplot2(MM,[ 1.8],'ViolinColor', [ 0.5 0.8 0; 0.9 0 0])
%    violinplot2(MM2,[ 2.3],'ViolinColor', [ 0.7 0.4 0.2; 0.9 0 0])
%    line([ 0.7 2.7], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
%    xlim([ 0.9 2.7]); %ylim([ -0.1 0.7])
%       print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'Violin_stimFrq_spike_vm_plV_' DBS_reg{area}(3:end) '_cond_' DBS_cond{cond}(3:end) '.pdf'])
%     [h,p,ci,stats] = ttest(MM,MM1)
%end
  

%% Violin plotSpike-Vm PLV at stimulation frequency
stats_log = [savepath 'PLV' f 'StimFreq_spike_vm_violin_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

for area=1:2;
    for cond=1:2
        M= [allPLVb{area,cond}] ; M1= [allPLVs{area,cond}] ;M2= [allPLVp{area,cond}] ;
        if cond == 1
            fsel=frs > 38 & frs < 42;  
        elseif cond == 2
            fsel=frs > 138 & frs < 142;  
        end
        
        MM= mean( M(fsel,:),1, 'omitnan')'; 
        MM1= mean( M1(fsel,:),1,'omitnan')'; 
        MM2= mean( M2(fsel,:),1, 'omitnan')';

        figure('COlor','w','Position', [ 300 400 300 500],'Renderer', 'painters')
        
        % Clean up NaNs from data
        Base_data = MM;
        Base_data(isnan(Base_data)) = [];
        Stim_data = MM1;
        Stim_data(isnan(Stim_data)) = [];
        
        num_bases = length(Base_data);
        num_stims = length(Stim_data);

        labels = [repmat({'Base'}, num_bases, 1); repmat({'Stim'}, num_stims, 1)];

        data = [Base_data; Stim_data];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        %line([ 0.7 2.7], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
        %xlim([ 0.9 2.7]); %ylim([ -0.1 0.7])
        title(['Violin stim frequency filtered spike-vm ' DBS_reg{area} ' ' DBS_cond{cond}], 'Interpreter', 'none');

        set(gca, 'Units', 'centimeters', 'Position', [3, 3, 3.214, 5.826], 'PositionConstraint', 'innerposition');

        print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'Violin_stimFrq_spike_vm_plV_' DBS_reg{area}(3:end) '_cond_' DBS_cond{cond}(3:end) '.pdf'])
        diary on;
        disp(['Violin stim frequency spike-vm ' DBS_reg{area} ' ' DBS_cond{cond}]);
        [p,h,stats] = signrank(MM, MM1)
        diary off;
    end
end

%%
 
 if 0
   cond=1
  M= [allPLVb2{1,cond},allPLVb2{2,cond}] ;
  M1= [allPLVs2{1,cond},allPLVs2{2,cond}] ;
    M2= [allPLVp2{1,cond},allPLVp2{2,cond}] ;  
    
     figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
     plot(frs, abs(nanmean(M,2)),'b','Linewidth',1); hold on,
         plot(frs, abs(nanmean(M1,2)),'g','Linewidth',1); 
  % plot(frs, abs(nanmean(M2,2)),'m','Linewidth',1); 
 axis tight  
 end
  
%% Spike-Vm PLV for theta (but at each individual frequency, not broad band filtered)
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
            axis tight  ; 
            %ylim([ -0.02 0.07]); 
            xlim([ 2 10])
            % set(gca,'Xscale','log')
            title([DBS_reg{area}(3:end) ' ' DBS_cond{cond}(3:end)], 'Interpreter', 'none');
            print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'PLV_THETA' DBS_reg{area}(3:end) '_area_high_' DBS_cond{cond}(3:end) '.pdf'])
    end
end
    
%% Violin plots for spike-vm phase for theta range
stats_log = [savepath 'PLV' f 'Theta_spike_vm_violin_stats'];
if exist(stats_log), delete(sprintf('%s', stats_log)), end;
diary(stats_log);
diary off

for area=1:2%:2;
    for cond=1:2;
        %M= [allPLVb{area,1} allPLVb{area,2}] ; M1= [allPLVs{area,1},allPLVs{area,2}] ;M2= [allPLVp{area,1} allPLVp{area,2}] ;
        M = [allPLVb{area,cond}];
        M1 = [allPLVs{area,cond}];
        fsel=frs>Multi_func.theta_range(1) & frs< Multi_func.theta_range(2);  
        MM= mean( M(fsel,:),1, 'omitnan')';
        MM1= mean( M1(fsel,:),1, 'omitnan')'; 
        
        % Remove Nans because they are there for some reason

        num_neurons = size(MM, 1);
        
        figure('COlor','w','Position', [ 300 400 300 500],'Renderer', 'painters')
        Base_data = MM;
        Base_data(isnan(Base_data)) = [];
        Stim_data = MM1;
        Stim_data(isnan(Stim_data)) = [];
        
        num_bases = length(Base_data);
        num_stims = length(Stim_data);

        labels = [repmat({'Base'}, num_bases, 1); repmat({'Stim'}, num_stims, 1)];

        data = [Base_data; Stim_data];
        
        ViolinOpts = Multi_func.get_default_violin();
        violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);
        violins(1).ViolinColor = {Multi_func.base_color};
        violins(2).ViolinColor = {Multi_func.stim_color};
        
        %line([ 0.7 2.7], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
        %xlim([ 0.9 2.7]); %ylim([ -0.1 0.7])
        title(['Violin theta spike-vm ' DBS_reg{area} ' ' DBS_cond{cond}], 'Interpreter', 'none');
        
        set(gca, 'Units', 'centimeters', 'Position', [3, 3, 3.214, 5.826], 'PositionConstraint', 'innerposition');
        print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'Violin_THETA_' DBS_reg{area} '_cond_' DBS_cond{cond} '.pdf']);
        diary on;
        disp(['Violin theta spike-vm ' DBS_reg{area} ' ' DBS_cond{cond}]);
        [p,h,stats] = signrank(MM, MM1)
        diary off;
    end
end

%%
      
           M=  [allPLVb{1,1},allPLVb{1,2} allPLVb{2,1},allPLVb{2,2} ] ; 
        M1=[allPLVs{1,1},allPLVs{1,2} allPLVs{2,1},allPLVs{2,2} ] ; 
                M2=[allPLVp{1,1},allPLVp{1,2} allPLVp{2,1},allPLVp{2,2} ] ; 

        %  M2= [allPLVp{1,cond},allPLVp{2,cond}] ;
        figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters')
        fill_error_area2(frs, nanmean(M,2),nanstd(M,[],2)./sqrt(sum(~isnan((M(1,:))))  ) ,[ .5 .5 .5] )
        fill_error_area2(frs, nanmean(M1,2),nanstd(M1,[],2)./sqrt(sum(~isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
        fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(~isnan((M2(1,:))))  ) ,[ .5 .5 .5] )
        plot(frs,nanmean(M,2),'b','Linewidth',1); hold on; plot(frs,nanmean(M1,2),'g','Linewidth',1)
          plot(frs,nanmean(M2,2),'m','Linewidth',1)
        axis tight  ; ylim([ -0.02 0.07]); 
       xlim([ 1 10])
        % set(gca,'Xscale','log')
        print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_THETA_COMB_ALL_''.pdf'])
    
    
        


function  [filt_sig]=filt_data(sig,frs, FS)
Fn = FS/2;

for steps=1:length(frs);


FB=[ frs(steps)*0.97 frs(steps)*1.03];

 [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
 filt_sig(:,steps)= hilbert(filtfilt(B,A,sig));
            
 
end
end



function [ PLV_output,PLV_output2]= spike_field_ppc_Pierre(wavD,spikes,timshift)
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
     M=exp(1i.*((squeeze(wavD(:,:,ind))))); end
      s=find(spikes{ind});s=s-timshift; s(s<=0)=[];
      
  mat=[mat,    M(s,:)'];
end

NT= sum(~isnan(mat(1,:)));

  Z=abs(nanmean(mat,2));  % MVL
   T=   Z.^2; 
PLV_output= (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
PLV_output2= mat; 
if NT<7
PLV_output=PLV_output.*NaN;
%PLV_output2=PLV_output2.*NaN;
end

end
