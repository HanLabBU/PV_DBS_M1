clear all
 
f = filesep;

% Linux computer paths
server_root_path = '~/handata_server/';
local_root_path = '~/Projects/';
pv_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
addpath([local_root_path f 'Pierre Fabris' f 'PV DBS neocortex' f 'Scripts' f]);

savepath = Multi_func.save_plot;


switch_Source=2; % 1=DBS, 2=Vm


for id=1:2; for id2=1:3; allPLVb2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVs2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVp2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allVmDBS{id,id2}=[];end;end

ses = dir([pv_path '*.mat']);
all_matfiles = {ses.name}
          
[region_matfiles] = Multi_func.find_region(all_matfiles);        
all_matfiles = {ses.name};
gg= fieldnames(region_matfiles)' ;   
region=0;

%% loop through each condition and calculate the PLV
for f_region = fieldnames(region_matfiles)' %gg(2)% 
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
    DBS_cond= {'f_40'; 'f_140' ; 'f_1000'}
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
            frs=[2:1:10 11:1:150 ];FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([ matfile{1}]);
            clear alls wavD  allVm
            SNR=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
                if ~isempty(  data.align.trial{tr}) 
                    if mean(data.align.trial{tr}.spike_info375.spike_snr{1}) > 4 % Trial included if spike SNR above 4SNR
                        tr_iter= tr_iter+1;
                        voltsig=data.align.trial{tr}.spike_info375.trace_ws;
                        [ fitbaseline, coeff]=Multi_func.exp_fit_Fx(voltsig',FS); %remove photobleaching
                        voltsig=(voltsig-fitbaseline);
                        
                        if  1
                            Camtime=  data.raw.trial{tr}.raw_camera_frame_time- data.raw.trial{tr}.raw_camera_start_time(1);
                            DBStim=data.raw.trial{tr}.raw_stimulation_time- data.raw.trial{tr}.raw_camera_start_time(1);
                            
                            traceDBS=zeros(1,length(data.raw.trial{tr}.raw_traces));
                            for id=1:length(DBStim)
                                z=find(Camtime > DBStim(id));try
                                    traceDBS(z(1))=1;end
                            end
                        end
                        if switch_Source==1
                        voltsig=traceDBS;end
                        allVm(:, tr_iter)=voltsig;
                        [filt_sig]=Multi_func.filt_data(voltsig,frs, FS); % filtering+Hilbert
                        wt=filt_sig';
                         %   [wt, f,f1] = cwt(voltsig, round(FS),'TimeBandwidth',100);
                          %  frs=f;
                    
                        
                        wavD(:,:, tr_iter)=wt';
                        SNR = [SNR;(data.align.trial{tr}.spike_info375.spike_snr{1})];
                        alls(:, tr_iter)=data.align.trial{tr}.spike_info375.roaster;
                   end 
               end
            end
       
            if tr_iter   >0 
                baseTimsel= [10:FS-10];
                StimTimsel= [FS:FS*2];
                PostTimsel= [FS*2+10:size(alls,1)];
   
                if  mean(SNR) >4  % mean spike SNR for the neuron to be over 4SNR (can be changed)
                    % Read notes, account for spike interference in Vm power value
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
    
    
%%
pheight=130;

% mm=0;  figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
% for region=1:2;
%     for cond=1:2
%         mm=mm+1;
%         subplot(2,2,mm)
%         plot(frs,nanmean(allPLVb{region,cond},2),'b','Linewidth',1); hold on
%         plot(frs,nanmean(allPLVs{region,cond},2),'r','Linewidth',1)
%         %plot(frs,nanmean(allPLVp{region,cond},2),'m','Linewidth',1)
%         axis tight;
%         set(gca,'Xscale','log')
%         ylim([ -0.05 0.33])
%     end
% end
%     
% for cond=1:2
%     M= [allPLVb{1,cond},allPLVb{1,cond}] ; M1= [allPLVs{1,cond},allPLVs{1,cond}] ;
%     M2= [allPLVp{1,cond},allPLVp{2,cond}] ;
%     figure('COlor','w','Position', [ 300 400 120 pheight],'Renderer', 'painters') 
%     fill_error_area2(frs, nanmean(M,2),nanstd(M,[],2)./sqrt(sum(isnan((M(1,:))))  ) ,[ .5 .5 .5] )
%     fill_error_area2(frs, nanmean(M1,2),nanstd(M1,[],2)./sqrt(sum(isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
%     %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] ) 
%     plot(frs,nanmean(M,2),'b','Linewidth',1); hold on; plot(frs,nanmean(M1,2),'g','Linewidth',1)
%     %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
%     axis tight
%     ylim([ -0.05 0.2]); xlim([ 2 20])
%     % set(gca,'Xscale','log')
%     print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_across_area_low_' num2str(cond) '.pdf'])
% end
% 
% % ACROSS AREAS  THETA
% for cond=2
%     M= [allPLVb{1,cond},allPLVb{2,cond},allPLVb{1,2},allPLVb{2,2},] ;
%     M1= [allPLVs{1,cond},allPLVs{2,cond},allPLVs{1,2},allPLVs{2,2}] ;
%     M2= [allPLVp{1,cond},allPLVp{2,cond}] ;M(M==0)=NaN ;M1(M1==0)=NaN;
%     figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
%     fill_error_area2(frs, fastsmooth(nanmean(M,2),1,1,1),nanstd(M,[],2)./sqrt(sum(isnan((M(1,:))))  ) ,[ .5 .5 .5] )
%     fill_error_area2(frs, fastsmooth(nanmean(M1,2),1,1,1),nanstd(M1,[],2)./sqrt(sum(isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
%     %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] ) 
%     plot(frs,fastsmooth(nanmean(M,2),1,1,1),'b','Linewidth',1); hold on; plot(frs, fastsmooth(nanmean(M1,2),1,1,1),'g','Linewidth',1)
%     %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
%     axis tight  ; ylim([ -0.05 0.12]); xlim([ 1 30])
%     % set(gca,'Xscale','log')
%   %  print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_across_area_' num2str(switch_Source) '_' num2str(cond) '.pdf'])
% end
    
%% PLV combined between both regions
sm=5;
% ACROSS AREAS
for cond=1:2
    % PLV baseline
    M= [allPLVb{1,cond},allPLVb{2,cond}] ; 
    %PLV stimulation 
    M1= [allPLVs{1,cond},allPLVs{2,cond}] ;
    M2= [allPLVp{1,cond},allPLVp{2,cond}] ;M(M==0)=NaN ;M1(M1==0)=NaN;
    figure('COlor','w','Position', [ 300 400 180 pheight],'Renderer', 'painters') 
    fill_error_area2(frs, fastsmooth(nanmean(M,2),sm,1,1),nanstd(M,[],2)./sqrt(sum(isnan((M(1,:))))  ) ,[ .5 .5 .5] )
    fill_error_area2(frs, fastsmooth(nanmean(M1,2),sm,1,1),nanstd(M1,[],2)./sqrt(sum(isnan((M1(1,:))))  ) ,[ .5 .5 .5] )
    %fill_error_area2(frs, nanmean(M2,2),nanstd(M2,[],2)./sqrt(sum(isnan((M2(1,:))))  ) ,[ .5 .5 .5] ) 
    plot(frs,fastsmooth(nanmean(M,2),sm,1,1),'b','Linewidth',1); hold on; plot(frs,fastsmooth(nanmean(M1,2),sm,1,1),'g','Linewidth',1)
    %  plot(frs,nanmean(M2,2),'m','Linewidth',1)
    axis tight  ; ylim([ -0.05 0.3]);% xlim([ 21 170])
    % set(gca,'Xscale','log')
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'PLV_across_area_' num2str(switch_Source) '_' num2str(cond) '.pdf'])
end
   

%% Violin plots PLV at theta frequency combining all regions
cond=1
M= [allPLVb{1,cond},allPLVb{2,cond},allPLVb{1,2},allPLVb{2,2}] ; 
M1= [allPLVs{1,cond},allPLVs{2,cond},allPLVs{1,2},allPLVs{2,2}] ;
fsel=frs>2 & frs< 12;   
M(M==0)=NaN;
M1(M1==0)=NaN;
MM1= nanmean( M(fsel,:),1)'; 
MM= nanmean( M1(fsel,:),1)';
num_neurons = size(MM1, 1);

figure('COlor','w','Position', [ 300 400 150 120],'Renderer', 'painters')
data = [MM1; MM];
labels = [repmat({'Base'}, num_neurons, 1); repmat({'Stim'}, num_neurons, 1)];
ViolinOpts = Multi_func.get_default_violin();
violins = violinplot(data, labels, 'GroupOrder', {'Base', 'Stim'}, ViolinOpts);

violins(1).ViolinColor = {Multi_func.base_color};
violins(2).ViolinColor = {Multi_func.stim_color};

% Original violin plot codes
%violinplot2(MM1,[1.3 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
%hold on,
%violinplot2(MM,[ 1.8],'ViolinColor', [ 0.5 0.8 0; 0.9 0 0])

line([ 0.7 2.3], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
xlim([ 0.9 2.1])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Violin_theta_' num2str(switch_Source) '_' num2str(cond) '.pdf'])
[h,p,ci,stats] =ttest2(MM,MM1)

%% Violin PLV combined area at respective stimulation frequency range
cond=1
M= [allPLVb{1,cond},allPLVb{2,cond}]; 
M1= [allPLVs{1,cond},allPLVs{2,cond}] ;
fsel=frs>38 & frs< 42;  
MM1= nanmean( M(fsel,:),1)'; 
MM= nanmean( M1(fsel,:),1)';

figure('COlor','w','Position', [ 300 400 150 120],'Renderer', 'painters')

violinplot2(MM1,[1.3 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
hold on,
violinplot2(MM,[ 1.8],'ViolinColor', [ 0.5 0.8 0; 0.9 0 0])
line([ 0.7 2.3], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
xlim([ 0.9 2.1]); ylim([ -0.2 0.9])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Violin_gamma_' num2str(switch_Source) '_' num2str(cond) '.pdf'])
[h,p,ci,stats] =ttest2(MM,MM1)

cond=2
M= [allPLVb{1,cond},allPLVb{2,cond}] ; M1= [allPLVs{1,cond},allPLVs{2,cond}] ;
fsel=frs>138 & frs< 142;  
MM1= nanmean( M(fsel,:),1)'; MM= nanmean( M1(fsel,:),1)';
figure('COlor','w','Position', [ 300 400 150 120],'Renderer', 'painters')
violinplot2(MM1,[1.3 ],'ViolinColor', [ 0 0 0.9; 0 0 0.9])
hold on,
violinplot2(MM,[ 1.8],'ViolinColor', [ 0.5 0.8 0; 0.9 0 0])
line([ 0.7 2.3], [ 0  0],'COlor', [ 1 1 1 ],'Linewidth',0.5)
xlim([ 0.9 2.1]); ylim([ -0.2 0.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Violin_140_' num2str(switch_Source) '_' num2str(cond) '.pdf'])
[h,p,ci,stats] =ttest2(MM,MM1)



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


FB=[ frs(steps)*0.95 frs(steps)*1.05];

 [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
 filt_sig(:,steps)= hilbert(filtfilt(B,A,sig));
            
 
end
end

function [ fitbaseline, coeff]=exp_fit_Fx(v,FS)
    %% fits an expoential to estimate the photobleaching
    v1=v;
    v1([FS-10:2000])= mean([ v([FS-30:FS-10 2000:2020 ])]);
    v1(1)=v1(2);
    F = @(x,xdata)x(1)+x(2)*exp(- xdata./x(3));%+ x(3)*exp(- xdata./x(4))  ;
    x0 = [mean(v1) 40 1.5   ] ;
    OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
    t = (1:length(v))./FS;
    tsel=1:length(t);
    [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);
    fitbaseline=xunc(1)+xunc(2)*exp(-t./xunc(3));
    coeff=xunc;
% figure,plot(t,fit_baseline)
%       hold on,plot(t,v)
end



%% This function does not seem to be called in this script
%function [ PLV_output,PLV_output2]= spike_field_ppc_Pierre(wavD,spikes,timshift)
%    if size(wavD,4)>1
%        [tr chs frs t]= size(wavD);
%    else
%        [t frs tr]= size(wavD);
%    end
%    mat=[];
%    for ind=1:tr
%        % Question: when is the size of wavD different?
%        if size(wavD,4)>1
%            M=exp(1i.*((squeeze(wavD(ind,CH,:,:))))); 
%        else
%            M=exp(1i.*((squeeze(wavD(:,:,ind))))); end
%            s=find(spikes{ind});
%            s=s-timshift; s(s<=0)=[];
%            mat=[mat,    M(s,:)'];
%        end
%    
%        NT= sum(~isnan(mat(1,:)));
%        Z=abs(nanmean(mat,2));  % MVL
%        T=   Z.^2; 
%        PLV_output= (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
%        PLV_output2= mat; 
%    if NT<10
%        PLV_output=PLV_output.*NaN;
%        %PLV_output2=PLV_output2.*NaN;
%    end
%end
