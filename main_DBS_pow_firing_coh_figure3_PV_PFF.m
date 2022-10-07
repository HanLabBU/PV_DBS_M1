addpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\DBS_Scripts\')

% power, entrainment


%%
clear all

 stim_type=[];
 % 1= 140Hz, other 40Hz
allsp=[]; allVm=[]; firing_conc=[];Vm_conc=[];
  allsamp=[];stim_type_tr=[];
  allangs40=[];allangs140=[];stim_type_sp=[];
mr=0;

for stim_freq=[ 0 1]

    if stim_freq==1
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\PV\140\')
    else
        cd('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\PV\40\')
    end
    
    
    ses=dir('*.mat');
    
    for ind=1:length(ses)
        Cpath= ses(ind).name;
        load(Cpath)  %loading
        
        trial_numb=unique(result.trial_vec);
        
        if length(find(result.trial_vec==trial_numb(1) )) <3500 & length(trial_numb)>2
        
            %%  denoise
            FS=828;
            Fn = FS/2;FB=[ 86.5 88.5];
            [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
            LFPg= ((filtfilt(B,A,    result.traces(:,1))));
            %result.traces(:,1)= result.traces(:,1)-LFPg; % Camera noise filtering
            %%
            lfp=[];
            clear allV allS alls1  allV40 allV140
            tr=0;
        
            for  ne=unique(result.trial_vec) % trials
                tr=tr+1;
                %% Vm
                v= result.traces(result.trial_vec==ne,1)  ;
                v=v(1:2500);
                
                vsub= result.resultS{ne}.trace_ws(1,:)'  ;
                vsub=vsub(1:2500);
            
                FS=828;
                Fn = FS/2;FB=[ 86.5 88.5];
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                LFPg= ((filtfilt(B,A,   vsub)));
                vsub=vsub-LFPg;
           
                Fn = FS/2;FB=[ 35 45];
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                Vsub40= angle(hilbert(filtfilt(B,A,     vsub)));
                Fn = FS/2;FB=[ 135 145];
                [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
                Vsub140= angle(hilbert(filtfilt(B,A,     vsub)));
            
                %v=v-mean(v);
                
                %  v=v-fastsmooth(v,2300,1,1);
          
                % exponential fitting to remove photobleaching
                [ fitbaseline, coeff]=exp_fit_Fx(v,FS);
      
                v=(v-fitbaseline');
                % v(1:10)=NaN;
                lfp.trial{tr}= ( vsub)';
                lfp.time{tr}= (1:size(v,1))./FS;
                allV(1:length(v),ne)=v;
                allV40(1:length(v),ne)=Vsub40;
                allV140(1:length(v),ne)=Vsub140;
            
                %%spikes
                strain= result.resultS{ne}.roaster(1:2500);
                allS(1:length(strain),tr)=strain;
                sid=result.resultS{ne}.spike_idx{1};sid=sid;
                sid(sid<10)=[];
                vect=zeros(1,size(allS,1));vect(sid)=1;
                alls1( :,tr)=vect(1:2500);
                spx= result.resultS{ne}.spike_idx{1}   ;
                samp= result.resultS{ne}.spike_amplitude{1};
                samp(spx <10)=NaN;
         
                allsamp=[allsamp ,samp];
            end
            
            for id=1%:size(vsig,2)
                lfp.label{id}= num2str(id);
            end  
            
            []; %block_type == cfg.blk
            cfg.method ='mtmconvol'; %'mvar';
            cfg.output ='fourier';
            cfg.taper='hanning';
            cfg.keeptapers ='yes';
            cfg.keeptrials ='yes';
            cfg.trials='all';cfg.tapsmofrq =5;%
            cfg.channel= [1]%; %chans=cfg.channel;
            cfg.foi= [2:1:150];
            cfg.toi=lfp.time{1}(1:1:end) ;
            cfg.width =5;
            cfg.t_ftimwin =[ones(1,length(cfg.foi))*0.5];
            freq2 = ft_freqanalysis(cfg, lfp);
      
            wavD = angle(squeeze(freq2.fourierspctrm));
            wavA = abs(squeeze(freq2.fourierspctrm));
    
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            mr=mr+1;
            allname{mr}=  Cpath;
    
            allpow(1:size(wavA,3),:,mr)= squeeze(nanmean(wavA,1))';
    
            allfiring(:,mr)= nanmean(allS.*FS,2);
            
            sm=10 ;%smoothing parameter (rectangular window)
            allfiringSM(:,mr)= nanfastsmooth(nanmean(allS.*FS,2),sm,1,1);
            allfiringSM(:,mr)= allfiringSM(:,mr)- nanmean(allfiringSM(500:800,mr));
            allVm(:,mr)= nanmean(allV,2)./nanmean(samp);
     
            sm=50 ;%smoothing parameter (rectangular window)
            allVmSM(:,mr)= nanfastsmooth(nanmean(allV,2),sm,1,1)./nanmean(samp);
            allVmSM(:,mr)= allVmSM(:,mr)- nanmean(allVmSM(500:800,mr));
            firing_conc=[ firing_conc,allS];
            Vm_conc=[Vm_conc,bsxfun(@rdivide,allV, nanmean(samp))];   
            stim_type=[stim_type, stim_freq];   
            stim_type_tr=[ stim_type_tr, ones(1,size(allS,2)).*stim_freq];
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
            baseTimsel= [250:750];
            StimTimsel= [829:1656];
            PostTimsel= [1690:2490];
    
            sta_base=STA_computeDBS(lfp.trial, alls1,1,60,60,0);
            sta_stim=STA_computeDBS(lfp.trial, alls1,1,60,60,1);
    
            % figure,plot(sta_stim)
    
            clear s
            M= wavD(:,:,baseTimsel); 
            
            for tr1=1:size(allS,2); 
                s{tr1}=alls1(baseTimsel,tr1)'; 
            end
                
            [ PLV_b]= spike_field_ppcDBS(M,s,1);
            clear s
            M= wavD(:,:,StimTimsel); 
      
            for tr1=1:size(allS,2); 
                s{tr1}=alls1(StimTimsel,tr1)'; 
            end        
    
            [ PLV_s]= spike_field_ppcDBS(M,s,1);
            M= wavD(:,:,PostTimsel); 
      
            for tr1=1:size(allS,2); 
                s{tr1}=alls1(PostTimsel,tr1)'; 
            end
            
            [ PLV_p]= spike_field_ppcDBS(M,s,1);
    
            allPLVb(:,mr)=PLV_b;
            allPLVs(:,mr)=PLV_s;
            allPLVp(:,mr)=PLV_p;
       
            %%%%%%%%%%%%%%
            %% PLV with filtered signals
            allS1= alls1;
            allS1([1:800 1650:size(allS1,2)],:)=0;
            
            allangs140= [allangs140; allV140(allS1==1)];
            allangs40= [allangs40; allV40(allS1==1)];
        
            stim_type_sp=[ stim_type_sp, ones(1,length(find(allS1==1))).*stim_freq];
            Z=abs(nanmean(exp(1i.*allV140(allS1==1))));
            Za=angle(nanmean(exp(1i.*allV140(allS1==1))));NT=sum(sum(allS1==1));
            if NT>5; 
                allPLVf140(mr)=Z; 
                allPLVf140a(mr)=Za;
            else 
                allPLVf140(mr)=NaN;
                allPLVf140a(mr)=NaN;
            end
           
            Z=abs(nanmean(exp(1i.*allV40(allS1==1))));
            Za=angle(nanmean(exp(1i.*allV40(allS1==1))));NT=sum(sum(allS1==1)); 
            
            if NT>5; 
                allPLVf40(mr)=Z; 
                allPLVf40a(mr)=Za;
            else 
                allPLVf40(mr)=NaN;
                allPLVf40a(mr)=NaN;
            end
         
            Z=abs(nanmean(exp(1i.*allV40(allS1(randperm(length(allS1)),:)==1))));
            NT=sum(sum(allS1==1)); 
          
            if NT>5; 
                allPLVfSH40(mr)=Z;
            else 
                allPLVfSH40(mr)=NaN;
            end
          
            Z=abs(nanmean(exp(1i.*allV140(allS1(randperm(length(allS1)),:)==1))));NT=sum(sum(allS1==1)); 
          
            if NT>5; 
                allPLVfSH140(mr)=Z;
            else 
                allPLVfSH140(mr)=NaN;
            end
        end  
    end % sessions
end % stimulation type

%allS(allS==0)=NaN;
allPLVb( allPLVb==0)=NaN;
allPLVs( allPLVs==0)=NaN;
allPLVp( allPLVp==0)=NaN;

baseTimsel= [150:750];
StimTimsel= [828:1600];
StimTimselTR= [828:1000];
StimTimselSU= [1000:1600];
PostTimsel= [1690:2490];

tim_axis=([1:size(allV,1)]-(FS-10))./FS;

%%%%
savepath='\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\PV_figure\'
pheight=150;
%%%%%%%%%%%%%%%%%%

if 0
    V1=nanmean(allVm(:, stim_type==0),2);
    V1s=nanstd(allVm(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
    V2=nanmean(allVm(:, stim_type==1),2);
    V2s=nanstd(allVm(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));
    %% NOn smoothed Vm , overlay 40 and 140
    figure('COlor','w','Position',[300 300 250 pheight],'Renderer', 'painters'),
    plot(tim_axis,V1,'b','Linewidth',1.5)
    hold on,plot(tim_axis,V2,'r','Linewidth',1.5)
    axis tight;xlim([-0.7 1.7])
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Nonsmooth_Vm_av_40_140t.pdf'])
    savefig(gcf, [ savepath 'Nonsmooth_Vm_av_40_140.fig'])
    
    %% Vm average plots (smoothed)
    V1=nanmean(allVmSM(:, stim_type==0),2);
    V1s=nanstd(allVmSM(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
    V2=nanmean(allVmSM(:, stim_type==1),2);
    V2s=nanstd(allVmSM(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));

    figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
    plot(tim_axis,V1,'k','Linewidth',1.5)
    fill_error_area2(tim_axis,V1,V1s, [ 0.5 0.5 0.5]);
    axis tight;; xlim([-0.7 1.7]);ylim([-0.3 1.2])
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'smooth_Vm_av_40.pdf'])
    savefig(gcf, [ savepath 'smooth_Vm_av_40.fig'])

    figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
    plot(tim_axis,V2,'k','Linewidth',1.5)
    fill_error_area2(tim_axis,V2,V2s, [ 0.5 0.5 0.5]);
    axis tight; xlim([-0.7 1.7]);ylim([-0.3 1.2])
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'smooth_Vm_av_140.pdf'])
    savefig(gcf, [ savepath 'smooth_Vm_av_140.fig'])

    %% Vm bar plots /quantifications
    
    V1b=nanmean(allVmSM(baseTimsel, stim_type==0),1);
    V2b=nanmean(allVmSM(baseTimsel, stim_type==1),1);
    V1s=nanmean(allVmSM(StimTimselTR, stim_type==0),1);
    V2s=nanmean(allVmSM(StimTimselTR, stim_type==1),1);
    V1s2=nanmean(allVmSM(StimTimselSU, stim_type==0),1);
    V2s2=nanmean(allVmSM(StimTimselSU, stim_type==1),1);
    V1p=nanmean(allVmSM(PostTimsel, stim_type==0),1);
    V2p=nanmean(allVmSM(PostTimsel, stim_type==1),1);
    [h,p,ci,stats] = ttest(V1b, V1s)
    [h,p,ci,stats] = ttest(V2b, V2s)
  
    figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
    M=[V1b',V1s', V1s2',V1p'];M2=[ V2b',V2s', V2s2',V2p'];
    b1=bar(1,nanmean(V1b),'Facecolor',[ 0 0 0.4]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(2,nanmean(V1s),'Facecolor',[ 0 0 0.9]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(3,nanmean(V1s2),'Facecolor',[ 0 0 0.9]);hold on,
    set(b1,'FaceAlpha',0.4)
    b1=bar(4,nanmean(V1p),'Facecolor',[ 0  0 0.2]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(6,nanmean(V2b),'Facecolor',[ 0.4 0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(7,nanmean(V2s),'Facecolor',[ 0.9 0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(8,nanmean(V2s2),'Facecolor',[ 0.9 0 0]);hold on,
    set(b1,'FaceAlpha',0.4)
    b1=bar(9,nanmean(V2p),'Facecolor',[ 0.2  0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    errorbar([1 2 3 4],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
    errorbar([6 7 8 9],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Vm_av_barplot.pdf'])
    savefig(gcf, [ savepath 'Vm_av_barplot.fig'])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% FIRING RATE %%
    
    
    V1=nanmean(allfiringSM(:, stim_type==0),2);
    V1s=nanstd(allfiringSM(:, stim_type==0),[],2)./sqrt(length(find(stim_type==0)));
    V2=nanmean(allfiringSM(:, stim_type==1),2);
    V2s=nanstd(allfiringSM(:, stim_type==1),[],2)./sqrt(length(find(stim_type==1)));
    SPM=firing_conc(:,stim_type_tr==0)==1;
    figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
    plot(tim_axis,V1,'k','Linewidth',1)
    fill_error_area2(tim_axis,V1,V1s, [ 0.5 0.5 0.5]);
    axis tight;; xlim([-0.7 1.7]);%ylim([-0.3 1.2])
    % for ind=1:size(SPM,2)
    %    plot(tim_axis(SPM(:,ind)),  ones(1,length(find(SPM(:,ind))))*ind+32,'r.');  hold on,
    % end
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_40_sm.pdf'])
    savefig(gcf, [ savepath 'Firing_40_sm.fig'])

    figure('COlor','w','Position',[300 300 300 pheight],'Renderer', 'painters'),
    plot(tim_axis,V2,'k','Linewidth',1.5)
    fill_error_area2(tim_axis,V2,V2s, [ 0.5 0.5 0.5]);
    axis tight; xlim([-0.7 1.7]);%ylim([-0.3 1.2])
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_140_sm.pdf'])
    savefig(gcf, [ savepath 'Firing_140_sm.fig'])

    %% BARPLOT FIRING
    V1b=nanmean(allfiringSM(baseTimsel, stim_type==0),1);
    V2b=nanmean(allfiringSM(baseTimsel, stim_type==1),1);
    V1s=nanmean(allfiringSM(StimTimselTR, stim_type==0),1);
    V2s=nanmean(allfiringSM(StimTimselTR, stim_type==1),1);
    V1s2=nanmean(allfiringSM(StimTimselSU, stim_type==0),1);
    V2s2=nanmean(allfiringSM(StimTimselSU, stim_type==1),1);
    V1p=nanmean(allfiringSM(PostTimsel, stim_type==0),1);
    V2p=nanmean(allfiringSM(PostTimsel, stim_type==1),1);
    
    %%
    %% Base vs trans
    [h,p,ci,stats] = ttest(V1b, V1s) % 40Hz
    [h,p,ci,stats] = ttest(V2b, V2s)% 140Hz
    %% Base vs sustained
    [h,p,ci,stats] = ttest(V1b, V1s2) % 40Hz
    [h,p,ci,stats] = ttest(V2b, V2s2)% 140Hz
    %% Base vs post
    [h,p,ci,stats] = ttest(V1b, V1p) % 40Hz
    [h,p,ci,stats] = ttest(V2b, V2p)% 140Hz
    %% trans between DBS cond
    [h,p,ci,stats] = ttest2(V2s-V2b, V1s-V1b) 
    %% sust between DBS cond
    [h,p,ci,stats] = ttest2(V2s2-V2b, V1s2-V1b) 
    %% post between DBS cond
    [h,p,ci,stats] = ttest2(V2p-V2b, V1p-V1b) 
    %%
  
    figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
    M=[V1b',V1s', V1s2',V1p'];M2=[ V2b',V2s', V2s2',V2p'];
    b1=bar(1,nanmean(V1b),'Facecolor',[ 0 0 0.4]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(2,nanmean(V1s),'Facecolor',[ 0 0 0.9]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(3,nanmean(V1s2),'Facecolor',[ 0 0 0.9]);hold on,
    set(b1,'FaceAlpha',0.4)
    b1=bar(4,nanmean(V1p),'Facecolor',[ 0  0 0.2]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(6,nanmean(V2b),'Facecolor',[ 0.4 0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(7,nanmean(V2s),'Facecolor',[ 0.9 0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    b1=bar(8,nanmean(V2s2),'Facecolor',[ 0.9 0 0]);hold on,
    set(b1,'FaceAlpha',0.4)
    b1=bar(9,nanmean(V2p),'Facecolor',[ 0.2  0 0]);hold on,
    set(b1,'FaceAlpha',0.7)
    errorbar([1 2 3 4],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
    errorbar([6 7 8 9],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
    print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Firing_barplot.pdf'])
    savefig(gcf, [ savepath 'Firing_barplot.fig'])
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%  Spectral power %%%%%%%%%%%%

Pow_40=nanmedian(allpow(:,:,stim_type==0),3);
%Pow_40=bsxfun(@rdivide, Pow_40, nanmean(Pow_40(baseTimsel,:)));
Pow_140=nanmedian(allpow(:,:,stim_type==1),3);

figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
imagesc(tim_axis,freq2.freq,(Pow_40.*repmat(freq2.freq.^0.5,length(tim_axis),1)  )');axis xy
colormap(jet); xlim([-0.7 1.7])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'TFR_40.pdf'])
savefig(gcf, [ savepath 'TFR_40.fig'])

figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
imagesc(tim_axis,freq2.freq,(Pow_140.*repmat(freq2.freq.^0.6,length(tim_axis),1)  )');axis xy
colormap(jet); xlim([-0.7 1.7])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'TFR_140.pdf'])
savefig(gcf, [ savepath 'TFR_140.fig'])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

baseF=find(tim_axis>-0.5 & tim_axis<0);

figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
imagesc(tim_axis,freq2.freq,bsxfun(@rdivide, Pow_40,nanmean(Pow_40(baseF,:),1)   )'-1);axis xy
colormap(jet); xlim([-0.7 1.7])
set(gca,'Clim',[-0.5 0.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'TFR_40REL.pdf'])


figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
imagesc(tim_axis,freq2.freq,bsxfun(@rdivide, Pow_140,nanmean(Pow_140(baseF,:),1)   )'-1);axis xy
colormap(jet); xlim([-0.7 1.7])
set(gca,'Clim',[-0.5 0.5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'TFR_140REL.pdf'])





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allpow2=bsxfun(@times, allpow, freq2.freq.^0);
fsel=find(freq2.freq>36  & freq2.freq<44); 
fsel3=find(freq2.freq>100 ); 
N1b=squeeze(nanmean(nanmean(allpow2(baseTimsel, fsel,stim_type==0),1),2));
N2b=squeeze(nanmean(nanmean(allpow2(baseTimsel, fsel3,stim_type==0),1),2));
N1s=squeeze(nanmean(nanmean(allpow2(StimTimselTR,fsel, stim_type==0),1),2));
N2s=squeeze(nanmean(nanmean(allpow2(StimTimselTR,fsel3, stim_type==0),1),2));
N1s2=squeeze(nanmean(nanmean(allpow2(StimTimselSU,fsel, stim_type==0),1),2));
N2s2=squeeze(nanmean(nanmean(allpow2(StimTimselSU,fsel3, stim_type==0),1),2));
N1p=squeeze(nanmean(nanmean(allpow2(PostTimsel,fsel, stim_type==0),1),2));
N2p=squeeze(nanmean(nanmean(allpow2(PostTimsel,fsel3, stim_type==0),1),2));
fsel=find(freq2.freq>136  & freq2.freq<144); 
V1b=squeeze(nanmean(nanmean(allpow2(baseTimsel, fsel3,stim_type==1),1),2));
V2b=squeeze(nanmean(nanmean(allpow2(baseTimsel, fsel,stim_type==1),1),2));
V1s=squeeze(nanmean(nanmean(allpow2(StimTimselTR,fsel3, stim_type==1),1),2));
V2s=squeeze(nanmean(nanmean(allpow2(StimTimselTR,fsel, stim_type==1),1),2));
V1s2=squeeze(nanmean(nanmean(allpow2(StimTimselSU,fsel3, stim_type==1),1),2));
V2s2=squeeze(nanmean(nanmean(allpow2(StimTimselSU,fsel, stim_type==1),1),2));
V1p=squeeze(nanmean(nanmean(allpow2(PostTimsel,fsel3, stim_type==1),1),2));
V2p=squeeze(nanmean(nanmean(allpow2(PostTimsel,fsel, stim_type==1),1),2));

%% 140Hz BAR pLOT POWER
figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
M=[N1s./N1b, N1s2./N1b];
M2=[V2s./V2b, V2s2./V2b];
%b1=bar(1,nanmean(N1b),'Facecolor',[ 0. 0 0.4]);hold on,
%set(b1,'FaceAlpha',0.7)
b1=bar(2,nanmean(N1s./N1b),'Facecolor',[ 0. 0 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(N1s2./N1b),'Facecolor',[ 0. 0 0.9]);hold on,
set(b1,'FaceAlpha',0.4)
%  b1=bar(4,nanmean(N1p),'Facecolor',[ 0.  0 0.2]);hold on,
%set(b1,'FaceAlpha',0.7)
%b1=bar(6,nanmean(V2b),'Facecolor',[ 0.4 0 0]);hold on,
%set(b1,'FaceAlpha',0.7)
b1=bar(5,nanmean(V2s./V2b),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(6,nanmean(V2s2./V2b),'Facecolor',[ 0.9 0 0]);hold on,
set(b1,'FaceAlpha',0.4)
%  b1=bar(9,nanmean(V2p),'Facecolor',[ 0.2  0 0]);hold on,
%set(b1,'FaceAlpha',0.7)
errorbar([ 2 3  ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
errorbar([  5 6  ],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
ylim([0.8 5])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'BARPLOT_140.pdf'])
savefig(gcf, [ savepath 'BARPLOT_140.fig'])

%%
[~,p,~,stats] =ttest(M(:,1),1) % 40Hz trans
% df=20,  1.8862e-04
[~,p,~,stats] =ttest(M(:,2),1) % 40Hz sustained
% df=20,  2.0015e-04
[~,p,~,stats] =ttest(M(:,1),M(:,2)) % 40Hz sustained
%df=20, 0.9682

[~,p,~,stats] =ttest(M2(:,1),1) % 140Hz trans
% df=25, 1.0157e-05
[~,p,~,stats] =ttest(M2(:,2),1) % 140Hz sustained
% df=25,   0.0020
[~,p,~,stats] =ttest(M2(:,1),M2(:,2)) % 40Hz sustained
%df=25, 5.6779e-05


[~,p,~,stats] =ttest2(M(:,1),M2(:,1)) % difference trans
% df=45,    5.2934e-04
[~,p,~,stats] =ttest2(M(:,2),M2(:,2)) % difference sust
% df=45,    4.8991e-05


fsel=find(freq2.freq>35  & freq2.freq<45); 
V1b=squeeze(nanmean(nanmean(allpow(baseTimsel, fsel,stim_type==0),1),2));
V2b=squeeze(nanmean(nanmean(allpow(baseTimsel, fsel,stim_type==1),1),2));
V1s=squeeze(nanmean(nanmean(allpow(StimTimselTR,fsel, stim_type==0),1),2));
V2s=squeeze(nanmean(nanmean(allpow(StimTimselTR,fsel, stim_type==1),1),2));
V1s2=squeeze(nanmean(nanmean(allpow(StimTimselSU,fsel, stim_type==0),1),2));
V2s2=squeeze(nanmean(nanmean(allpow(StimTimselSU,fsel, stim_type==1),1),2));
V1p=squeeze(nanmean(nanmean(allpow(PostTimsel,fsel, stim_type==0),1),2));
V2p=squeeze(nanmean(nanmean(allpow(PostTimsel,fsel, stim_type==1),1),2));
%% 40Hz BAR pLOT POWER
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
M=[V1b,V1s, V1s2,V1p];M2=[ V2b,V2s, V2s2,V2p];
b1=bar(1,nanmean(V1b),'Facecolor',[ 0 0 0.4]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(2,nanmean(V1s),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V1s2),'Facecolor',[ 0 0 0.9]);hold on,
set(b1,'FaceAlpha',0.4)
b1=bar(4,nanmean(V1p),'Facecolor',[ 0  0 0.2]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2 3 4],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'BARPLOT_40.pdf'])
savefig(gcf, [ savepath 'BARPLOT_40.fig'])

%%%%%%%%%%%%%%%%%%%%%
fsel=find(freq2.freq>135  & freq2.freq<145); 
fsel2=find(freq2.freq>35  & freq2.freq<40); 
V1b=squeeze(nanmean(nanmean(allpow(baseTimsel, fsel2,stim_type==0),1),2));
V2b=squeeze(nanmean(nanmean(allpow(baseTimsel, fsel,stim_type==1),1),2));
V1s=squeeze(nanmean(nanmean(allpow(StimTimselTR,fsel2, stim_type==0),1),2));
V2s=squeeze(nanmean(nanmean(allpow(StimTimselTR,fsel, stim_type==1),1),2));
V1s2=squeeze(nanmean(nanmean(allpow(StimTimselSU,fsel2, stim_type==0),1),2));
V2s2=squeeze(nanmean(nanmean(allpow(StimTimselSU,fsel, stim_type==1),1),2));
V1p=squeeze(nanmean(nanmean(allpow(PostTimsel,fsel2, stim_type==0),1),2));
V2p=squeeze(nanmean(nanmean(allpow(PostTimsel,fsel, stim_type==1),1),2));

%% Base vs trans
[h,p,ci,stats] = ttest(V1b, V1s) % 40Hz
[h,p,ci,stats] = ttest(V2b, V2s)% 140Hz
%% Base vs sustained
[h,p,ci,stats] = ttest(V1b, V1s2) % 40Hz
[h,p,ci,stats] = ttest(V2b, V2s2)% 140Hz
%% Base vs post
[h,p,ci,stats] = ttest(V1b, V1p) % 40Hz
[h,p,ci,stats] = ttest(V2b, V2p)% 140Hz
%% trans between DBS cond
[h,p,ci,stats] = ttest2(V2s./V2b, V1s./V1b) %
%% sust between DBS cond
[h,p,ci,stats] = ttest2(V2s2./V2b, V1s2./V1b) %

%%
Pow_140b=nanmean(nanmean(allpow(baseTimsel,:,stim_type==0),3),1);
Pow_140s=nanmean(nanmean(allpow(StimTimsel,:,stim_type==0),3),1);
Pow_140p=nanmean(nanmean(allpow(PostTimsel,:,stim_type==0),3),1);
Pow_140bS=nanstd(squeeze(nanmean(allpow(StimTimsel,:,stim_type==0),1))./squeeze(nanmean(allpow(baseTimsel,:,stim_type==0),1)),[],2)./sqrt(length(stim_type==0));
Pow_140sS=nanmean(nanstd(allpow(StimTimsel,:,stim_type==1),[],3),1)./sqrt(length(stim_type==1));
Pow_140pS=nanmean(nanstd(allpow(PostTimsel,:,stim_type==1),[],3),1)./sqrt(length(stim_type==1));
figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
plot(freq2.freq,Pow_140s./Pow_140b,'k','COlor',[ 0 0 0.8]);hold on,
%plot(freq2.freq,Pow_140s,'k','COlor',[ 0.8 0 0])
%plot(freq2.freq,Pow_140p,'k','COlor',[ 0.4 0.3 0])
fill_error_area2(freq2.freq,Pow_140s./Pow_140b,Pow_140bS', [ 0.5 0.5 0.5]);
%fill_error_area2(freq2.freq,Pow_140s,Pow_140sS, [ 0.5 0.5 0.5]);
%fill_error_area2(freq2.freq,Pow_140p,Pow_140pS, [ 0.5 0.5 0.5]);
axis tight
%set(gca,'Yscale','log')
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Powplots_40.pdf'])
savefig(gcf, [ savepath 'Powplots_40.fig'])

%%
Pow_140b=nanmean(nanmean(allpow(baseTimsel,:,stim_type==1),3),1);
Pow_140s=nanmean(nanmean(allpow(StimTimsel,:,stim_type==1),3),1);
Pow_140p=nanmean(nanmean(allpow(PostTimsel,:,stim_type==1),3),1);
Pow_140bS=nanstd(squeeze(nanmean(allpow(StimTimsel,:,stim_type==1),1))./squeeze(nanmean(allpow(baseTimsel,:,stim_type==1),1)),[],2)./sqrt(length(stim_type==1));
Pow_140sS=nanmean(nanstd(allpow(StimTimsel,:,stim_type==1),[],3),1)./sqrt(length(stim_type==1));
Pow_140pS=nanmean(nanstd(allpow(PostTimsel,:,stim_type==1),[],3),1)./sqrt(length(stim_type==1));
figure('COlor','w','Position', [ 300 400 250 pheight],'Renderer', 'painters') 
plot(freq2.freq,Pow_140s./Pow_140b,'k','COlor',[ 0.8 0 0]);hold on,
%plot(freq2.freq,Pow_140s,'k','COlor',[ 0.8 0 0])
%plot(freq2.freq,Pow_140p,'k','COlor',[ 0.4 0.3 0])
fill_error_area2(freq2.freq,Pow_140s./Pow_140b,Pow_140bS', [ 0.5 0.5 0.5]);
%fill_error_area2(freq2.freq,Pow_140s,Pow_140sS, [ 0.5 0.5 0.5]);
%fill_error_area2(freq2.freq,Pow_140p,Pow_140pS, [ 0.5 0.5 0.5]);
axis tight
%set(gca,'Yscale','log')
set(gca,'Xscale','log')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Powplots_140.pdf'])
savefig(gcf, [ savepath 'Powplots_140.fig'])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% BARPLOTS POWER


%%% Spike locking %%%%%%%%%%%%%%%%%%%%
if 0
    figure,
    plot(freq2.freq,nanmean(allPLVs(:,stim_type==0),2))
    hold on,plot(freq2.freq,nanmean(allPLVb(:,stim_type==0),2))
    
    figure,
    plot(freq2.freq,nanmean(allPLVs(:,stim_type==1),2))
    hold on,plot(freq2.freq,nanmean(allPLVb(:,stim_type==1),2))
end



V1=allPLVf40(stim_type==0);
V1b=allPLVfSH40(stim_type==0);
V2b=allPLVfSH140(stim_type==1);
V2=allPLVf140(stim_type==1);
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
M=[V1',V1b'];M2=[ V2',V2b'];
b1=bar(1,nanmean(V1),'Facecolor',[ 0 0 0.7]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(2,nanmean(V1b),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(3,nanmean(V2),'Facecolor',[ 0.7 0. 0.]);hold on,
set(b1,'FaceAlpha',0.7)
b1=bar(4,nanmean(V2b),'Facecolor',[ 0.3 0.3 0.3]);hold on,
set(b1,'FaceAlpha',0.7)
errorbar([1 2 ],nanmean(M,1), nanstd(M)./sqrt(size(M,1)),'.k')
errorbar([3 4],nanmean(M2,1), nanstd(M2)./sqrt(size(M2,1)),'.k')
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV_barplot.pdf'])
savefig(gcf, [ savepath 'PLV_barplot.fig'])

[h,p,ci,stats] = ttest(V2, V2b)
% df: 22,  1.2903e-04;  tstat: 4.6316
[h,p,ci,stats] = ttest2(V2, V1)
% df,41,   2.6763e-06, tstat: -5.4425
[h,p,ci,stats] = ttest(V1, V1b)
% df, 19, 1.6595e-08,  tstat: 9.3044
figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
polarhistogram(allangs140(stim_type_sp==1), 20,'Normalization','probability','FaceColor',[ 0.7 0 0] )
rlim([ 0 0.22])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'polarhisto140.pdf'])
savefig(gcf, [ savepath 'polarhisto140.fig'])

figure('COlor','w','Position', [ 300 400 200 pheight],'Renderer', 'painters') 
polarhistogram(allangs40(stim_type_sp==0), 20,'Normalization','probability','FaceColor',[ 0 0 0.7])
rlim([ 0 0.22])
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'polarhisto40.pdf'])
savefig(gcf, [ savepath 'polarhisto40.fig'])
