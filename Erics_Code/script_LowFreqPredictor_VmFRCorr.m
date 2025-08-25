clear all
f = filesep;

%% Addpaths
%addpath(genpath('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\'))
%cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\PV_Data\')

% Linux computer paths
server_root_path = '~/handata_server/';
local_root_path = '~/Projects/';
pv_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];

savepath = Multi_func.save_plot;

addpath([local_root_path f 'Pierre Fabris' f 'PV DBS neocortex' f 'Scripts' f]);


%% Reset variables
switch_Source=2; % 1=DBS, 2=Vm

for id=1:2; for id2=1:3; allPLVb2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVs2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVp2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allVmDBS{id,id2}=[];end;end

allFIRING=[];allDEP=[];allCH=[];allSA=[];aCOND=[];   allFIR_ratio=[];
ses = dir([pv_path '*.mat']);
all_matfiles = {ses.name};
          
[region_matfiles] = Multi_func.find_region(all_matfiles);        
all_matfiles = {ses.name};
gg= fieldnames(region_matfiles)' ;   
region=0; mr=0;

%% Just look at M1 region
for f_region =gg(2)% fieldnames(region_matfiles)'
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
    for dbs =1:2 % 1 is 40Hz, 2 is 140Hz
        field = DBS_cond{dbs};  
        cond=cond+1;
        matfiles = matfile_stim.(field).names;    
     
        % Initialize field to keep track of neurons and trials
        cond_stats.(field) = struct();
        cond_stats.(field).trial_nums = [];
        cond_stats.(field).num_fovs = 0;
    
        % Loop through each matfile of the current stimulation condition
       
        fnames=matfiles;
        for iter=1:length(fnames)
            matfile = matfiles(iter)
            frs=[2:1:10 11:1:150 ];
            FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([pv_path matfile{1}]);
            clear alls wavD  allVm
            SNR=[];SA=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
                if ~isempty(  data.align.trial{tr}) 
                    if mean(data.align.trial{tr}.spike_info375.spike_snr{1})>0 % Trial included if spike SNR above 4SNR
                        tr_iter= tr_iter+1;
                        voltsig=data.align.trial{tr}.spike_info375.trace_ws;
                        [fitbaseline, coeff]=Multi_func.exp_fit_Fx(voltsig',FS); %remove photobleaching
                        voltsig=(voltsig-fitbaseline);
                        
                        if  1
                            Camtime=  data.raw.trial{tr}.raw_camera_frame_time- data.raw.trial{tr}.raw_camera_start_time(1);
                            DBStim=data.raw.trial{tr}.raw_stimulation_time- data.raw.trial{tr}.raw_camera_start_time(1);
                            
                            if length(Camtime) <100
                                load([ server_root_path 'EricLowet'  f 'Camtime.mat'])
                            end

                            traceDBS=zeros(1,length(data.raw.trial{tr}.raw_traces));
                            for id=1:length(DBStim)
                                z=find(Camtime > DBStim(id));
                                try
                                    traceDBS(z(1)-4)=1;
                                end
                            end
                        end

                        if switch_Source==1
                            voltsig=traceDBS;
                        end
                        allVm(:, tr_iter)=voltsig(1:2496);
                        [filt_sig]=Multi_func.filt_data( allVm(:, tr_iter),frs, FS); % filtering+Hilbert
                        wt=filt_sig';
                         %   [wt, f,f1] = cwt(voltsig, round(FS),'TimeBandwidth',100);
                          %  frs=f;
                    
                        % Hilbert transform of filtered data
                        wavD(:,:, tr_iter)=wt';
                        SNR = [SNR;(data.align.trial{tr}.spike_info375.spike_snr{1})];
                        % Spike Amplitude
                        SA=[SA;data.align.trial{tr}.spike_info375.spike_amplitude{1}'];
                        alls(:, tr_iter)=[ data.align.trial{tr}.spike_info375.roaster];
                   end 
               end
            end
       
            if tr_iter   >0 
                baseTimsel= [200:FS-100];
                StimTimsel= [FS+80:FS*2-100];
                PostTimsel= [FS*2+10:size(alls,1)];
   
                if  mean(SNR) >4  % mean spike SNR for the neuron to be over 4SNR (can be changed)
                    % Read notes, account for spike interference in Vm power value
                    shift_time=0; % shift from spike backwards
                    clear s M
                    Mbase=nanmean(nanmean(abs(wavD(baseTimsel,:,:)),3),1);
                
                    clear s M;
                    Mstim=nanmean(nanmean(abs(wavD(StimTimsel,:,:)),3),1);
                   
                    mr=mr+1

                    % Normalized Spectra Power Change
                    allP(:,mr)=(Mstim-Mbase)./(Mstim+Mbase);
                           
                    % Average firing rate for the baseline
                    spB= nanmean(nanmean(alls(baseTimsel,:),1)).*FS  ;
                    % Average firing rate for stimulation
                    spST= nanmean(nanmean(alls(StimTimsel,:),1)).*FS  ;
                    voltsig=    fastsmooth(nanmean(allVm,2),1,1,1)';
                    [ fitbaseline, coeff]=Multi_func.exp_fit_Fx(voltsig',FS); %remove photobleaching
                    %   voltsig=(voltsig-fitbaseline);
                
                    allDEP=[allDEP,  (voltsig./nanmean(SA))'];
                    aCOND=[aCOND,dbs];
                    % Firing rate change from stimulation to baseline
                    allFIR_ratio= [allFIR_ratio, spST-spB];
                    % allVmDBS{region,cond}(1:length(allVm),mr) = nanmean(allVm,2);
                    % 
                    % allPLVb{region,cond}(1:length(PLV_b),mr)=PLV_b;
                    % allPLVs{region,cond}(1:length(PLV_s),mr)=PLV_s;
                    % allPLVp{region,cond}(1:length(PLV_p),mr)=PLV_p;
                    % 
                    % allPLVb2{region,cond}=[allPLVb2{region,cond},  PLV_b2];
                    % allPLVs2{region,cond}=[allPLVs2{region,cond},  PLV_s2];
                    % allPLVp2{region,cond}=[allPLVp2{region,cond},  PLV_p2];
                end  
            end
        end
    end
end

%% Final calculation from consolidated data
pheight=400;

freq2.freq=frs;
Vm= zscore(nanmean(allDEP(StimTimsel,:),1));

% Average Normalized VM during stimulation time
VmA= (nanmean(allDEP(StimTimsel,:),1));

% Selecting frequency range within 2-10
fsel=freq2.freq>2 & freq2.freq<10 ;%| freq2.freq>6 & freq2.freq<10

LF=zscore(nanmean(allP(fsel,:),1));

% Averaging across the lower frequencies, 2-10Hz
LFA=nanmean(allP(fsel,:),1);

fsel=freq2.freq>138 & freq2.freq<142;
fsel2=freq2.freq>38 & freq2.freq<42
SF=([zscore(nanmean(allP(fsel2,aCOND==1),1))  zscore(nanmean(allP(fsel,aCOND==2),1))  ]);

% Normlized Stimulation Frequency Power Change
SFA=([(nanmean(allP(fsel2,aCOND==1),1))  (nanmean(allP(fsel,aCOND==2),1))  ]);

M_Vm_40=mean(VmA(aCOND==1)) ;
M_Vm_140=mean(VmA(aCOND==2)) ;

M_Fr_40med=median(allFIR_ratio(aCOND==1));
M_Fr_140med=median(allFIR_ratio(aCOND==2))  ;
M_Fr_40=mean(allFIR_ratio(aCOND==1));
M_Fr_140=mean(allFIR_ratio(aCOND==2));

M_SF_40=mean(SFA(aCOND==1)) ;
M_SF_140=mean(SFA(aCOND==2)) ;

%  save('\\engnas.bu.edu\research\eng_research_handata\EricLowet\DBS\M1_quants','M_Vm_40','M_Vm_140','M_Fr_40med','M_Fr_140med','M_Fr_40','M_Fr_140','M_SF_40','M_SF_140')

%% Plot of entrainment to low frequency power
figure('COlor','w','Position', [ 300 400 600 pheight],'Renderer', 'painters'); 
plot(SFA(aCOND==2),(LFA(aCOND==2)),'.','Markersize',10, 'COlor',[ 0 0.4 0.5], 'DisplayName', '140Hz');
hold on;
plot(SFA(aCOND==1),LFA(aCOND==1),'.b','Markersize',10, 'COlor',[ 0.6 0.4  0], 'DisplayName', '40Hz');
axis tight;
xdata1=SFA(aCOND==1);
ydata1=LFA(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0.6 0.4  0],'Linewidth',1.5);
xdata1=SFA(aCOND==2);
ydata1=LFA(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0.4 0.5],'Linewidth',1.5);
legend();
title('Entrainment(Power) vs Low Frequency Power Change', 'Interpreter', 'none');
%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Raster_plot_SF' '.pdf'])
saveas(gcf, [savepath 'Correlation' f 'Entrainment(Power) vs Low Frequency Power.png' ]);

%% Plotting Vm during stimulation vs low frequency power
figure('COlor','w','Position', [ 300 400 600 pheight],'Renderer', 'painters');
plot(VmA(aCOND==2),(LFA(aCOND==2)),'.','Markersize',10, 'COlor',[ 0.5 0 0], 'DisplayName', '140Hz');
hold on;
plot(VmA(aCOND==1),LFA(aCOND==1),'.r','Markersize',10, 'COlor',[ 0 0.4  0.5], 'DisplayName', '40Hz');
axis tight;
xdata1=VmA(aCOND==2);
ydata1=LFA(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.02:max(xdata1));

% Plot the fit
fitLine1 = plot(min(xdata1):0.02:max(xdata1),yplot1,'DisplayName','140Hz linear','XLimInclude','off',...
    'Tag','linear',...
    'MarkerSize',6,...
    'Color',[0.5 0 0]);
xdata1=VmA(aCOND==1);
ydata1=LFA(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.02:max(xdata1));

% Plot the fit
fitLine1 = plot(min(xdata1):0.02:max(xdata1),yplot1,'DisplayName','40Hz linear','XLimInclude','off',...
    'Tag','linear',...
    'MarkerSize',6,...
    'COlor',[0 0.4  0.5]);
legend();
title('Vm during stim vs. 2-10Hz Frequency Power Change');
saveas(gcf, [savepath 'Correlation' f 'Vm during stim vs. low frequency power change.png']);

%print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Raster_plot_Vm' '.pdf'])

%% Calculate the variance explained of the multiple regressions of Vm and power at stim frequency
X = [ones(size(Vm')) Vm' SF' ];
%[b,bint,r,rint,stats]

% LF is the zscore of the power of low frequency Vm
%X is 
b= regress(LF(aCOND==1)',X(aCOND==1,:)).^2 
b2= regress(LF(aCOND==2)',X(aCOND==2,:)).^2 

figure('COlor','w','Position', [ 300 400 600 pheight],'Renderer', 'painters') 
bar([ 1 ],[b(2) ].*100,0.8,'FaceColor', [0.5 0 0], 'DisplayName', 'Vm 40Hz'); 
hold on;
bar([ 4],[b2(2) ].*100,0.8,'FaceColor', [0.99 0 0], 'DisplayName', 'Vm 140Hz'); 
hold on;
bar([ 2 ],[b(3) ].*100,0.8,'FaceColor', [0 0 0.5 ], 'DisplayName', 'ENT 40Hz');
bar([ 5 ],[b2(3) ].*100,0.8,'FaceColor', [0 0 0.99 ], 'DisplayName', 'ENT 140Hz');
xlim([0 6])
title('Variance explained of LF power vs VM and ENT');
legend();
saveas(gcf, [savepath 'Correlation' f 'LF_vs_Vm_and_Ent_variance' '.png']);
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Correlation' f 'LF_vs_Vm_and_Ent_variance' '.pdf'])


%% Vma-LFA, Vm during stim vs. Normalized Lower frequency Power Change
figure('COlor','w','Position', [ 300 400 600 pheight],'Renderer', 'painters') 
% 40Hz
subplot(1,2,1)
plot(VmA(aCOND==1),-LFA(aCOND==1),'.b','Markersize',12, 'COlor',[ 0.6 0.4  0]);
hold on;
axis tight;
xdata1=VmA(aCOND==1);ydata1=-LFA(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0  0],'Linewidth',1.5);
xlabel('Vm Change');
ylabel('Power Change');
xlim([-0.7 1.5]);
ylim([-.6 1]);
line([ 0 0 ], [-.6 1],'Color','k');
title('40Hz', 'Interpreter', 'none');

% 140Hz
subplot(1,2,2);
plot(VmA(aCOND==2),(-LFA(aCOND==2)),'.','Markersize',12, 'COlor',[ 0 0.4 0.5]);
hold on;
axis tight;
xdata1=VmA(aCOND==2);
ydata1=-LFA(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0  0 0],'Linewidth',1.5);
xlim([-0.7 1.5]); ylim([-.6 1]); line([ 0 0 ], [-.6 1],'Color','k') ;
title('140Hz', 'Interpreter', 'none');
sgtitle('Vm during stimulation vs. Low Frequency Vm Power change (neg)');
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'Correlation' f 'M1_Vm_LF' '.pdf'])


%% SFA-LFA
figure('COlor','w','Position', [ 300 400 350 pheight],'Renderer', 'painters') 
subplot(1,2,1)
plot(SFA(aCOND==1),-LFA(aCOND==1),'.b','Markersize',12, 'COlor',[ 0.6 0.4  0]);
hold on;
axis tight;
xdata1=SFA(aCOND==1);ydata1=-LFA(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0 0],'Linewidth',1.5);
xlim([-0.2 0.7]);  ylim([-.6 1]); line([ 0 0 ], [-.6 1],'Color','k') ;
title('40Hz');

subplot(1,2,2)
plot(SFA(aCOND==2),(-LFA(aCOND==2)),'.','Markersize',12, 'COlor',[ 0 0.4 0.5])
hold on ,
axis tight;
xdata1=SFA(aCOND==2);ydata1=-LFA(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0 0],'Linewidth',1.5);
xlim([-0.2 0.7]); ylim([-.6 1]); line([ 0 0 ], [-.6 1],'Color','k') ;
title('140Hz');
sgtitle('38-42Hz Vm Power Change vs 2-10Hz Vm Power Change (neg)');
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'M1_SF_LF' '.pdf'])

%% Vma-SFA
figure('COlor','w','Position', [ 300 400 350 pheight],'Renderer', 'painters') 
% 40Hz
subplot(1,2,1)
plot(VmA(aCOND==1),SFA(aCOND==1),'.b','Markersize',12, 'COlor',[ 0.6 0.4  0]);
hold on;
axis tight;
xdata1=VmA(aCOND==1);ydata1=SFA(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0  0 0 ],'Linewidth',1.5);
xlim([-0.7 1.5]); ylim([-.3 1]);line([ 0 0 ], [-0.3 1],'Color','k');

% 140Hz
subplot(1,2,2);
plot(VmA(aCOND==2),(SFA(aCOND==2)),'.','Markersize',12, 'COlor',[ 0 0.4 0.5]);
hold on;
axis tight;
xdata1=VmA(aCOND==2);ydata1=SFA(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0 0],'Linewidth',1.5);
xlim([-0.7 1.5]); ylim([-.3 1]);line([ 0 0 ], [-0.3 1],'Color','k');

print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'M1_Vm_SFA' '.pdf']);



%% Vm change to Firing Rate
figure('COlor','w','Position', [ 300 400 800 pheight],'Renderer', 'painters');
% 40Hz
subplot(1,2,1);
plot(VmA(aCOND==1),allFIR_ratio(aCOND==1),'.b','Markersize',12, 'COlor',[ 0.6 0.4  0]);
hold on;
axis tight;
xdata1=VmA(aCOND==1);
ydata1=allFIR_ratio(aCOND==1);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0.0 0.0  0],'Linewidth',1.5);
xlim([-0.7 1.5]); ylim([-2 8]); line([ 0 0 ], [-2 8],'Color','k');
title('40Hz');

% 140Hz
subplot(1,2,2);
plot(VmA(aCOND==2),(allFIR_ratio(aCOND==2)),'.','Markersize',12, 'COlor',[ 0 0.4 0.5]);
hold on;
axis tight;
xdata1=VmA(aCOND==2);
ydata1=allFIR_ratio(aCOND==2);
fitResults1 = polyfit(xdata1,ydata1,1);
yplot1 = polyval(fitResults1,min(xdata1):0.01:max(xdata1));
fitLine1 = plot(min(xdata1):0.01:max(xdata1),yplot1,'DisplayName','linear','XLimInclude','off', 'Tag','linear','MarkerSize',6,'Color',[ 0 0.0 0.0],'Linewidth',1.5);
xlim([-0.7 1.5]); ylim([-2 8]); line([ 0 0 ], [-2 8],'Color','k');
title('140Hz');
sgtitle('Vm during stimulation vs. Firing rate change');

saveas(gcf, [savepath 'PLV' f 'M1_Vm_FIR' '.png']);
print(gcf, '-dpdf' , '-r300' ,'-painters', [ savepath 'PLV' f 'M1_Vm_FIR' '.pdf']);


%% Testing to do other stuff, but no longer
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
    


function [ fitbaseline, coeff]=exp_fit_Fx(v,FS)
    % fits an expoential to estimate the photobleaching
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




function [ PLV_output,PLV_output2]= spike_field_ppc_Pierre(wavD,spikes,timshift)
    if size(wavD,4)>1
        [tr chs frs t]= size(wavD);
    else
        [t frs tr]= size(wavD);
    end
    mat=[];
    for ind=1:tr
        % Question: when is the size of wavD different?
        if size(wavD,4)>1
            M=exp(1i.*((squeeze(wavD(ind,CH,:,:))))); 
        else
            M=exp(1i.*((squeeze(wavD(:,:,ind))))); end
            s=find(spikes{ind});
            s=s-timshift; s(s<=0)=[];
            mat=[mat,    M(s,:)'];
        end
    
        NT= sum(~isnan(mat(1,:)));
        Z=abs(nanmean(mat,2));  % MVL
        T=   Z.^2; 
        PLV_output= (((1/(NT-1))*((T.*NT-1))));  %adjusted MLV (PPC)
        PLV_output2= mat; 
    if NT<10
        PLV_output=PLV_output.*NaN;
        %PLV_output2=PLV_output2.*NaN;
    end
end
