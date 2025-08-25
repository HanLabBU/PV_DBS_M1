clear all
f = filesep;

%addpath(genpath('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\'))
%cd('\\engnas.bu.edu\research\eng_research_handata\eng_research_handata3\Pierre Fabris\PV Project\PV_Data\')

% Linux computer paths
server_root_path = '~/handata_server/';
local_root_path = '~/Projects/';
pv_path = [server_root_path 'eng_research_handata3' f 'Pierre Fabris' f 'PV Project' f 'PV_Data' f];
addpath([local_root_path f 'Pierre Fabris' f 'PV DBS neocortex' f 'Scripts' f]);

savepath = Multi_func.save_plot;

switch_Source=2; % 1=DBS, 2=Vm

%% Loop through each region
allAREA=[];

for id=1:2; for id2=1:3; allPLVb2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVs2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allPLVp2{id,id2}=[];end;end
for id=1:2; for id2=1:3; allVmDBS{id,id2}=[];end;end


allFIRING=[];
allDEP=[];
allCH=[];
allSA=[];
allPOW = [];
ses = dir([pv_path '*.mat']);
all_matfiles = {ses.name}
          
[region_matfiles] = Multi_func.find_region(all_matfiles);        
all_matfiles = {ses.name};
     
region=0; %area=1--> V1, area=2--> M1
fnams= fieldnames(region_matfiles)'
for f_regiong =2
    region=region+1;
    f_region = fnams{f_regiong};% f_region{1};

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
    for dbs =2
        field = DBS_cond{dbs};  
        cond=cond+1;
        matfiles = matfile_stim.(field).names;    
     
        % Initialize field to keep track of neurons and trials
        cond_stats.(field) = struct();
        cond_stats.(field).trial_nums = [];
        cond_stats.(field).num_fovs = 0;
    
        % Loop through each matfile of the current stimulation condition
        mr=0;
        mats=matfiles;
        for matfile =[1:length(mats)]
            frs=[2:1:10 11:1:170 ];FS=round(1000./1.2);
            % Read in the mat file of the current condition
            data = load([pv_path mats{matfile}]);
            clear alls wavD  allVm
            SNR=[];SA=[];
            tr_iter=0;
            for tr=1:length(data.align.trial)  % trials
                if ~isempty(  data.align.trial{tr}) 
                    if mean(data.align.trial{tr}.spike_info2.spike_snr{1})>0 % Trial included if spike SNR above 4SNR
                        tr_iter= tr_iter+1;
                        voltsig=data.align.trial{tr}.spike_info2.trace_ws;
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
                        SNR = [SNR;(data.align.trial{tr}.spike_info2.spike_snr{1})];
                        SA=[SA;data.align.trial{tr}.spike_info2.spike_amplitude{1}'];
                        alls(:, tr_iter)=data.align.trial{tr}.spike_info2.roaster;
                   end 
               end
            end
       
            if tr_iter   >8
                baseTimsel= [10:FS-10];
                StimTimsel= [FS:FS*2];
                PostTimsel= [FS*2+10:size(alls,1)];
   
                if  mean(SNR) >4  % mean spike SNR for the neuron to be over 4SNR (can be changed)
                    % Read notes, account for spike interference in Vm power value
                    shift_time=0; % shift from spike backwards
                  mr=mr+1
                   
                  % Firing Rate is binned at 5 points
                  FR=fastsmooth(nanmean(alls,2),5,1,1);
                 allFIRING=[allFIRING,    FR-nanmean(FR(100:FS-10))];
            %      allDEP=[allDEP,  nanmean(allVm,2)];
              voltsig=    fastsmooth(nanmedian(allVm,2),1,1,1)';
                          [ fitbaseline, coeff]=Multi_func.exp_fit_Fx(voltsig',FS); %remove photobleaching
           %   voltsig=(voltsig-fitbaseline);

                
                % DEBUG
                %pow_size = size(allPOW)
                %dep_size = size(allDEP)

                allPOW(:,:,mr)= nanmean(abs(wavD),3);
                allDEP=[allDEP,  (voltsig./nanmean(SA))'];
                
                allCH=[allCH, matfile];
                allSA=[allSA,nanmean(SA)];
                allAREA=[allAREA, f_regiong];
               
                end 
                  
            end
        end
    end
end
    
pheight=430;
%figure('COlor','w','Position', [ 300 400 150 120],'Renderer', 'painters')
% figure('COlor','w','Position', [ 300 400 380 150],'Renderer', 'painters')
% fill_error_area2((1:size(  allFIRING,1))./827-1,nanmean((  allFIRING(:,allAREA==1)),2),nanstd((  allFIRING(:,allAREA==1)),[],2)./sqrt(size(allDEP,2)),[.5 .5 .5])
% plot((1:size(  allFIRING,1))./827-1,nanmean((  allFIRING(:,allAREA==1)),2),'b','Linewidth',1); hold on,
% fill_error_area2((1:size(  allFIRING,1))./827-1,nanmean((  allFIRING(:,allAREA==2)),2),nanstd((  allFIRING(:,allAREA==2)),[],2)./sqrt(size(allDEP,2)),[.5 .5 .5])
% plot((1:size(  allFIRING,1))./827-1,nanmean((  allFIRING(:,allAREA==2)),2),'r','Linewidth',1); hold on,
% 
% axis tight
% xlim([-0.1 0.2])
% 
% 

% This is sorting the median of the Vm during the first half of stimulation
[A B]=sort(median(allDEP(1000:1500,:),1));


%%  Plotting all of the Vm depolarizations
col=colormap(jet(size(allDEP,2)));
figure('COlor','w','Position', [ 300 400 840 450],'Renderer', 'painters');
%fill_error_area2((1:size(allDEP,1))./828-1,nanmedian((allDEP(:,allAREA==2)),2),nanstd((allDEP(:,allAREA==2)),[],2)./sqrt(size(allDEP,2)),[.5 .5 .5])
for id=1:size(allDEP,2)
    plot((1:size(allDEP,1))./828-1,fastsmooth((allDEP(:,B(id))),100,1,1),'k','Linewidth',1,'Color', col((id),:)); 
    hold on;
end
hold on;
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);
axis tight;
xlim([-0.5 1.8]);
title('Plotting the smoothed Vm of each neurons');


%% Plotting the Vms with color indicating whether hyperpolarized, depolarized, or neutral
col=colormap(jet(size(allDEP,2)));
figure('COlor','w','Position', [ 300 400 840 450],'Renderer', 'painters');
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);
hold on;
%fill_error_area2((1:size(allDEP,1))./828-1,nanmedian((allDEP(:,allAREA==2)),2),nanstd((allDEP(:,allAREA==2)),[],2)./sqrt(size(allDEP,2)),[.5 .5 .5])
for id=1:size(allDEP,2)
    if A(id)< -0.2
        plot((1:size(allDEP,1))./828-1,fastsmooth((allDEP(:,B(id))),150,1,1),'b','Linewidth',1); 
        hold on;
    elseif A(id) > 0.2
        plot((1:size(allDEP,1))./828-1,fastsmooth((allDEP(:,B(id))),150,1,1),'r','Linewidth',1); 
        hold on;
    else
        plot((1:size(allDEP,1))./828-1,fastsmooth((allDEP(:,B(id))),150,1,1),'k','Linewidth',1); 
        hold on;
    end
end
axis tight;
xlim([-0.5 1.8]);
title('Plotting Vm with color indicating depolarization, hyperpolarization, and neutral');
saveas(gcf, [savepath 'VmChange' f 'Summary_Vm_change_category.png']);
saveas(gcf, [savepath 'VmChange' f 'Summary_Vm_change_category.pdf']);
% figure('COlor','w','Position', [ 300 400 280 470],'Renderer', 'painters')
% subplot(2,1,1)
% imagesc((1:size(allDEP(100:2300,:),1))./828-1,[],zscore(smooth2a(nanmedian(allPOW(100:2300,:,B(A>0.4)),3),200,1),[],1)')
% axis xy
% colormap(jet);%set(gca,'CLim', [ 0 4])
% subplot(2,1,2)
% imagesc((1:size(allDEP(100:2300,:),1))./828-1,[],zscore(smooth2a(nanmedian(allPOW(100:2300,:,B(A<-0.4)),3),200,1),[],1)')
% axis xy
% colormap(jet);%set(gca,'CLim', [ 0 4])

%% Plotting the lower frequency power across time from each Vm change category
% Find the median of the Vm Power across most of the trace period avereaged across 2-9Hz frequency
%V1 is the hyperpolarized Vm
%V2 is the depolarized
%V3 is neutral Vm
V1=nanmedian(nanmean(allPOW(300:2300,2:10,B(A<-0.2)),2),3);
V2=nanmedian(nanmean(allPOW(300:2300,2:10,B(A>0.2)),2),3);
V3=nanmedian(nanmean(allPOW(300:2300,2:10,B(A<0.2  & A>-0.2)),2),3);
figure('COlor','w','Position', [ 300 400 840 450],'Renderer', 'painters');

subplot(1,3,1);
plot((1:size(allDEP(300:2300,:),1))./828-0.64,V1- mean(V1(1:200)),'b','Linewidth',1.5);
axis tight;
hold on;
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);
xlim([-0.5 1.8]);
ylim([-2 2]);
title('Hyperpolarized'); 
ylabel(['Low Freq pow']);

subplot(1,3,2);
plot((1:size(allDEP(300:2300,:),1))./828-0.64,V3- mean(V3(1:200)),'k','Linewidth',1.5);
axis tight;
hold on;
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);xlim([-0.5 1.8]);ylim([-4 4]);
title('No/weak change');
ylabel(['Low Freq pow']);


subplot(1,3,3);
plot((1:size(allDEP(300:2300,:),1))./828-0.64,V2- mean(V2(1:200)),'r','Linewidth',1.5);
axis tight;
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);xlim([-0.5 1.8]);
hold on;
ylim([-4 4]);
title('Depolarized'); 
ylabel(['Low Freq pow']);

saveas(gcf, [savepath 'VmChange' f 'Low_Freq_by_category.pdf']);
saveas(gcf, [savepath 'VmChange' f 'Low_Freq_by_category.png']);

% figure,plot(nanmedian(nanmean(allPOW(100:2300,100:end,B(A<-0.2)),2),3));hold on,plot(nanmedian(nanmean(allPOW(100:2300,100:end,B(A>0.2)),2),3))

%% Plotting Firing rate change across time for each category
% Firing rate is baseline subtracted
V1=fastsmooth(nanmedian(nanmean(allFIRING(300:2300,B(A<-0.2)),2),3),100,1,1).*828;
V2=fastsmooth(nanmedian(nanmean(allFIRING(300:2300,B(A>0.2)),2),3),100,1,1).*828;
V3=fastsmooth(nanmedian(nanmean(allFIRING(300:2300,B(A< 0.2 & A>-0.2)),2),3),100,1,1).*828;
figure('COlor','w','Position', [ 300 400 840 450],'Renderer', 'painters')
subplot(1,3,1)
,plot((1:size(allDEP(300:2300,:),1))./828-0.6,V1- mean(V1(1:400)),'b','Linewidth',1.5);axis tight
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);xlim([-0.5 1.8]);%ylim([-2 2])
title('Hyperpolarized'); ylabel(['Firing Hz change']);ylim([-3 3])
subplot(1,3,2)
plot((1:size(allDEP(300:2300,:),1))./828-0.6,V3- mean(V3(1:400)),'k','Linewidth',1.5)
axis tight
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);xlim([-0.5 1.8]);ylabel(['Firing Hz change']);ylim([-3 3])
title('No/weak change');
subplot(1,3,3)
plot((1:size(allDEP(300:2300,:),1))./828-0.6,V2- mean(V2(1:400)),'r','Linewidth',1.5)
axis tight
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1);xlim([-0.5 1.8])
%ylim([-4 4])
title('Depolarized'); ylabel(['Firing Hz change']);ylim([-3 3])

saveas(gcf, [savepath 'VmChange' f 'Firing_Rate_by_category.pdf']);
saveas(gcf, [savepath 'VmChange' f 'Firing_Rate_by_category.png']);

%% Plotting todo
figure('COlor','w','Position', [ 300 400 380 150],'Renderer', 'painters')
%fill_error_area2((1:size(allDEP,1))./828-1,nanmedian((allDEP(:,allAREA==2)),2),nanstd((allDEP(:,allAREA==2)),[],2)./sqrt(size(allDEP,2)),[.5 .5 .5])
for id=1:size(allDEP,2)
plot((1:size(allDEP,1))./828-1,fastsmooth((allDEP(:,B(id))),100,1,1),'k','Linewidth',1,'Color', col((id),:)); hold on,
end
line([-0.5 1.8 ], [ 0 0],'COlor', [ 0 0 0],'Linewidth',1)
axis tight
xlim([-0.5 1.8])


% 
% figure,plot(nanmean(allDEP(850:1800,:),1),log(nanmean(allFIRING(850:1800,:),1)),'ok')
% 
% V1=nanmean(allFIRING(850:1800,:),1)-nanmean(allFIRING(100:800,:),1);
% [a1 b1]=sort(V1)
% figure
% subplot(1,2,1)
% imagesc(allFIRING(:,b1)');axis xy
% colormap(jet)
% subplot(1,2,2)
% imagesc(allDEP(1:end-100,b1)');axis xy
% colormap(jet)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

          
            
% function [ fitbaseline, coeff]=exp_fit_Fx(v,FS)
%     %% fits an expoential to estimate the photobleaching
%     v1=v;
%     v1=v1([250:FS-10]);%= mean([ v([FS-30:FS-10 2000:2020 ])]);
%     v1(1)=v1(2);
%     F = @(x,xdata)x(1)+x(2)*(- xdata./x(3));%+ x(3)*exp(- xdata./x(4))  ;
%     x0 = [mean(v1) 40 1.5   ] ;
%     OPTIONS = optimoptions('lsqcurvefit','Algorithm','levenberg-marquardt');
%     t = (1:length(v1))./FS; t2 = (1:length(v))./FS;
%     tsel=1:length(t);
%     [xunc,RESNORM,RESIDUAL] = lsqcurvefit(F, x0, t(tsel)', v1(tsel),[],[], OPTIONS);
%     fitbaseline=xunc(1)+xunc(2)*(-t2./xunc(3));
%     coeff=xunc;
% % figure,plot(t,fit_baseline)
% %       hold on,plot(t,v)
% end



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
