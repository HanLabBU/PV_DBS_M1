% Housekeeping
clc;
clear all;
close all;

% Setup OS specific path stuff
f = filesep;

%% USER modification
% Pierre's Linux defined rootpath
%server_rootpath = '/home/pierfier/handata_server/';

server_rootpath = 'Z:\';

% Filepath where the matfiles will be saved to
savepath = [server_rootpath 'Pierre Fabris' f 'PV DBS neocortex project' f 'PV_Data' f];

% Flags for the code todo during extractio
save_matfile = 1;
save_videos = 0;
show_figures = 0;

% Add paths to specific functions for new capabilities
if 1
    addpath(genpath([server_rootpath 'EricLowet' f 'Scripts']));
    addpath(genpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'DCIMG_Functions']));
    addpath(genpath([server_rootpath 'Pierre Fabris' f 'Imaging Scripts' f 'sharedLabCode' f 'Camera Read']));
    addpath(genpath([server_rootpath 'Pierre Fabris' f 'DMD Project' 'DMD Scripts' 'In Vitro Analysis Scripts']));
end

%%User modify the save path

% Pierre's data save folder
%savepath = [server_rootpath 'Pierre Fabris' f '' ]

% If data is from old HCImage version (i.e. metadata indeces have changed)
% Switch "code_new" to 0
code_new = 1;

% To navigate quickly to correct folder, switch desired variable to 1 and
% set others to 0.
cd([server_rootpath 'Pierre Fabris' f 'PV DBS neocortex project' f 'Mice voltage recordings' f]);

%DBS_new = 1;
%DBS_old = 0;
%Ultrasound = 0;
%
%if DBS_new
%    cd('C:\DBS_TRANSFER\');
%elseif Ultrasound
%    cd('C:\Ultrasound_Voltage');
%elseif DBS_old
%    cd('D:\DBS')
%else
%    cd('C:\');
%end


[file, path1] = uigetfile('*.dcimg');

cd(path1)
file_split = split(file, '_');
ses = [strjoin({file_split{1:end-1}}, '_'), '_'];
trial = split(file_split{end}, '.');
trialT = trial{1}(end);
filname = [path1 '\' file];

path_split = split(path1, '\');
rec = path_split{end-1}(end);
mouse = path_split{end-2};

trial_ext = split(file_split{end}, '.');
start_trial = str2num(trial_ext{1}(end-1:end));

fov_str = file_split{1}
rec_param = strjoin({file_split{2:3}}, '_')

cell=1;

clear lfpath
lfpath = dir('2_2020-02-07_18*')

stimpath='allstim3_737936.7605_fov2';

stim=dir('allstim3_*');

missing_tr=[];

all_vall=[];trial_vec=[]; all_shifts=[];
trialn=0;
%%

dcimgfile_names = dir([path1 '*.dcimg']);    
isfile_idx = ~[dcimgfile_names.isdir];    
dcimgfile_names = {dcimgfile_names.name};    
dcimgfile_names = dcimgfile_names(isfile_idx);

num_trials = 0;

% Get trial of each matching session defined dcimg file and
% save the last trial
for i=1:length(dcimgfile_names)
    if contains(dcimgfile_names{i}, ses) 
        f_split = split(dcimgfile_names{i}, '_');
        f_split = split(f_split{end}, '.');
        
        if num_trials < str2num(f_split{1})
            num_trials = str2num(f_split{1})
        end
    end
end

for trialT =start_trial : start_trial+num_trials -1
    trialn = trialn+1;
    trialT
    %path1='C:\hipopto\606716\rec10\'
    %path1='C:\hippo_opto\613198\rec3\'
    %path1='C:\hippo_opto\611128\rec3\'
    
    try
        clear all_stim
        for fg=1:length(stim)
            stim_path= [path1 '\' mouse '\' ses '\' stim(fg).name ];
            load(stim_path); all_stim{fg}=allstim;
        end
    catch
    end
    
    %savepath= ['Z:\EricLowet\hippo_opto_main\'  mouse '_rec'  num2str(rec) '_' 'fov' num2str(cell) ''   ];
    
    %savepath= ['Z:\EricLowet\'  mouse '_rec'  num2str(rec) '_' 'fov' num2str(cell) ''   ];
    
    %Z:\EricLowet\hippo_opto_main\RESULTS\DC
    if trialT<10
        %filname= [path1 '\' mouse '\' ses '\' 'fov' num2str(cell) '0000'  num2str(trialT) '.dcimg'];
        filname=  [path1 '\' ses '' '0000'  num2str(trialT) '.dcimg']
    
    elseif trialT>=10  & trialT< 100
        % filname= [path1 '\' mouse '\' ses '\' 'fov' num2str(cell) '000'  num2str(trialT) '.dcimg'];
        filname=  [path1 '\' ses '' '000'  num2str(trialT) '.dcimg'];
    
    elseif trialT>99
        % filname= [path1 '\' mouse '\' ses '\' 'fov' num2str(cell) '000'  num2str(trialT) '.dcimg'];
        filname=  [path1 '\' ses '' '00'  num2str(trialT) '.dcimg']; 
    end

    %filname= [path1 'control_placefielkd_stim00014.dcimg'];
    %ilname= [path1 'try00022.dcimg'];
    
    %  stack = read_dcimg(filname,10);
    %  figure,imagesc(stack(:,:,1))
    
    if code_new
        header_end=2821; % point of data collection -- important! 
        betI=257; %points in between image aquisition
        meta1 = 148;
        meta2 = 147;
        meta3 = 10;
    else
        header_end = 489;
        betI = 17;
        meta1 = 48
        meta2 = 47;
        meta3 = 44;
    end
     
    %489;  % 489 for new point of data collectin %% important! 
    % might change based on camera configuaration
    %% if not correct, check first 1000 valuesof mystack to check when header ends
    
    %   points in between image aqu
    
    fid = fopen(filname,'rb');
    % read the entire stack into memory
    mystack = fread(fid,inf,'uint16=>uint16'); % 16-bit images
    fclose(fid);%
    % read just the first bit of the file (the header plus a bit more)
    fid2 = fopen(filname,'rb');
    metadata = fread(fid2,29450,'uint32=>uint32'); % 32-bit header information
    fclose(fid2);
    
    Height = double(metadata(meta1)); 
    Width = ceil(double(metadata(meta2))/64).*64; % 32=for old voltage scope
    TotalFrames = double(metadata(meta3));
    
    vall=(zeros(Width, Height, TotalFrames,'uint16' )); % allocation
    
    for trial= 1:TotalFrames
        if trial==1
            vv=mystack(header_end:header_end-1+  (Height*Width)  ) ;
            vv(end-4:end)=median(vv);
            vv=reshape(vv(1:end),[Width Height]);
            vv1= [ vv(end-3:end,1:end); vv(1:end-4,1:end)];
            vall(:,:,trial)= vv1;
            lastid= header_end-1+  (Height*Width)  ;
        else
            vv=mystack(lastid+betI:   lastid+(Height*Width)-1+betI );
            vv(end-4:end)=median(vv);
            vv=reshape(vv(1:end),[Width Height]);
            vv1= [ vv(end-3:end,1:end); vv(1:end-4,1:end)];
            vall(:,:,trial)= vv1;
            lastid= lastid+(Height*Width)-1+betI   ;     
        end  
    end
    
    vall(vall==0)=NaN;
    vall=vall(:,:,1:end);

    % Save pre-motion corrected video
    if save_videos && ~exist([path1 'videos'])
        mkdir([path1 'videos']);
    end
    if save_videos && trialT == 1
        save_frames_to_avi(vall, [path1 'videos' f [filname(length(path1)+2:end - 6) '_raw.avi'] ], 1/(1.2*10^-3));
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear XP YP
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if trialn==1 % 
        frame = nanmean(vall(:,:,300:600),3);
        frame(zscore(frame(:))>10)=NaN;
        f1=figure;imagesc(frame);axis image; 
        axis off; colormap(gray);drawnow;
        ROI = imrect;
        
        if isvalid(ROI)
            bw = createMask(ROI);
            pos = round(getPosition(ROI));
            close(f1);
        end
    
        Xwin = pos(2):pos(2)+pos(4); Ywin= pos(1):pos(1)+pos(3);
        roiWindow = vall(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:); 
    else
        roiWindow = vall(pos(2):pos(2)+pos(4), pos(1):pos(1)+pos(3),:); 
        %roiWindow = vall(30:end-30, :,:); % crop a window for motion correction
    end
    
    Y = single(roiWindow);
    SmoothY=(imboxfilt3(Y,[1 1 9]));
    h = fspecial('gaussian',50, 1) - fspecial('gaussian',50, 25);  % high pass the image
    filteredY = gather(imfilter(SmoothY, h, 'replicate', 'same'));
    % run motion estimation
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',30,'max_shift',50,'us_fac',1); % do not use upsampling
    [~,shifts,~] = normcorre(filteredY,options_rigid);
    % apply shift on the image
    %Y = apply_shifts(stack,shifts,options_rigid);
    Y = zeros(size(roiWindow), 'uint16');
    
    for i = 1:size(roiWindow, 3)
        Y(:,:,i) =  circshift(roiWindow(:,:,i), round(shifts(i).shifts));
    end
    
    vall = Y;
    
    
    % Save motion corrected video
    if save_videos && ~exist([path1 'videos'])
        mkdir([path1 'videos']);
    end
    if save_videos && trialT == 1
        save_frames_to_avi(vall, [path1 'videos' f [filname(length(path1)+2:end - 6) '_motion.avi'] ], 1/(1.2*10^-3));
    end

    all_vall=[all_vall ; permute(vall, [ 3 1 2])];
    all_shifts=[all_shifts ;[shifts.shifts]'];
    trial_vec=[trial_vec,ones(1,size(vall,3) ).*trialT];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
all_vall=permute(all_vall,[2 3 1]);

clear allIm,idN=0;
for id=unique(trial_vec)
    idN=idN+1;
    allIm(:,:,idN) = nanmean(all_vall(:,:,trial_vec==id),3);
end

options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',1,'max_shift',50,'us_fac',1); % do not use upsampling
[~,shiftsA,~] = normcorre(allIm,options_rigid);

for i=1:size(all_vall,3)  
    i
    all_vall(:,:,i) =  circshift(all_vall(:,:,i), round(shiftsA(trial_vec(i)-min(trial_vec)+1).shifts));
end

%figure,imagesc(nanmean(all_vall,3))
% 
% 
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%
if 1
    %% select neuron ROI
    averageFrame = mean(vall(:,:,1:end), 3); % take an average frame
    %averageFrame (zscore(averageFrame (:))>5)=NaN;
    cmax = max(averageFrame(:));cmin = min(averageFrame(:));
    averageFrame (isnan(averageFrame))= median(averageFrame(:));
    averageF= smoothn(averageFrame ,1000);
    if 1 
        ROIs= {};
        f1=figure;imagesc(averageFrame(:,:), [cmin, cmax]);%axis image
        axis off;colormap(gray);drawnow;
        mask = zeros(size(averageFrame));
        ROI = drawpolygon;
        
        while isvalid(ROI)
            ROIs{end+1} = createMask(ROI);
            bw = createMask(ROI);
            mask = double(mask | bw);
            imagesc(averageFrame.*(1 - mask.*0.1), [cmin, cmax]);axis image;axis off;colormap(gray);drawnow
            ROI = drawpolygon;
        end
        %close(f1);
    else
        ROIs = {};
        %     % Check if there is an existing ROI.mat file
        %     if isfile([path 'ROIs_' num2str(1) '.mat'])
        %         load([path 'ROIs_' num2str(1)], 'ROIs');
        %         disp("Old ROIs loaded!!!!!!!!!!!!!!!");
        %     end
        %     
        mask = zeros(size(averageFrame));
    
        % Create mask from ROIs array regardless of state
        for i=1:length(ROIs)
            mask = double(mask | ROIs{i});
        end    
            
        figure;imagesc(averageFrame, [cmin, cmax]);
        B = bwboundaries(mask);
        visboundaries(B);%axis image
        axis off;colormap(gray);
        imcontrast(gca);
        zoom on 
        % title('zoom in/contrast')
        pause; drawnow;%title('Select ROIs') %%% Zoom to neuron; select optimal contrast
        ROI = drawpolygon; %  select boundary of neuron ROI

        while isvalid(ROI)
            hold off;
            ROIs{end+1} = createMask(ROI);
            bw = createMask(ROI);
            mask = double(mask | bw);
            imagesc(averageFrame.*(1 - mask.*0.1), [cmin, cmax]);
            axis image;axis off;
       
            for ind=1:length(ROIs)
                [x1]=contourc(double(ROIs{ind}),[0.5, 0.5]);
                hold on,plot(x1(1,2:end),x1(2,2:end),'.r','Markersize',1)
            end;
    
            hold off;
            colormap(gray);
            imcontrast(gca); zoom on;
            title('zoom in/contrast');
            pause;drawnow;
            title('Select ROIs');
            ROI = drawpolygon;pause
        end
        
        close all;    
    end
    
    % %%%%%%%% DENOISING %%%%%%%%%
    % X1=double( vall);r=8;
    % for id=1:size(vall,3)
    %     id
    %      [U,S,V] = svd(X1(:,:,id),'econ'); 
    %    %   X(:,:,ind)= medfilt2(X(:,:,ind),[ 2 4]);       
    %      X1(:,:,id)=  U(:,1:r)*S(1:r,1:r)*V(:,1:r)';
    %       X1(:,:,id)= medfilt2(X1(:,:,id),[ 3 6]);
    % end
    %    vall =uint32(X1); 
    N = length(ROIs);
    traces = zeros(size(all_vall, 3), N);
    for i = 1:N
        neuron = all_vall.*uint16(ROIs{i});
        traces(:,i) = sum(neuron,[1, 2])./sum(ROIs{i},'all');
    end
    
    RR=0;
    
    for id=1:size(traces,2)
        RR=RR+ROIs{id};   
    end
    
    RR=RR>0;
    neuron = all_vall.*uint16(~RR);
    tracesB = squeeze(sum(neuron,[1, 2])./sum(~RR,'all'));
    vmov=all_shifts;vmov=vmov(1:end,:);
    vmovA= fastsmooth(abs(hilbert(sqrt(vmov(:,1).^2  + vmov(:,2).^2))),1,1,1);
    
    if 0
        figure('Color','w')
        subplot(3,1,1)
        for id=1:size(traces,2)
           Vx=( traces(:,id)./tracesB);
           Vx=Vx-fastsmooth(Vx,1000,1,1);
           Vx2=Vx;Vx2(vmovA<30)=NaN;
           plot((Vx)+0.5.*id,'k'); hold on
           plot((Vx2)+0.5.*id,'r'); hold on
           axis tight
        end
        title('Voltage trace')
        subplot(3,1,2)
        
        for id=1:size(traces,2)
           Vx=( traces(:,id)./tracesB);
           Vx=(Vx-fastsmooth(Vx,30,1,1));
           plot((Vx)+0.5.*id,'k'); hold on
           axis tight
        end
        title('Voltage trace')
        subplot(3,1,3)
        plot(vmovA,'k');axis tight
        title('Image movement (combined)')
    end
    
    
    clear REL_SP
    for id=1:size(traces,2)
        Vx=( traces(:,id)./tracesB);;
        Vx4=(Vx-fastsmooth(Vx,40,1,1));
        SP=(zscore(Vx4)>2.7).*0.05;
        SPim=mean(all_vall(:,:,SP>0),3);SPNOim=mean(all_vall(:,:,SP==0),3);
        REL_SP{id}=(SPim - SPNOim)./SPNOim ;REL_SP{id}(REL_SP{id}<0.00)=0;
    end
    
    figure,  imagesc( REL_SP{1} )
    colormap('gray')
    
    clear tracesSP
    for id=1:size(traces,2)
        neuron = all_vall.*uint16(REL_SP{id}.*100).*uint16(ROIs{id});
        tracesSP(:,id) = sum(neuron,[1, 2]);
    end
    
    clear allU
    figure,
    
    for id2=1:1:size(tracesSP,2)
        %figure,plot(zscore(tracesSP./tracesB)+0), hold on,plot(zscore(traces./tracesB))
        Vx4=(traces(1:end,id2)./tracesB(1:end)-fastsmooth(traces(1:end,id2)./tracesB(1:end),40,1,1));
        SP=(zscore(Vx4)>2.35); %Spike threshold, binarized trace
        SP=single(SP);
        SP(1:20)=0;
        SP(end-20:end)=0;
        SP2=SP;
        SP(SP==0)=NaN; 
        subplot(1,round(size(tracesSP,2)./1),id2);
        clear allspM
        
        for id = unique(trial_vec)
            vv =(traces(trial_vec==id,id2)./tracesB(trial_vec==id));
            plot(zscore(vv-fastsmooth(vv,300,1,1))+ 5*id,'k'), hold on,
            plot(SP(trial_vec==id) + 3*id + 5*max(trial_vec)+4 ,'k.','Markersize',6)
            SS=SP2(trial_vec==id);
            SS(1:30)=0;
            SS(end-30:end)=0;
            allspM(1:length(SS),id)= fastsmooth(SS,90,1,1); %calculate spike rate for each trial
            %allU(:,id,id2)= vv;
        end
    
        plot(nanmean(allspM,2).*max(trial_vec)*50 + 3 + 5*max(trial_vec)+4,'r') % plot average spike rate over trials
        axis tight
    end 
        
     
    if 0
        figure('COlor','w')
        subplot(1,1,1)
        
        for id=1:size(traces,2)
            Vx=( tracesSP(:,id)./tracesB);
            Vx=Vx;
            Vx2=Vx;Vx2(vmovA<30)=NaN;    
            X=[ones(length(Vx),1), vmov(:,1)];
            [b,bint,r,rint,stats] = regress(Vx,X);
            Vx=Vx- (vmov(:,1).*b(2)  +    b(1)) ;
            X=[ones(length(Vx),1), vmov(:,2)];
            [b,bint,r,rint,stats] = regress(Vx,X);
            Vx=Vx-vmov(:,2).*b(2) +b(1);Vx3=Vx-fastsmooth(Vx,1500,1,1);
            Vx4=(Vx3-fastsmooth(Vx3,40,1,1));
            SP=(zscore(Vx4)>2.7).*0.05;
            plot((Vx3./1600)+0.5.*(id-1),'k'); hold on
            plot(SP-0.1,'r')
            axis tight
        end
        
        title('Voltage trace')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
    lfp=[]; mm=0;
        
    for kl=1:length(lfpath)  
        lfp_path= [path1 '\' mouse '\' ses '\'  lfpath(kl).name  '\'];
    
        cd( lfp_path)
        %  prenum=124;
        n=dir;
        prenum=(n(10).name(1:3));
        [adc7, timestamps, info] = load_open_ephys_data([ (prenum) '_ADC7.continuous']);  % cMOS time stmaps
        
        if  median(adc7) < 0
            tr_sel=find(adc7<3);
        else
            %tr_sel=find(adc7>0);
            tr_sel=find(adc7>-3);
        end;
        % check polarity of trigger
        tr_ons=[ diff([0;  tr_sel])];
        tr_sel=tr_sel(tr_ons>5);
        tim_cMOS=round( tr_sel./1);
        
        try
            [data, timestamps, info] = load_open_ephys_data([ (prenum) '_Ch1.continuous']);
            data=data(tr_sel); % subsampling to 1000Hz
            [data2, timestamps, info] = load_open_ephys_data([ (prenum) '_Ch2.continuous']);
            data2=data2(tr_sel); % subsampling to 1000Hz
            [data9, timestamps, info] = load_open_ephys_data([ (prenum) '_Ch9.continuous']);
            data9=data9(tr_sel); % subsampling to 1000Hz
        end
     
        [adc5, timestamps, info] = load_open_ephys_data([ (prenum) '_ADC5.continuous']);  %   opto signal if present
        
        opto=adc5(tr_sel); % subsampling to 1000Hz
    
     
        [adc2, timestamps, info] = load_open_ephys_data([ (prenum) '_ADC2.continuous']);  % start CMS trigger
        tr_sel=find(adc2<0);
        tr_ons=[ diff([0;  tr_sel])];
        tr_sel=tr_sel(tr_ons>6);
        start_cMOS=round( tr_sel./1);
     
        %figure,plot(adc2(1:1:end)-4); hold on,plot(adc7(1:1:end))
      
        [adc8, timestamps, info] = load_open_ephys_data([  (prenum) '_ADC8.continuous' ]);  % Motion rec time stmaps
        tr_sel=find(adc8>4);
        tr_ons=[ diff([0;  tr_sel])];
        tr_sel=tr_sel(tr_ons>6);
        tim_MOV=round( tr_sel./1);
    
    
        FS=mean(10000./diff(tim_cMOS));
        
        for trial=1:length(start_cMOS)
            if trial < length(start_cMOS)
                mm=mm+1;
                
                try 
               
                    lfp.trial{mm}(1,:)=  data(tim_cMOS>start_cMOS(trial) & tim_cMOS<start_cMOS(trial+1)  )  ;
                    lfp.trial{mm}(2,:)=  data2(tim_cMOS>start_cMOS(trial) & tim_cMOS<start_cMOS(trial+1)  )  ;
                    lfp.trial{mm}(3,:)=  data9(tim_cMOS>start_cMOS(trial) & tim_cMOS<start_cMOS(trial+1)  )  ;
                end
                
                lfp.opto{mm}(1,:)=  opto(tim_cMOS>start_cMOS(trial) & tim_cMOS<start_cMOS(trial+1)  )  ;
            else
                mm=mm+1;
                
                try
                    lfp.trial{mm}(1,:)=  data(tim_cMOS>start_cMOS(trial) & tim_cMOS<=max(tim_cMOS)  )  
                    lfp.trial{mm}(2,:)=  data2(tim_cMOS>start_cMOS(trial) & tim_cMOS<=max(tim_cMOS)  )  
                    lfp.trial{mm}(3,:)=  data9(tim_cMOS>start_cMOS(trial) & tim_cMOS<=max(tim_cMOS)  )  
                end
                
                lfp.opto{mm}(1,:)=  opto(tim_cMOS>start_cMOS(trial)& tim_cMOS<=max(tim_cMOS) )  ;
            end
            lfp.time{mm}= (1:length(lfp.trial{trial}(1,:)))./FS;
        end
        
        lfp.mov{kl} = tim_MOV;
        lfp.path{kl}= lfpath;
        lfp.tim_cMOS{kl}=tim_cMOS;
        lfp.cmos_start{kl}= start_cMOS;
        lfp.sampling_rate{kl}=FS;
        lfp.label=  {'A'; 'B';'C'};    
    end
     
    if  0
        clear alllfp
        for id=1:length(lfp.trial)
            alllfp(:,id,1)=  abs(hilbert( zscore(lfp.trial{id}(1:4000)-fastsmooth(lfp.trial{id}(1:4000),2600,1,1))));
            alllfp(:,id,2)=   zscore(lfp.opto{id}(1:4000));
        end
            
        figure,plot(nanmean(alllfp(:,:,1),2));
        hold on, plot(nanmean(alllfp(:,:,2),2));
    end
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    result=[];
    result.traces=traces;
    result.traces_weighted=tracesSP;
    result.trial_vec=trial_vec;
    result.ROI= ROIs;
    result.REL_SP=REL_SP;
    result.tracesB=tracesB;
    result.shifts=all_shifts;
    result.dim= [ Height Width];
    
    try
        result.mis=missing_tr;
        result.cond=cond;
        result.lfp=lfp;
        result.allstim=all_stim;
    end
    
    % Save the matfile with the mouse, recording session, FOV number, frequency, and amperage
    if save_matfile
        save([savepath [mouse '_rec' rec '_' fov_str '_' rec_param ]],'result', '-v7.3');
    end
end

% frame = nanmean(vall(:,:,:),3);
% %frame(zscore(frame(:))>10)=0;
% frame=frame./max(frame(:));
%   N = length(ROIs);
% tracesW = zeros(size(vall, 3), N);
% for i = 1:N
%     neuron = vall.*uint16(frame.*ROIs{i});
%     tracesW(:,i) = sum(neuron,[1, 2])./sum(ROIs{i},'all');
% end
% 
% figure,plot(tracesW)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SP=find(zscore(traces-fastsmooth(traces,100,1,1))>3);
% 
% figure,
% imagesc(nanmean(vall(10:end-10,20:end-20,SP),3) -nanmean(vall(10:end-10,20:end-20,:),3)   );

%path2='\\engnas.cbu.edu\research\eng_research_handata\eng_research_handata2\Hua-an Tseng\Data\SomArchon_Striatum\Linear_movement_PinnedBall\Eric_data\'
if  0 
    
    figure
    for ind=1:5:size(vall,3)
        imagesc(all_vall(:,:,ind))
        %imagesc(filteredY(:,:,ind))
        set(gca,'Clim', [300 1700])
        %   set(gca,'Clim', [0 400])
        pause(0.01)  
    end
end

if  1
    
%%%%%%%%%%%%%%%%%%%%%%%% 
clear allF allV

for  ne=unique(result.trial_vec);
    % ne
    %% extract phases
    FS=830;
    %   FS=1000;
    Fn = FS/2;FB=[ 86.4 88.7];
    [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
    LFPg= ((filtfilt(B,A,    result.traces(:,1))));
    result.traces(:,1)= result.traces(:,1)-LFPg;
 
    v= result.traces(result.trial_vec==ne,1)  ;
    v=v-fastsmooth(v,1200,1,1);
    
    allV(1:length(v),ne)=v;
    [ s,w1,t]= spectrogram(zscore(v), 440,438,[1:170],FS);
    
    allF(:,1:size(s,2),ne)= abs(s);end

    S=nanmean(allF,3);
    B=nanmean(nanmean(allF(:,10:50,:),3),2);

    figure('Color','w')
    imagesc(t,w1,smooth2a(nanmean(allF(:,:,1:1:end),3),1,1).*repmat(w1.^0.6,1,length(t)) );axis xy
    colormap(jet)
    
    % set(gca, 'Clim',[-70 300])
  
    figure('Color','w')
    imagesc(t,w1,S./B);axis xy
    colormap(jet)
    ylabel('Frequency')
    xlabel('Time')
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Color','w');clear allV allVS
    T=[];
    
    for  ne=unique(result.trial_vec)
        v= result.traces(result.trial_vec==ne); %eb edited 20211111
        %   v= result.traces(result.trial_vec==ne,3);%./result.tracesB(result.trial_vec==ne)   ;
        v=v-fastsmooth(v,1400,1,1); 
        vS=double(zscore(v-fastsmooth(v,7,1,1))>2.3);
        vS(vS==0)=NaN;
        vS(1:50)=NaN;allV(:,ne)=v;
        vS(end-50:end)=NaN;
        v=zscore(v); 
        v(1:50)=NaN;;
        v(end-50:end)=NaN;;
        allVS(:,ne)=vS;
        plot((1:length(v))./FS,((v-nanmean(v))./nanstd(v))./6 + ne+1,'k'); 
        T=[T;(v)];
        hold on,
        plot((1:length(v))./FS,vS+ ne*1+0.6,'.r','Markersize',11)
    end
   
    axis tight 
    allVS(isnan(allVS))=0;
    
    figure,plot(nanmean(allV(:,1:1:end),2))
    hold on,plot(nanmean(allV(:,1:1:end),2))
    plot(nanmean(allV(:,1:1:end),2))
    hold on,plot(nanmean(allV(:,1:1:end),2))
       
    figure('COlor','w'),
    plot(((1:length(T))./828)./60,T,'k')   
    axis tight
    
    %   [ s,w,t]= spectrogram(zscore(mean(allV,2)), 400,350,[1:70],830)
    %   
    %   figure
    %   imagesc(t,w,abs(s));axis xy
    %   colormap(jet)
      
      
    % figure('COlor','w'),plot((1:length(v))./828,fastsmooth(nanmean(allVS(:,1:4:end),2),50,1,1));hold on,
    % plot((1:length(v))./828,fastsmooth(nanmean(allVS(:,2:4:end),2),50,1,1))
    % plot((1:length(v))./828,fastsmooth(nanmean(allVS(:,3:4:end),2),50,1,1))
    % plot((1:length(v))./828,fastsmooth(nanmean(allVS(:,4:4:end),2),50,1,1))
    % axis tight
    
end
