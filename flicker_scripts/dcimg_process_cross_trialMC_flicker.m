%% Script for Extracting Voltage Traces from DCIMGs
% adapted from Eric Lowet by Emma Bortz
% Detects spikes using spike_detect_SNR_v3b
% Need to change read_dcimg function as appropriate for your camera
clear; clc;
close all;
f = filesep;
%% Read file
% addpath(genpath('\\engnas.bu.edu\research\eng_research_handata\EricLowet\'))
addpath(genpath(['~/handata_server' f 'Emma_Bortz' f 'scripts' f]));
% addpath(genpath('Z:\EricLowet\'))
addpath(genpath('./'));

% source folder for dcimg files
datapath = ['~/handata_server/eng_research_handata3/Yangyang_Wang/PV_V1_LED_SomArchon/TODO/169337_V1/'];

[fname,fdir] = uigetfile([datapath '*.dcimg'],'MultiSelect','on');

if ~iscell(fname)
    len = 1;
    savename = fname(1:end-6);
else
    len = length(fname);
    savename = [fname{1}(1:end-11),'_all'];
end

% savename = 'FOV2_all'

all_vall = []; trial_vec = []; all_shifts = [];
roi_list = struct;

% Sampling frequency for dual scope camera
Fs = 1/(2*10^-3);
for T = 1:len % motion correction over each trial
    
    if iscell(fname)
        file = fname{T};
        roi_list.file(T).name = file;
    else
        file = fname;
        roi_list.file.name = file;
    end
    
    cur_filename = fname{T};
    cur_filename = cur_filename(1:strfind(cur_filename, '.dcimg') - 1)

    path = fdir;
    fprintf(['Processing Trial ' num2str(T),'\n'])
%     stack = read_dcimg_dualscope_cam2([path, file]);% calcium camera
%     stack = read_dcimg_dualscope_cam1([path, file]); % new voltage camera
%     Height = size(stack,1);
%     Width = size(stack,2);

%   for somArchon
%    [stack,Height,Width,TotalFrames] = read_dcimg_kk_scope_new([path, file]); % Originally 
 
    % For electra
    stack = read_dcimg_dualscope_cam2([path, file]); % Originally read_dcimg_kk_scope_new
     Height = size(stack,1);
     Width = size(stack,2);
    
    %     [stack,Height,Width,TotalFrames] = read_dcimg_kk_scope([path, file]);

%     stack = read_dcimg_voltageroom([path, file]); % for Rebecca room camera

    % pick an ROI for motion correction
    if T == 1 % 
        frame = mean(stack(:,:,300:600),3,'omitnan');
        frame(zscore(frame(:))>10) = NaN; % any pixels too high variance get NaN
        
        f1 = figure;imagesc(frame);axis image; 
        axis off; colormap(gray); drawnow;
        ROI = imrect;

        if isvalid(ROI)
            bw = createMask(ROI);
            pos = round(getPosition(ROI));
            close(f1);
        end

        Xwin = pos(2):pos(2)+pos(4); 
        Ywin = pos(1):pos(1)+pos(3);

        roiWindow = stack(Xwin, Ywin,:); 
    else
        roiWindow = stack(Xwin, Ywin,:); 
    end
    
    Y = single(roiWindow);

    % Y = stack; % without imrect for motion correction

    % smooth the image stack with 3D box filtering
    try
        SmoothY = gpuArray(imboxfilt3(Y,[1 1 11])); % use GPU
    catch
        SmoothY = imboxfilt3(Y,[1 1 11]); % don't use GPU
    end
    h = fspecial('gaussian',50, 1) - fspecial('gaussian',50, 25);  % high pass the image
    filteredY = gather(imfilter(SmoothY, h, 'replicate', 'same'));

    % run motion estimation
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',30,'max_shift',50,'us_fac',1); % do not use upsampling
    [~,shifts,~] = normcorre(filteredY,options_rigid);

    Y = zeros(size(stack), 'uint16');
    for i = 1:size(stack, 3)
        Y(:,:,i) =  circshift(stack(:,:,i), shifts(i).shifts);
    end
    stack = Y; % apply motion shifts to stack

    % Save motion corrected video
    if not(isfolder([path 'motion_corrected' f]))
        mkdir([path 'motion_corrected' f]);
    end

    save_frames_to_avi(Y, [path 'motion_corrected' f 'motcorrected_' cur_filename], Fs);
    disp('Saved motion corrected video');
    

    all_vall = [all_vall ; permute(stack, [ 3 1 2])];
    all_shifts = [all_shifts ;[shifts.shifts]'];
    trial_vec = [trial_vec,ones(1,size(stack,3) ).*T];

    roi_list.file(T).name = file;
    if iscell(fname)
        roi_list.file(T).shifts = shifts;
    else
        roi_list.file.shifts = shifts;
    end
    beg_inds(T) = (T - 1)*size(stack, 3) + 1; % trial beg. indices
    end_inds(T) = (T)*size(stack, 3); % trial end indices
end

all_vall = permute(all_vall, [2 3 1]);

% motion correction over trials
clear allIm,idN=0;
for id = unique(trial_vec)
    idN = idN + 1;
    allIm(:,:,idN) = mean(all_vall(:, :, trial_vec==id), 3,'omitnan');
end
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',1,'max_shift',50,'us_fac',1); % do not use upsampling
[~,shiftsA,~] = normcorre(allIm, options_rigid);

for i=1:size(all_vall,3)  
    all_vall(:,:,i) =  circshift(all_vall(:,:,i), round(shiftsA(trial_vec(i)-min(trial_vec)+1).shifts));
end
%% select neuron ROI
averageFrame = mean(stack(:,:,1:end),3); % take an average frame
cmax = max(averageFrame(:))*1;
cmin = min(averageFrame(:))*1;

averageFrame(isnan(averageFrame)) = median(averageFrame(:));

ROIs= {};
figure;imagesc(averageFrame, [cmin, cmax]);axis image;axis off;colormap(gray);drawnow;
imcontrast;
mask = zeros(size(averageFrame));
ROI = drawpolygon;
while isvalid(ROI)
    ROIs{end+1} = createMask(ROI);
    bw = createMask(ROI);
    mask = double(mask | bw);
    imcontrast;
    imagesc(averageFrame.*(1 - mask.*0.1), [cmin, cmax]);axis image;axis off;colormap(gray);drawnow
    ROI = drawpolygon;
end % will stop when you close the window
close all;
%%
upthresh = 4; 
downthresh = 3;

nskip = 20;
t_starts = find(diff([0 trial_vec])); %starts of trials
t_starts = cumsum([t_starts; ones(nskip-1, size(t_starts,2))]);
t_starts = reshape(t_starts,[1,numel(t_starts)]);

t_ends = find(diff([trial_vec 0])); % ends of trials
t_ends = cumsum([(t_ends-nskip-1); ones(nskip-1, size(t_ends,2))]);
t_ends = reshape(t_ends,[1,numel(t_ends)]);

N = length(ROIs);
traces = zeros(size(all_vall, 3), N);
for i = 1:N
    neuron = all_vall.*uint16(ROIs{i});
    traces(:,i) = sum(neuron,[1, 2])./sum(ROIs{i},'all');
end

RR = 0; % Create mask with all ROIs
for id = 1:size(traces, 2)
    RR = RR + ROIs{id};   
end
RR = RR > 0; 

neuron = all_vall.*uint16(~RR); % create background mask
tracesB = squeeze(sum(neuron,[1, 2])./sum(~RR,'all')); % get background trace

vmov = all_shifts; vmov = vmov(1:end,:);
vmovA = fastsmooth(abs(hilbert(sqrt(vmov(:,1).^2  + vmov(:,2).^2))),1,1,1);
% stack = reshape(Y,[],size(Y,3));

clear REL_SP
figure;
Vx_all = (traces./repmat(tracesB,1,size(traces,2))); % normalize to background
result = spike_detect_SNR_v3b(Vx_all,upthresh,downthresh);
for rr = 1:size(traces,2) % loop through ROIs    
    SP = result.roaster(rr,:);
   
    SP = single(SP); SP(t_starts) = 0; SP(t_ends) = 0; SP2 = SP; % spikes
    
    Vx4 = (Vx_all(:,rr)-fastsmooth(Vx_all(:,rr), 40, 1, 1)); % fastsmooth
    
    SPim = mean(all_vall(:, :, SP > 0),3); % Average spike Image?
    SPNOim = mean(all_vall(:, :, SP == 0),3); % average NO spike Image?

    REL_SP{rr} = (SPim - SPNOim)./SPNOim ; % normalized spike image
    REL_SP{rr}(REL_SP{rr} < 0.00) = 0;
    SP(SP==0) = NaN; 
    subplot(1,round(size(Vx_all,2)./1), rr); clear allspM
    for tr = unique(trial_vec)
        vv = (traces(trial_vec == tr, rr)./tracesB(trial_vec == tr)); 
        plot(zscore(vv - fastsmooth(vv, 300, 1, 1)) + 5*tr, 'k'); hold on;
        plot(SP(trial_vec == tr) + 3*tr + 5*max(trial_vec) + 4 , 'k.', 'Markersize', 6)
        SS = SP2(trial_vec == tr); SS(1:30) = 0; SS(end-30:end) = 0;
        allspM(1:length(SS), tr) = fastsmooth(SS, 90, 1, 1);
    end
    % plot PSTH
    plot(mean(allspM, 2, 'omitnan').*max(trial_vec)*50 + 3 + 5*max(trial_vec) + 4, 'r')
    axis tight

end

figure; % This is plotting neuron hot spots during spikes
imagesc( REL_SP{1} )
colormap('gray')
axis image

clear tracesSP
for id = 1:size(traces,2)
    neuron = all_vall.*uint16(REL_SP{id}.*100).*uint16(ROIs{id});
    tracesSP(:,id) = sum(neuron,[1, 2]);
end

clear allU
figure;

for id2 = 1:size(tracesSP,2)
    Vx4(:,id2) = (tracesSP(1:end,id2)./tracesB(1:end)-fastsmooth(tracesSP(1:end,id2)./tracesB(1:end),40,1,1));
end
result2 = spike_detect_SNR_v3b(Vx4,upthresh,downthresh);

for id2 = 1:size(tracesSP,2) % loop through ROIs
    SP = result2.roaster(id2,:);
    SP = single(SP); SP(t_starts) = 0; SP(t_ends) = 0; SP2 = SP; % spikes
    SP(SP==0) = NaN; 
    subplot(1,round(size(tracesSP,2)./1), id2); clear allspM
    for id = unique(trial_vec)
        vv = (traces(trial_vec == id, id2)./tracesB(trial_vec == id)); 
        plot(zscore(vv - fastsmooth(vv, 300, 1, 1)) + 5*id, 'k'); hold on;
        plot(SP(trial_vec == id) + 3*id + 5*max(trial_vec) + 4 , 'k.', 'Markersize', 6)
        SS = SP2(trial_vec == id); SS(1:30) = 0; SS(end-30:end) = 0;
        allspM(1:length(SS), id) = fastsmooth(SS, 90, 1, 1);
        %   allU(:,id,id2)= vv;
    end
    % plot PSTH
    plot(mean(allspM, 2, 'omitnan').*max(trial_vec)*50 + 3 + 5*max(trial_vec) + 4, 'r')
    axis tight
 end 
%% create struct to save data: 
for iii = 1:len
    roi_list.file(iii).trace = traces(beg_inds(iii):end_inds(iii));
end
roi_list.traces = traces;
roi_list.traces_weighted = tracesSP;
roi_list.trial_vec = trial_vec;
roi_list.ROI = ROIs;
roi_list.REL_SP = REL_SP;
roi_list.tracesB = tracesB;
roi_list.shifts = all_shifts;
roi_list.dim = [ Height Width];
roi_list.result = result;
roi_list.result_weighted = result2;

cd(path)
save(['RawTracesV2_' savename '.mat'],'roi_list')

%% plots
if  1
    clear allF allV
    for  ne = unique(roi_list.trial_vec)
        % extract phases
        FS=500;
%           FS=1000/1.5;
        Fn = FS/2;FB=[ 86.4 88.7];
        [B, A] = butter(2, [min(FB)/Fn max(FB)/Fn]);
        LFPg= ((filtfilt(B,A, roi_list.traces(:,1))));
        %   roi_list.traces(:,1)= roi_list.traces(:,1)-LFPg;
        v= roi_list.traces(roi_list.trial_vec==ne)  ;
        v=v-fastsmooth(v,1200,1,1);

        allV(1:length(v),ne)=v;
        [ s,w1,t]= spectrogram(zscore(v), 340,338,[1:170],FS);

        allF(:,1:size(s,2),ne)= abs(s);
    end

    S = mean(allF,3,'omitnan');
    B = mean(mean(allF(:,10:50,:),3,'omitnan'),2,'omitnan');
% 
%     figure('Color','w')
%     imagesc(t,w1,smooth2a(nanmedian(allF(:,:,1:1:end),3),1,1).*repmat(w1.^0.6,1,length(t)) );axis xy
%     colormvap(jet)
%     ylabel('Frequency')
%     xlabel('Time (s)')
%     % set(gca, 'Clim',v[-70 300])

    figure('Color','w')
    imagesc(t,w1,S./B);axis xy
    colormap(jet)
    ylabel('Frequency')
    xlabel('Time (s)')

%    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure('Color','w');clear allV allVS allspM
    T=[];
%     subplot(2,1,2)
goodtrials = 1:10;
    for ne = unique(roi_list.trial_vec) % loop through trials
%     for ne = goodtrials 
        v = roi_list.traces(roi_list.trial_vec == ne, 1)./roi_list.tracesB(roi_list.trial_vec == ne);
        v = v-fastsmooth(v,1400,1,1); 
%         vS = double(zscore(v-fastsmooth(v,7,1,1))>2.6); % get spikes
        vS = roi_list.result.roaster(1, roi_list.trial_vec == ne);
        vS(vS==0) = NaN; vS(1:nskip) = NaN; 
        allV(:,ne) = v; vS(end-nskip:end) = NaN; 
        v = zscore(v); v(1:nskip) = NaN; v(end-nskip:end) = NaN;
        allVS(:, ne) = vS;        
        allVS(isnan(allVS)) = 0;
        allspM(1:length(vS), ne) = fastsmooth(allVS(:, ne), 90, 1, 1);
        plot((1:length(v))./FS,((v-mean(v,'omitnan'))./nanstd(v))./6 + ne+1,'k'); 
        T = [T; (v)];
        hold on,
        plot((1:length(v))./FS, vS + ne*1 + 0.6, '.r','Markersize',6)
    end
%     for id = unique(trial_vec)
%         vv = (traces(trial_vec == id, id2)./tracesB(trial_vec == id)); 
%         plot(zscore(vv - fastsmooth(vv, 300, 1, 1)) + 5*id, 'k'); hold on;
%         plot(SP(trial_vec == id) + 3*id + 5*max(trial_vec) + 4 , 'k.', 'Markersize', 6)
%         SS = SP2(trial_vec == id); SS(1:30) = 0; SS(end-30:end) = 0;
%         allspM(1:length(SS), id) = fastsmooth(SS, 90, 1, 1);
%     end
%     plot PSTH
%     plot(mean(allspM, 2, 'omitnan').*max(trial_vec)*50 + 3 + 5*max(trial_vec) + 4, 'r')
 
%     subplot(2,1,1) % PSTH
%     plot((1:length(v))./FS, vS + ne*1 + 0.6, '.r','Markersize',11)
    hold on;
    plot((1:length(v))./FS,mean(allspM, 2, 'omitnan').*max(trial_vec)*10 + ne*1 + 2, 'r')
 
    axis tight  
    allVS(isnan(allVS)) = 0;
    xlabel('Time (s)')
 
%     figure; plot(mean(allV(:,1:1:end),2,'omitnan'))
%     hold on; plot(mean(allV(:,1:1:end),2,'omitnan'))
%     plot(mean(allV(:,1:1:end),2,'omitnan'))
%     hold on,plot(mean(allV(:,1:1:end),2,'omitnan'))

%     figure('Color','w'),
%     plot(((1:length(T))./828)./60,T,'k')   
%     axis tight
%     xlabel('Time (s)')
end
