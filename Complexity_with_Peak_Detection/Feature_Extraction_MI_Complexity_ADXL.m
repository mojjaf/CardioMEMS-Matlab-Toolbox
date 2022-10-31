function [ features,features_avg,features_fusion_median,features_acc,features_gyro,labID_ext] = Feature_Extraction_MI_Complexity_ADXL(data,fs,win,step,lab, id,segment)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(segment, 'systole')
    start=1;
    stop=400;
elseif strcmp(segment, 'diastole')
    start=400;
    stop=800;
elseif strcmp(segment, 'fullcycle')
    start=1;
    stop=1000;%1000
end
    
    
signal=data.accX;
windowLength = round(win * fs);
step = round(step * fs);

curPos = 1;
L = length(signal);

% compute the total number of frames:
numOfFrames = floor((L-windowLength)/step) + 1;
for i=1:numOfFrames % for each frame
    % get current frame:
     ploc=0; avpeak=0; ploc_ECG=0;

fs=400;

[avpeak,ploc_ECG] = averagedpeaks(detrend(-data.ekg(curPos:curPos+windowLength-1))',fs); % ORIG
% ploc_ECG=ploc_ECG(2:end)-0;

    aX  = data.accX(curPos:curPos+windowLength-1);
    aY  = data.accY(curPos:curPos+windowLength-1);
    aZ  = data.accZ(curPos:curPos+windowLength-1);
    gX  = data.accXLF(curPos:curPos+windowLength-1);
    gY  = data.accYLF(curPos:curPos+windowLength-1);
    gZ  = data.accZLF(curPos:curPos+windowLength-1);
   
% figure;
% plot(detrend(-data.ekg(curPos:curPos+windowLength-1))'); hold on;
% plot(ploc_ECG,ploc_ECG.*0+2,'*','Color','r','MarkerSize',14);
% figure;
ploc_ECG(1)=[];

shifter=50;
%%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
[alltraces_axHF,avpeak_accx,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aX),fs, ploc_ECG-shifter); 
[alltraces_ayHF,avpeak_accy,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aY),fs, ploc_ECG-shifter); 
[alltraces_azHF,avpeak_accz,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZ),fs, ploc_ECG-shifter); 
[alltraces_axLF,avpeak_gyrx,FS_accxLF] = averagedpeaks_ECG_REF_mj(detrend(gX),fs, ploc_ECG-shifter); 
[alltraces_ayLF,avpeak_gyry,FS_accyLF] = averagedpeaks_ECG_REF_mj(detrend(gY ),fs, ploc_ECG-shifter); 
[alltraces_azLF,avpeak_gyrz,FS_acczLF] = averagedpeaks_ECG_REF_mj(detrend(gZ),fs, ploc_ECG-shifter); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FS=1000;
 Features_axHF(i,:) = MI_FeatureExtraction_FrameWise(alltraces_axHF(:,start:stop),FS_accx);
 Features_ayHF(i,:) = MI_FeatureExtraction_FrameWise(alltraces_ayHF(:,start:stop),FS_accy);
 Features_azHF(i,:) = MI_FeatureExtraction_FrameWise(alltraces_azHF(:,start:stop),FS_accz);
 Features_axLF(i,:)= MI_FeatureExtraction_FrameWise(alltraces_axLF(:,start:stop),FS_accxLF);
 Features_ayLF(i,:)= MI_FeatureExtraction_FrameWise(alltraces_ayLF(:,start:stop),FS_accyLF);
 Features_azLF(i,:) = MI_FeatureExtraction_FrameWise(alltraces_azLF(:,start:stop),FS_acczLF);
%     [coeff_tot,scores_tot,latent,~,PTvar_tot]=pca([acc;gyro]);

%     Mechanical_Complexity_Ratios(3,i)=latent(2)/latent(1);

    curPos = curPos + step;
    
end


if lab
labels=ones(size(Features_axLF,1),1);
% labels_x=ones(6,1);
labels_x=ones(size(Features_axLF,1),1);
 else
   labels=zeros(size(Features_axLF,1),1);
   labels_x=zeros(size(Features_axLF,1),1);
 end
idx=zeros(size(labels));
idx(:,:)=id;

ID_X=zeros(size(labels_x));
ID_X(:,:)=id;
features=[Features_axHF,Features_ayHF,Features_azHF,Features_axLF,Features_ayLF,Features_azLF,labels,idx];
Features_3D(1,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_axHF;
Features_3D(2,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_ayHF;
Features_3D(3,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_azHF;
Features_3D(4,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_axLF;
Features_3D(5,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_ayLF;
Features_3D(6,1:size(Features_axHF,1),1:size(Features_axHF,2))=Features_azLF;
features_avg=median(Features_3D,1);
% features_avg=[F_gx,F_gy,F_gz,F_ax,F_ay,F_az,labels_x,ID_X];
features_avg = reshape(features_avg,size(features_avg,1)*size(features_avg,2),size(features_avg,3));


features_fusion_median=[median(Features_axHF),median(Features_ayHF),median(Features_azHF),median(Features_axLF),median(Features_ayLF),median(Features_azLF)];

features_acc=[Features_axLF,Features_ayLF,Features_azLF];
features_gyro=[Features_axHF,Features_ayHF,Features_azHF];
labID_ext=[labels,idx];
% TMC=mean(Mechanical_Complexity_Ratios(3,:));;
% TMC_sd=std(Mechanical_Complexity_Ratios(3,:));;
end

