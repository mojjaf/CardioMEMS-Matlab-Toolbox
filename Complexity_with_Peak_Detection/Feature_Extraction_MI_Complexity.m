function [ features,features_avg,features_fusion_median,features_acc,features_gyro,labID_ext] = Feature_Extraction_MI_Complexity(data,fs,win,step,lab, id,segment)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(segment, 'systole')
    start=1;
    stop=500; %400%250%50
elseif strcmp(segment, 'diastole')
    start=500;%400%250%50
    stop=950;%800%650%150
elseif strcmp(segment, 'fullcycle')
    start=1;
%     stop=1000; %strech on
    if fs<400 %sterch off
         stop=300;
    else
         stop=500;
    end
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

fs=200;

[avpeak,ploc_ECG] = averagedpeaks(detrend(-data.ekg(curPos:curPos+windowLength-1))',fs); % ORIG
% ploc_ECG=ploc_ECG(2:end)-0;

    aX  = data.accX(curPos:curPos+windowLength-1);
    aY  = data.accY(curPos:curPos+windowLength-1);
    aZ  = data.accZ(curPos:curPos+windowLength-1);
    gX  = data.gyroX(curPos:curPos+windowLength-1);
    gY  = data.gyroY(curPos:curPos+windowLength-1);
    gZ  = data.gyroZ(curPos:curPos+windowLength-1);
   
% figure;
% plot(detrend(-data.ekg(curPos:curPos+windowLength-1))'); hold on;
% plot(ploc_ECG,ploc_ECG.*0+2,'*','Color','r','MarkerSize',14);
% figure;
ploc_ECG(1)=[];
shifter=20;
%%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
[alltraces_ax,avpeak_accx,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aX),fs, ploc_ECG-shifter); 
[alltraces_ay,avpeak_accy,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aY),fs, ploc_ECG-shifter); 
[alltraces_az,avpeak_accz,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZ),fs, ploc_ECG-shifter); 
[alltraces_gx,avpeak_gyrx,FS_gyrx] = averagedpeaks_ECG_REF_mj(detrend(gX),fs, ploc_ECG-shifter); 
[alltraces_gy,avpeak_gyry,FS_gyry] = averagedpeaks_ECG_REF_mj(detrend(gY),fs, ploc_ECG-shifter); 
[alltraces_gz,avpeak_gyrz,FS_gyrz] = averagedpeaks_ECG_REF_mj(detrend(gZ),fs, ploc_ECG-shifter); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  FS=1000;
 Features_Gx(i,:) = MI_FeatureExtraction_FrameWise(alltraces_ax(:,start:stop),FS_accx);
 Features_Gy(i,:) = MI_FeatureExtraction_FrameWise(alltraces_ay(:,start:stop),FS_accy);
 Features_Gz(i,:) = MI_FeatureExtraction_FrameWise(alltraces_az(:,start:stop),FS_accz);
 Features_Ax(i,:)= MI_FeatureExtraction_FrameWise(alltraces_gx(:,start:stop),FS_gyrx);
 Features_Ay(i,:)= MI_FeatureExtraction_FrameWise(alltraces_gy(:,start:stop),FS_gyry);
 Features_Az(i,:) = MI_FeatureExtraction_FrameWise(alltraces_gz(:,start:stop),FS_gyrz);
%     [coeff_tot,scores_tot,latent,~,PTvar_tot]=pca([acc;gyro]);

%     Mechanical_Complexity_Ratios(3,i)=latent(2)/latent(1);

    curPos = curPos + step;
    
end


if lab
labels=ones(size(Features_Ax,1),1);
% labels_x=ones(6,1);
labels_x=ones(size(Features_Ax,1),1);
 else
   labels=zeros(size(Features_Ax,1),1);
   labels_x=zeros(size(Features_Ax,1),1);
 end
idx=zeros(size(labels));
idx(:,:)=id;

ID_X=zeros(size(labels_x));
ID_X(:,:)=id;
features=[Features_Gx,Features_Gy,Features_Gz,Features_Ax,Features_Ay,Features_Az,labels,idx];
Features_3D(1,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Gx;
Features_3D(2,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Gy;
Features_3D(3,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Gz;
Features_3D(4,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Ax;
Features_3D(5,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Ay;
Features_3D(6,1:size(Features_Gx,1),1:size(Features_Gx,2))=Features_Az;
features_avg=median(Features_3D,1);
% features_avg=[F_gx,F_gy,F_gz,F_ax,F_ay,F_az,labels_x,ID_X];
features_avg = reshape(features_avg,size(features_avg,1)*size(features_avg,2),size(features_avg,3));


features_fusion_median=[median(Features_Gx),median(Features_Gy),median(Features_Gz),median(Features_Ax),median(Features_Ay),median(Features_Az)];

features_acc=[Features_Ax,Features_Ay,Features_Az];
features_gyro=[Features_Gx,Features_Gy,Features_Gz];
labID_ext=[labels,idx];
% TMC=mean(Mechanical_Complexity_Ratios(3,:));;
% TMC_sd=std(Mechanical_Complexity_Ratios(3,:));;
end

