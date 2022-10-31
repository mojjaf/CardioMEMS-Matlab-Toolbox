function [features_out] = angle_features(data, fs, win, step,segment)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
if strcmp(segment, 'systole')
    start=1;
    stop=450;
elseif strcmp(segment, 'diastole')
    start=350;
    stop=750;
elseif strcmp(segment, 'fullcycle')
    start=1;
    stop=300;%1000
end
    

signal=data.accX;
windowLength = round(win * fs);
step = round(step * fs);

curPos = 1;
L = length(signal);

% compute the total number of frames:
numOfFrames = floor((L-windowLength)/step) + 1;
SCG_cmplx = zeros(1, numOfFrames);
SCG_mob = zeros(1, numOfFrames);
GCG_cmplx = zeros(1, numOfFrames);
GCG_mob = zeros(1, numOfFrames);

for i=1:numOfFrames % for each frame
    % get current frame:
    aX  = data.accX(curPos:curPos+windowLength-1);
    aY  = data.accY(curPos:curPos+windowLength-1);
    aZ  = data.accZ(curPos:curPos+windowLength-1);
    gX  = data.gyroX(curPos:curPos+windowLength-1);
    gY  = data.gyroY(curPos:curPos+windowLength-1);
    gZ  = data.gyroZ(curPos:curPos+windowLength-1);


    ploc=0; avpeak=0; ploc_ECG=0;

fs=200;

[avpeak,ploc_ECG] = averagedpeaks(detrend(-data.ekg(curPos:curPos+windowLength-1))',fs); % ORIG
% ploc_ECG=ploc_ECG(2:end)-20;
% figure;
% plot(detrend(-data.ekg(curPos:curPos+windowLength-1))'); hold on;
% plot(ploc_ECG,ploc_ECG.*0+2,'*','Color','r','MarkerSize',14);
% figure;
[alltraces_ax,avpeak_accx,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aX),fs, ploc_ECG); 
[alltraces_ay,avpeak_accy,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aY),fs, ploc_ECG); 
[alltraces_az,avpeak_accz,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZ),fs, ploc_ECG); 
[alltraces_gx,avpeak_gyrx,FS_gyrx] = averagedpeaks_ECG_REF_mj(detrend(gX),fs, ploc_ECG); 
[alltraces_gy,avpeak_gyry,FS_gyry] = averagedpeaks_ECG_REF_mj(detrend(gY ),fs, ploc_ECG); 
[alltraces_gz,avpeak_gyrz,FS_gyrz] = averagedpeaks_ECG_REF_mj(detrend(gZ),fs, ploc_ECG); 

%%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles

 alpha_acc=(acos(alltraces_ax(:,start:stop)./sqrt(alltraces_ax(:,start:stop).^2+alltraces_ay(:,start:stop).^2+alltraces_az(:,start:stop).^2)));
 beta_acc=(acos(alltraces_ay(:,start:stop)./sqrt(alltraces_ax(:,start:stop).^2+alltraces_ay(:,start:stop).^2+alltraces_az(:,start:stop).^2)));
 gamma_acc=(acos(alltraces_az(:,start:stop)./sqrt(alltraces_ax(:,start:stop).^2+alltraces_ay(:,start:stop).^2+alltraces_az(:,start:stop).^2)));
 alpha_gyr=( acos(alltraces_gx(:,start:stop)./sqrt(alltraces_gx(:,start:stop).^2+alltraces_gy(:,start:stop).^2+alltraces_gz(:,start:stop).^2)));
 beta_gyr= (acos(alltraces_gy(:,start:stop)./sqrt(alltraces_gx(:,start:stop).^2+alltraces_gy(:,start:stop).^2+alltraces_gz(:,start:stop).^2)));
 gamma_gyr= (acos(alltraces_gz(:,start:stop)./sqrt(alltraces_gx(:,start:stop).^2+alltraces_gy(:,start:stop).^2+alltraces_gz(:,start:stop).^2)));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

angle_aX(i)=median(tpr(alpha_acc'));
angle_aY(i)=median(tpr(beta_acc'));
angle_aZ(i)=median(tpr(gamma_acc'));
angle_gX(i)=median(tpr(alpha_gyr'));
angle_gY(i)=median(tpr(beta_gyr'));
angle_gZ(i)=median(tpr(gamma_gyr'));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hjorth based complexity and mobility

%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curPos = curPos + step;
    
end





features_out=[angle_aX',angle_aY',angle_aZ',angle_gX',angle_gY',angle_gZ'];
end





