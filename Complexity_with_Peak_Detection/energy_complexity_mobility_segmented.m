function [features_out] = energy_complexity_mobility_segmented(data, fs, win, step,segment)
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
    
SCG=sqrt(aX.^2+aY.^2+aZ.^2);

    ploc=0; avpeak=0; ploc_ECG=0;

fs=200;

[avpeak,ploc_ECG] = averagedpeaks(detrend(-data.ekg(curPos:curPos+windowLength-1))',fs); % ORIG
% ploc_ECG=ploc_ECG(2:end)-20;
% figure;
% plot(detrend(-data.ekg(curPos:curPos+windowLength-1))'); hold on;
% plot(ploc_ECG,ploc_ECG.*0+2,'*','Color','r','MarkerSize',14);
% figure;

%%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
[all_SCG_traces,avpeak_accx,ploc_accx] = averagedpeaks_ECG_REF(detrend(SCG),fs,3, ploc_ECG); 

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hjorth Based Complixty and Mobility
   [total_SCG_comlpx,total_SCG_mob] = complxity_and_mobility_segmented(all_SCG_traces(:,start:stop));
 
   
    SCG_cmplx(i)=median(abs(total_SCG_comlpx));
    SCG_mob(i)=median(abs(total_SCG_mob));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    gX  = data.gyroX(curPos:curPos+windowLength-1);
    gY  = data.gyroY(curPos:curPos+windowLength-1);
    gZ  = data.gyroZ(curPos:curPos+windowLength-1);
%       gyro=[gX;gY;gZ];
      GCG=sqrt(gX.^2+gY.^2+gZ.^2);
      [all_GCG_traces,avpeak_gyrx,ploc_gyrx] = averagedpeaks_ECG_REF(detrend(GCG),fs,3, ploc_ECG); 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Hjorth based complexity and mobility
[total_GCG_comlpx,total_GCG_mob] = complxity_and_mobility_segmented(all_GCG_traces(:,start:stop));

GCG_cmplx(i)=median(abs(total_GCG_comlpx));
GCG_mob(i)=median(abs(total_GCG_mob));

Mechanical_Energy_Complexity_Dispersion(i)=std(total_GCG_comlpx).^2+std(total_SCG_comlpx).^2;
Mechanical_Energy_Mobility_Dispersion(i)=std(total_GCG_mob).^2+std(total_SCG_mob).^2;
Total_Energy_Complxity(i)=sqrt(GCG_cmplx(i).^2+SCG_cmplx(i).^2);
Total_Energy_Mobility(i)=sqrt(SCG_mob(i).^2+GCG_mob(i).^2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    curPos = curPos + step;
    
end



Trans_Energy_complexity=(SCG_cmplx);
Angular_Energy_complexity=(GCG_cmplx);
Trans_Energy_mobility=(SCG_mob);
Angular_Energy_mobility=(GCG_mob);

features_out=[Trans_Energy_complexity;Trans_Energy_mobility;Angular_Energy_complexity;Angular_Energy_mobility;Total_Energy_Complxity;Total_Energy_Mobility];
end



