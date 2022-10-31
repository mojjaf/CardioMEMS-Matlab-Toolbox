function [ features,features_avg,features_fusion_median,labID_ext] = Feature_Extraction_MI_Complexity_HIGHFREQ_ENVELOPE(data,fs,win,step,lab, id)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

    
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

[avpeak,ploc_ECG] = averagedpeaks(detrend(abs(data.ekg(curPos:curPos+windowLength-1)))',fs); % ORIG
% ploc_ECG=ploc_ECG(2:end)-0;


% ploc_ECG=ploc_ECG(2:end)-0;

    aXHF  = data.accX(curPos:curPos+windowLength-1);
    aYHF  = data.accY(curPos:curPos+windowLength-1);
    aZHF  = data.accZ(curPos:curPos+windowLength-1);
    aXLF  = data.accXLF(curPos:curPos+windowLength-1);
    aYLF  = data.accYLF(curPos:curPos+windowLength-1);
    aZLF  = data.accZLF(curPos:curPos+windowLength-1);
   
%%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
[alltraces_axHF,avpeak_accx,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aXHF),fs, ploc_ECG); 
[alltraces_ayHF,avpeak_accy,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aYHF),fs, ploc_ECG); 
[alltraces_azHF,avpeak_accz,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZHF),fs, ploc_ECG); 
[alltraces_axLF,avpeak_gyrx,FS_accxLF] = averagedpeaks_ECG_REF_mj(detrend(aXLF),fs, ploc_ECG); 
[alltraces_ayLF,avpeak_gyry,FS_accyLF] = averagedpeaks_ECG_REF_mj(detrend(aYLF ),fs, ploc_ECG); 
[alltraces_azLF,avpeak_gyrz,FS_acczLF] = averagedpeaks_ECG_REF_mj(detrend(aZLF),fs, ploc_ECG); 
%%%%%%%%%%%%%

[up_axhf_en,~] = envelope(alltraces_axHF.^2);
[up_ayhf_en,~] = envelope(alltraces_ayHF.^2);
[up_azhf_en,~] = envelope(alltraces_azHF.^2);
[up_axlf_en,~] = envelope(alltraces_axLF.^2);
[up_aylf_en,~] = envelope(alltraces_ayLF.^2);
[up_azlf_en,~] = envelope(alltraces_azLF.^2);



tot_envprod_enhf=zscore((up_axhf_en).*(up_ayhf_en).*(up_azhf_en));
tot_envsum_enhf=zscore(sqrt((up_axhf_en).^2+(up_ayhf_en).^2+(up_azhf_en).^2));
tot_envprod_enlf=zscore((up_axlf_en).*(up_aylf_en).*(up_azlf_en));
tot_envsum_enlf=zscore(sqrt((up_axlf_en).^2+(up_aylf_en).^2+(up_azlf_en).^2));
tot_env_hflfpro=tot_envprod_enhf+tot_envprod_enlf;
tot_env_hflf_sum=tot_envsum_enhf+tot_envsum_enlf;


% [up_axhf,lo_axhf] = envelope(alltraces_axHF);
% [up_ayhf,lo_ayhf] = envelope(alltraces_ayHF);
% [up_azhf,lo_azhf] = envelope(alltraces_azHF);
% [up_axlf,lo_axlf] = envelope(alltraces_axLF);
% [up_aylf,lo_aylf] = envelope(alltraces_ayLF);
% [up_azlf,lo_azlf] = envelope(alltraces_azLF);
% 
% 
% 
% tot_envprod_hf=median(up_axhf).*median(up_ayhf).*median(up_azhf);
% tot_envsum_hf=sqrt(median(up_axhf).^2+median(up_ayhf).^2+median(up_azhf).^2);
% tot_envprod_lf=median(up_axlf).*median(up_aylf).*median(up_azlf);
% tot_envsum_lf=sqrt(median(up_axlf).^2+median(up_aylf).^2+median(up_azlf).^2);
% 
% tot_envprod_hflo=median(lo_axhf).*median(lo_ayhf).*median(lo_azhf);
% tot_envsum_hflo=sqrt(median(lo_axhf).^2+median(lo_ayhf).^2+median(lo_azhf).^2);
% tot_envprod_lflo=median(lo_axlf).*median(lo_aylf).*median(lo_azlf);
% tot_envsum_lflo=sqrt(median(lo_axlf).^2+median(lo_aylf).^2+median(lo_azlf).^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FS=1000;
 Features_envHFpro(i,:) = MI_FeatureExtraction_FrameWise(tot_envprod_enhf,FS_accx);
 Features_envHFsum(i,:) = MI_FeatureExtraction_FrameWise(tot_envsum_enhf,FS_accx);
 Features_envLFpro(i,:) = MI_FeatureExtraction_FrameWise(tot_envprod_enlf,FS_accxLF);
 Features_envLFsum(i,:)= MI_FeatureExtraction_FrameWise(tot_envsum_enlf,FS_accxLF);
 Features_envHFLFpro(i,:)= MI_FeatureExtraction_FrameWise(tot_env_hflfpro,FS_accx);
 Features_envHFLFsum(i,:) = MI_FeatureExtraction_FrameWise(tot_env_hflf_sum,FS_accxLF);
%     [coeff_tot,scores_tot,latent,~,PTvar_tot]=pca([acc;gyro]);

%     Mechanical_Complexity_Ratios(3,i)=latent(2)/latent(1);

    curPos = curPos + step;
    
end


if lab
labels=ones(size(Features_envLFsum,1),1);
% labels_x=ones(6,1);
labels_x=ones(size(Features_envLFsum,1),1);
 else
   labels=zeros(size(Features_envLFsum,1),1);
   labels_x=zeros(size(Features_envLFsum,1),1);
 end
idx=zeros(size(labels));
idx(:,:)=id;

ID_X=zeros(size(labels_x));
ID_X(:,:)=id;
features=[Features_envHFpro,Features_envHFsum,Features_envLFpro,Features_envLFsum,Features_envHFLFpro,Features_envHFLFsum,labels,idx];
Features_3D(1,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envHFpro;
Features_3D(2,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envHFsum;
Features_3D(3,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envLFpro;
Features_3D(4,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envLFsum;
Features_3D(5,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envHFLFpro;
Features_3D(6,1:size(Features_envHFpro,1),1:size(Features_envHFpro,2))=Features_envHFLFsum;
features_avg=median(Features_3D,1);
% features_avg=[F_gx,F_gy,F_gz,F_ax,F_ay,F_az,labels_x,ID_X];
features_avg = reshape(features_avg,size(features_avg,1)*size(features_avg,2),size(features_avg,3));

features_fusion_median=[median(Features_envHFpro),median(Features_envHFsum),median(Features_envLFpro),median(Features_envLFsum),median(Features_envHFLFpro),median(Features_envHFLFsum)];

labID_ext=[labels,idx];

end

