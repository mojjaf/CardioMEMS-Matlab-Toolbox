function [averaged_signal_lsm,averaged_signal_adxl] = ensemble_morph_sim_9axis(data)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
%  if isfield(data,'accX_adxl')
    fs_lsm=200;
        fs_adxl=400;
        jitter=floor(min(length(data.accX)/fs_lsm,length(data.accX_adxl)/fs_adxl));
        
        gyroX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,3,30,0);
        gyroY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,3,30,0);
        gyroZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,3,30,0);
        accX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,5,30,0);
        accY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,5,30,0);
        accZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,5,30,0);
        ekg_lf=fft_filter(data.ekg(1:jitter*fs_lsm),fs_lsm,1,49,0);
        
        
        [avpeak,ploc_ECG] = averagedpeaks(ekg_lf,fs_lsm);
[alltraces_ax,avpeak_accx,~] = averagedpeaks_ECG_REF_mj(accX,fs_lsm,3, ploc_ECG);
[alltraces_ay,avpeak_accy,~] = averagedpeaks_ECG_REF_mj(accY,fs_lsm,3, ploc_ECG);
[alltraces_az,avpeak_accz,~] = averagedpeaks_ECG_REF_mj(accZ,fs_lsm,3, ploc_ECG);
[alltraces_gx,avpeak_gyrx,~] = averagedpeaks_ECG_REF_mj(gyroX,fs_lsm,3, ploc_ECG);
[alltraces_gy,avpeak_gyry,~] = averagedpeaks_ECG_REF_mj(gyroY,fs_lsm,3, ploc_ECG);
[alltraces_gz,avpeak_gyrz,~] = averagedpeaks_ECG_REF_mj(gyroZ,fs_lsm,3, ploc_ECG);
averaged_signal_lsm=struct('Ens_accX',avpeak_accx,'Ens_accY',avpeak_accy,'Ens_accZ',avpeak_accz,'Ens_gyroX',avpeak_gyrx,'Ens_gyroY',avpeak_gyry,'Ens_gyroZ',avpeak_gyrz);


        accX_hf=fft_filter(data.accX_adxl(1:jitter*fs_adxl),fs_adxl,30,80,0);
        accY_hf=fft_filter(data.accY_adxl(1:jitter*fs_adxl),fs_adxl,30,80,0);
        accZ_hf=fft_filter(data.accZ_adxl(1:jitter*fs_adxl),fs_adxl,30,80,0);
        accX_lf=fft_filter(data.accX_adxl(1:jitter*fs_adxl),fs_adxl,5,30,0);
        accY_lf=fft_filter(data.accY_adxl(1:jitter*fs_adxl),fs_adxl,5,30,0);
        accZ_lf=fft_filter(data.accZ_adxl(1:jitter*fs_adxl),fs_adxl,5,30,0);
        ekg_hf=fft_filter(data.ekg_hf(1:jitter*fs_adxl),400,1,49,0);
                
%         data_adxl=struct('ekg',ekg_hf,'accX',accX_hf,'accY',accY_hf,'accZ',accZ_hf,'accXLF',accX_lf,'accYLF',accY_lf,'accZLF',accZ_lf); 
     
     [avpeak_ecgHF,ploc_ECGHF] = averagedpeaks(ekg_hf',fs_adxl);


[alltraces_axHF,avpeak_accxHF,~] = averagedpeaks_ECG_REF_mj(accX_hf,fs_adxl,3, ploc_ECGHF);
[alltraces_ayHF,avpeak_accyHF,~] = averagedpeaks_ECG_REF_mj(accY_hf,fs_adxl,3, ploc_ECGHF);
[alltraces_azHF,avpeak_acczHF,~] = averagedpeaks_ECG_REF_mj(accZ_hf,fs_adxl,3, ploc_ECGHF);
[alltraces_axLF,avpeak_accxLF,~] = averagedpeaks_ECG_REF_mj(accX_lf,fs_adxl,3, ploc_ECGHF);
[alltraces_ayLF,avpeak_accyLF,~] = averagedpeaks_ECG_REF_mj(accY_lf,fs_adxl,3, ploc_ECGHF);
[alltraces_azLF,avpeak_acczLF,~] = averagedpeaks_ECG_REF_mj(accZ_lf,fs_adxl,3, ploc_ECGHF);

averaged_signal_adxl=struct('Ens_accX_HF',avpeak_accxHF,'Ens_acc_YHF',avpeak_accyHF,'Ens_acc_ZHF',avpeak_acczHF,'Ens_accX_LF',avpeak_accxLF,'Ens_accY_LF',avpeak_accyLF,'Ens_accZ_LF',avpeak_acczLF);

     
%  else
%      
% [avpeak,ploc_ECG] = averagedpeaks(detrend(-data.ekg)',fs);
% [alltraces_ax,avpeak_accx,ploc_accx] = averagedpeaks_ECG_REF(detrend(data.accX),fs,3, ploc_ECG);
% [alltraces_ay,avpeak_accy,ploc_accy] = averagedpeaks_ECG_REF(detrend(data.accY),fs,3, ploc_ECG);
% [alltraces_az,avpeak_accz,ploc_accz] = averagedpeaks_ECG_REF(detrend(data.accZ),fs,3, ploc_ECG);
% [alltraces_gx,avpeak_gyrx,ploc_gyrx] = averagedpeaks_ECG_REF(detrend(data.gyroX),fs,3, ploc_ECG);
% [alltraces_gy,avpeak_gyry,ploc_gyry] = averagedpeaks_ECG_REF(detrend(data.gyroY ),fs,3, ploc_ECG);
% [alltraces_gz,avpeak_gyrz,ploc_gyrz] = averagedpeaks_ECG_REF(detrend(data.gyroZ),fs,3, ploc_ECG);
% 
% 
% averaged_signal=struct('Ens_accX',avpeak_accx,'Ens_accY',avpeak_accy,'Ens_accZ',avpeak_accz,'Ens_gyroX',avpeak_gyrx,'Ens_gyroY',avpeak_gyry,'Ens_gyroZ',avpeak_gyrz);

 end








