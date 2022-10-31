function [data_traces,data_samplingfreqs] = cycle_extractor(data,fs)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
if isfield(data,'gyroX')
    [~,ploc_ECG] = averagedpeaks(detrend(-data.ekg)',fs); % ORIG
    % ploc_ECG=ploc_ECG(2:end)-0;
    
    aX  = data.accX;
    aY  = data.accY;
    aZ  = data.accZ;
    gX  = data.gyroX;
    gY  = data.gyroY;
    gZ  = data.gyroZ;
    
    %%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
    [alltraces_ax,~,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aX),fs, ploc_ECG);
    [alltraces_ay,~,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aY),fs,ploc_ECG);
    [alltraces_az,~,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZ),fs,ploc_ECG);
    [alltraces_gx,~,FS_gyrx] = averagedpeaks_ECG_REF_mj(detrend(gX),fs,ploc_ECG);
    [alltraces_gy,~,FS_gyry] = averagedpeaks_ECG_REF_mj(detrend(gY ),fs,ploc_ECG);
    [alltraces_gz,~,FS_gyrz] = averagedpeaks_ECG_REF_mj(detrend(gZ),fs,ploc_ECG);
    [alltraces_ekg,~,FS_ekg] = averagedpeaks_ECG_REF_mj(detrend(data.ekg),fs,ploc_ECG);
    
    data_traces=struct('ekg_traces',alltraces_ekg,'accX_traces',alltraces_ax,'accY_traces',alltraces_ay,'accZ_traces',alltraces_az,'gyroX_traces',alltraces_gx,'gyroY_traces',alltraces_gy,'gyroZ_traces',alltraces_gz);
    data_samplingfreqs=struct('ekg_fs',FS_ekg,'accX_fs',FS_accx,'accY_fs',FS_accy,'accZ_fs',FS_accz,'gyroX_fs',FS_gyrx,'gyroY_fs',FS_gyry,'gyroZ_fs',FS_gyrz);
    
else
    [~,ploc_ECG] = averagedpeaks(detrend(-data.ekg)',fs); % ORIG
    % ploc_ECG=ploc_ECG(2:end)-0;
    
    aX  = data.accX;
    aY  = data.accY;
    aZ  = data.accZ;
    
    
    %%%%%%%%%%%%%%%% Ensemble Average and Streching of Cardiac cycles
    [alltraces_ax,~,FS_accx] = averagedpeaks_ECG_REF_mj(detrend(aX),fs, ploc_ECG);
    [alltraces_ay,~,FS_accy] = averagedpeaks_ECG_REF_mj(detrend(aY),fs,ploc_ECG);
    [alltraces_az,~,FS_accz] = averagedpeaks_ECG_REF_mj(detrend(aZ),fs,ploc_ECG);
    [alltraces_ekg,~,FS_ekg] = averagedpeaks_ECG_REF_mj(detrend(data.ekg),fs,ploc_ECG);
    
    data_traces=struct('ekg_traces',alltraces_ekg,'accX_traces',alltraces_ax,'accY_traces',alltraces_ay,'accZ_traces',alltraces_az);
    data_samplingfreqs=struct('ekg_fs',FS_ekg,'accX_fs',FS_accx,'accY_fs',FS_accy,'accZ_fs',FS_accz);
end



end

