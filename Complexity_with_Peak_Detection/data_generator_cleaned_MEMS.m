
hfiles = '\\utu.fi\verkkolevyt\MEDICALMEMS\Clinical_Data\STEMI\EXCEL_STEMI_ver_0_2\DATA_VER_0_2_AUTOCLEANED/';
filepattern=fullfile(hfiles,'*.mat')
matfiles=dir(filepattern)

    
    
 for k=1:length(matfiles)
       %             try

        baseFileName = matfiles(k).name;
        data=load([fullfile(hfiles, baseFileName)]);
        CVD_status=data.data_cleaned.info.STEMI_TYPE;
    if data.data_cleaned.info.USABLE_ECG ==1    
%         
%         fs_lsm=data.data_cleaned.info.SAMPLING_FREQ.others;
%         
%         gyroX=data.data_cleaned.gyroX';
%         gyroY=data.data_cleaned.gyroY';
%         gyroZ=data.data_cleaned.gyroZ';
%         accX=data.data_cleaned.accX_lsm';
%         accY=data.data_cleaned.accY_lsm';
%         accZ=data.data_cleaned.accZ_lsm';
%         ecg=data.data_cleaned.ekg';
%         
%         accdata=[accX,accY,accZ];
%         gyrodata=[gyroX,gyroY,gyroZ];
        if CVD_status==1
       
        full_name='asro_STEMI_';
%         save([full_name num2str(k)], 'accdata','gyrodata','ecg');
        save([full_name num2str(k)], 'data');

        elseif CVD_status==2
       
        full_name='asro_PostSTEMI_';
%         save([full_name num2str(k)], 'accdata','gyrodata','ecg')
         save([full_name num2str(k)], 'data')
        elseif CVD_status==3

         full_name='asro_control_';
%         save([full_name num2str(k)], 'accdata','gyrodata','ecg');
          save([full_name num2str(k)], 'data');
        end
    end
 end

% %%% 
%  for k=1:length(matfiles)
%        %             try
% 
%         baseFileName = matfiles(k).name;
%         data=load([fullfile(hfiles, baseFileName)]);
%         CVD_status=data.data_cleaned.info.STEMI_TYPE;
%     if data.data_cleaned.info.USABLE_ECG    
%         
%         fs_lsm=data.data_cleaned.info.SAMPLING_FREQ.others;
%         fs_adxl=data.data_cleaned.info.SAMPLING_FREQ.ADXL355_hf;
%  
%         accX_hf=fft_filter(data.data_cleaned.ADXL355_hf.accX(1:jitter*fs_adxl),fs_adxl,10,80);
%         accY_hf=fft_filter(data.data_cleaned.ADXL355_hf.accY(1:jitter*fs_adxl),fs_adxl,10,80);
%         accZ_hf=fft_filter(data.data_cleaned.ADXL355_hf.accZ(1:jitter*fs_adxl),fs_adxl,10,80);
%         accX_lf=fft_filter(data.data_cleaned.ADXL355_hf.accX(1:jitter*fs_adxl),fs_adxl,1,10);
%         accY_lf=fft_filter(data.data_cleaned.ADXL355_hf.accY(1:jitter*fs_adxl),fs_adxl,1,10);
%         accZ_lf=fft_filter(data.data_cleaned.ADXL355_hf.accZ(1:jitter*fs_adxl),fs_adxl,1,10);
%         ekg_hf=upsample(ekg_lf,2);
%         ekg_hf=ekg_hf(1:jitter*fs_adxl);
%                 
%         %         len=min([numel(gyroX) numel(accX)]);
% %         accdata=[accX,accY,accZ];
% %         gyrodata=[gyroX,gyroY,gyroZ];
% %         full_name='asro_control_';
% %         save([full_name num2str(k)], 'accdata','gyrodata','ecg');