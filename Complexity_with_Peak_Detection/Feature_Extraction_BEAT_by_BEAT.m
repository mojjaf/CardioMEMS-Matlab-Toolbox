clear all
close all
for i=1:2
    switch i
        case 1
            
            % %             NSTEMI / STEMI  
             hfiles = 'D:\Seafile\Seafile\Saeed_And_Mojtaba_Underwork\MI_ASRO\NEW CLEANED _ 15012019\STEMI_HQekg_noAF/';
          
            Is_MI=1;
            filepattern=fullfile(hfiles,'*.mat')
            matfiles=dir(filepattern)
            
        case 2
            %             all controls
                        hfiles = 'D:\Seafile\Seafile\Saeed_And_Mojtaba_Underwork\MI_ASRO\NEW CLEANED _ 15012019\CONTROL_HQekg_noAF/';
            
            Is_MI=0;
            filepattern=fullfile(hfiles,'*.mat')
            matfiles=dir(filepattern)
        case 3
            Is_MI=-1;
            
        otherwise
            disp('Unknown path.')
            
    end
    count=0
    
    for k=1:length(matfiles)
        %             try
%      for k=1:1
        fs=200; % sampling freqeuncy
        gr=0; %graphic
        filt=1;  % 1= normal BPF 2= Envelope 3= SSA
        %           data=import_smartphone_data(hfiles, k,filt,fs, gr );% load smartphone data
        baseFileName = matfiles(k).name;
        load([fullfile(hfiles, baseFileName)]);
        %         tt=(1:numel(data.accX))/fs;
        
        
        fs_lsm=200;
        fs_adxl=400;
        jitter=floor(min(length(data.accX)/fs_lsm,length(data.accX_adxl)/fs_adxl));
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LOW FREQ DATA
        gyroX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,3,30);
        gyroY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,3,30);
        gyroZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,3,30);
        accX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,5,30);
        accY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,5,30);
        accZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,5,30);
        ekg_lf=fft_filter(data.ekg(1:jitter*fs_lsm),fs_lsm,1,49);
        
        data_lsm=struct('ekg',ekg_lf,'accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyroX,'gyroY',gyroY,'gyroZ',gyroZ);
        
        [data_traces_LF,data_samplingfreqs_LF] = cycle_extractor(data,fs);
         save(['LSM_6axis_cycles' baseFileName '.mat'],'data_traces_LF')
                 save(['LSM_6axis_FS' baseFileName '.mat'],'data_samplingfreqs_LF')

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% HIGH FREQ DATA
        
        accX_hf=fft_filter(data.accX_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        accY_hf=fft_filter(data.accY_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        accZ_hf=fft_filter(data.accZ_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        ekg_hf=fft_filter(data.ekg_hf(1:jitter*fs_adxl),400,1,49);
                
        data_adxl=struct('ekg',ekg_hf,'accX',accX_hf,'accY',accY_hf,'accZ',accZ_hf);
        [data_traces_HF,data_samplingfreqs_HF] = cycle_extractor(data_adxl,fs);
        
        
         save(['ADXL_3axis_cycles_' baseFileName '.mat'],'data_traces_HF')
                 save(['ADXL_3axis_FS_' baseFileName '.mat'],'data_samplingfreqs_HF')
                 
                 count=count+size(data_traces_LF.accX_traces,1)

%         trace_signals=[data_traces_LF,data_samplingfreqs_LF,data_traces_HF,data_samplingfreqs_HF];
      

        %        adxl_energy(data_adxl,fs_adxl,1,k,Is_MI)
       close all
%          pause
%          continue

    end
end