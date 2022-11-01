close all;
clear all;
% filename_scg = 'Acc_Tue Aug 25 125439 GMT+0300 2015.txt';%%TYKS AF
% filename_gcg ='Gyro_Tue Aug 25 125439 GMT+0300 2015.txt';
% filename_scg = 'Acc_Tue Aug 25 141818 GMT+0300 2015.txt';%TYKS 
% filename_gcg ='Gyro_Tue Aug 25 141818 GMT+0300 2015.txt';
% filename_scg = 'Acc_Wed Aug 26 133517 GMT+0300 2015.txt';%%TYKS
% filename_gcg ='Gyro_Wed Aug 26 133517 GMT+0300 2015.txt';

% filename_scg = 'Acc_Mon Aug 24 10_47_06 GMT+03_00 2015.txt';%% bad signal
% filename_gcg ='Gyro_Mon Aug 24 10_47_06 GMT+03_00 2015.txt';

% filename_scg = 'Acc_Wed Mar 11 140802 GMT+0000 2015.txt';
% filename_gcg ='Gyro_Wed Mar 11 140802 GMT+0000 2015.txt';
% filename_scg = 'Acc_Thu Mar 12 173636 GMT+0000 2015.txt';
% filename_gcg ='Gyro_Thu Mar 12 173636 GMT+0000 2015.txt';
% filename_scg = 'Acc_Mon Mar 16 104138 GMT+0000 2015.txt';
% filename_gcg ='Gyro_Mon Mar 16 104138 GMT+0000 2015.txt';
%%%%%%%%%%%%%%%%%%%%%%%%
% filename_scg = 'Acc_Tue Aug 11 10_26_32 GMT+03_00 2015.txt';
% filename_gcg ='Gyro_Tue Aug 11 10_26_32 GMT+03_00 2015.txt';
% filename_scg = 'Acc_Mon Aug 24 10_47_06 GMT+03_00 2015.txt';
% filename_gcg ='Gyro_Mon Aug 24 10_47_06 GMT+03_00 2015.txt';
% filename_scg = 'Acc_Tue Sep 08 11_35_59 GMT+03_00 2015.txt';
% filename_gcg ='Gyro_Tue Sep 08 11_35_59 GMT+03_00 2015.txt';
% filename_gcg ='Gyro_Wed Aug 26 133517 GMT+0300 2015.txt';
% filename_scg = 'Acc_Wed Mar 11 140802 GMT+0000 2015.txt';
% filename_gcg ='Gyro_Wed Mar 11 140802 GMT+0000 2015.txt';
filename_scg ='Acc_Mon Jul 27 144741 GMT+0300 2015.txt';
filename_gcg ='Gyro_Mon Jul 27 144741 GMT+0300 2015.txt';
% filename_scg ='Acc_Mon Jul 27 150812 GMT+0300 2015.txt';
% filename_gcg ='Gyro_Mon Jul 27 150812 GMT+0300 2015.txt';

% 
% delimiterIn = ' ';
% headerlinesIn = 2;
% A = importdata(filename_scg,delimiterIn,headerlinesIn);
% B=importdata(filename_gcg,delimiterIn,headerlinesIn);
[accX,accY,accZ,timestamp] = ImportSmartPhoneData(filename_scg);
[GyrX,GyrY,GyrZ,timestamp] = ImportSmartPhoneData(filename_gcg);

% accX=A.data(:,1);
% accY=A.data(:,2);
% accZ=A.data(:,3);
% GyrX=B.data(:,1);
% GyrY=B.data(:,2);
% GyrZ=B.data(:,3);
if length(accX)<length(GyrY)
    len=length(accX);
else
    len=length(GyrY);
end
accX=accX(1:len);
accY=accY(1:len);
accZ=accZ(1:len);
GyrX=GyrX(1:len);
GyrY=GyrY(1:len);
GyrZ=GyrZ(1:len);
% Signals=struct('accX',accX(1000:6000),'accY',accY(1000:6000),'accZ',accZ(1000:6000),...
%     'gyroX',GyrX(1000:6000),'gyroY',GyrY(1000:6000),'gyroZ',GyrZ(1000:6000)); % only for dogs
Signals=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',GyrX,'gyroY',GyrY,'gyroZ',GyrZ);
Fs=200

coef=[1 1];

[AO_amp_gyr,AO_i_gyr,heartrate_gyr,t]=APD2(Signals,Fs,coef,2,1);

[AO_amp_acc,AO_i_acc,heartrate_acc,t]=APD2(Signals,Fs,coef,1,1);

[p_g_a,s_g_a,f_g_a] = Numeval( AO_i_gyr, AO_i_acc);

[p_a_g,s_a_g,f_a_g] = Numeval( AO_i_acc, AO_i_gyr);

if f_a_g<0.90 && f_g_a <0.90
    
    coef=[1 0];
            gr=1;

    [AO1_amp_gyr,AO_1_gyr,heartrate_gyr1,t]=APD2(Signals,Fs,coef,2,gr);
    
    [AO1_amp_acc,AO_1_acc,heartrate_acc1,t]=APD2(Signals,Fs,coef,1,gr);
    
    [p1_g_a,s_g_a,f1_g_a] = Numeval( AO_1_gyr, AO_1_acc);
    
    [p1_a_g,s_a_g,f1_a_g] = Numeval( AO_1_acc, AO_1_gyr);
    
    AO_amp_gyr=AO1_amp_gyr;
    
    AO_i_gyr=AO_1_gyr;
    
    heartrate_gyr=heartrate_gyr1;
    
    AO_amp_acc=AO1_amp_acc;
    
    AO_i_acc=AO_1_acc;
    
    heartrate_acc=heartrate_acc1;
    
    if f1_g_a<0.90 && f1_a_g <0.90
        
        coef=[0 1];
                gr=1;
        
        [AO2_amp_gyr,AO_2_gyr,heartrate_gyr2,t]=APD2(Signals,Fs,coef,2,gr);
        
        [AO2_amp_acc,AO_2_acc,heartrate_acc2,t]=APD2(Signals,Fs,coef,1,gr);
        
        [p2_g_a,s_g_a,f2_g_a] = Numeval( AO_2_gyr, AO_2_acc);
        
        [p2_a_g,s_a_g,f2_a_g] = Numeval( AO_2_acc, AO_2_gyr);
        
        AO_amp_gyr=AO2_amp_gyr;
        
        AO_i_gyr=AO_2_gyr;
        
        heartrate_gyr=heartrate_gyr2;
        
        AO_amp_acc=AO2_amp_acc;
        
        AO_i_acc=AO_2_acc;
        
        heartrate_acc=heartrate_acc2;
        
    else
        gr=1;
        [AO_amp_gyr,AO_i_gyr,heartrate_gyr]=APD2(Signals,Fs,coef,2,gr);
        
        [AO_amp_acc,AO_i_acc,heartrate_acc]=APD2(Signals,Fs,coef,1,gr);
        
    end
    
end


[AO_amp_acc,AO_i_acc,HR_avg_acc,t]=Modified_PanTompkins(Signals,200,1,1);
[AO_amp_gyr,AO_i_gyr,HR_avg_gyr,t]=Modified_PanTompkins(Signals,200,2,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b2b_acc=diff(AO_i_acc);
b2b_gyr=diff(AO_i_gyr);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_HR = table(heartrate_acc,heartrate_gyr,...
    'VariableNames',{'SCG_HR','GCG_HR'})
Read=0;
if Read
savdir1 ='D:\old computer\Signal Processing Database\6axis analysis\Mobile Measurements\ithlete_measurements_with_6axis\Kubios HRV\Raw RR intervals\Acc' ;
savdir2 ='D:\old computer\Signal Processing Database\6axis analysis\Mobile Measurements\ithlete_measurements_with_6axis\Kubios HRV\Raw RR intervals\Gyro'   ;   
    
                             b2b_acc=(b2b_acc'/Fs);
                             b2b_gyr=(b2b_gyr'/Fs);
                             
                             
                save(fullfile(savdir1,sprintf('AAscg_hil_sub_%02d.txt',15)),'b2b_acc','-ascii');
                save(fullfile(savdir2,sprintf('AAgcg_hil_sub_%02d.txt',15)),'b2b_gyr','-ascii');
                
end

