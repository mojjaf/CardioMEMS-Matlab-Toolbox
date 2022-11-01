
% clear all, close all;

for ind = 1:5;
    close all;
%     plotandroidfile

  Signals=load_SPfiles( ind )
%     Signals=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyroX,'gyroY',gyroY,'gyroZ',gyroZ);
Fs=200

coef=[1 1];


[AO_amp_gyr,AO_i_gyr,heartrate_gyr,t]=APD2(Signals,Fs,coef,2,1);

[AO_amp_acc,AO_i_acc,heartrate_acc,t]=APD2(Signals,Fs,coef,1,1);
[AO_amp_acc,AO_i_acc,HR_avg_acc,t]=Modified_PanTompkins(Signals,Fs,1,1);
[AO_amp_gyr,AO_i_gyr,HR_avg_gyr,t]=Modified_PanTompkins(Signals,Fs,2,1);

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
pause;
end