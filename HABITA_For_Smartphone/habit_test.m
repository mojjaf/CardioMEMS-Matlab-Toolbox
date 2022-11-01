close all
clear all

smartphoneread=1;
j=1;
if smartphoneread
    
    for ind = 1:15
        close all;
        try
        %     plotandroidfile
        count(j)=j;
        AFib=0;
        data=import_smartphone_data( ind );
        fs=200;
        gr=1;
        coef=[1 1]
        sensor=1 % 1= acc 2 = gyro
        %HABITA
        [AO_amp_raw,AO_i_raw,heartrate_avg,t]=HABITA(data,fs,coef,sensor,gr);
        sensor=2;
        [AO_amp_raw,AO_i_raw,heartrate_avg,t]=HABITA(data,fs,coef,sensor,gr);
                
        % HURNANEN-PAN-TOMPKINS
        [AO_amp_acc,AO_i_acc,HR_avg_acc,t]=Modified_PanTompkins(data,fs,1,1);
        [AO_amp_gyr,AO_i_gyr,HR_avg_gyr,t]=Modified_PanTompkins(data,fs,2,1);
                pause;
        end
        
         
    end
   
end