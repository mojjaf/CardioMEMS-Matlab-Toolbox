% close all;
% clear all;
% filename_scg = 'Acc_Tue Aug 25 125439 GMT+0300 2015.txt';%%TYKS AF
% filename_gcg ='Gyro_Tue Aug 25 125439 GMT+0300 2015.txt';
% % filename_scg = 'Acc_Tue Aug 25 141818 GMT+0300 2015.txt';%TYKS 
% % filename_gcg ='Gyro_Tue Aug 25 141818 GMT+0300 2015.txt';
% % filename_scg = 'Acc_Wed Aug 26 133517 GMT+0300 2015.txt';%%TYKS
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

delimiterIn = ' ';
headerlinesIn = 2;
A = importdata(filename_scg,delimiterIn,headerlinesIn);
B=importdata(filename_gcg,delimiterIn,headerlinesIn);

accX=A.data(:,1);
accY=A.data(:,2);
accZ=A.data(:,3);
GyrX=B.data(:,1);
GyrY=B.data(:,2);
GyrZ=B.data(:,3);
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
%Signals=[accX,accY,accZ,GyrX,GyrY,GyrZ];
Signals=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',GyrX,'gyroY',GyrY,'gyroZ',GyrZ);
fs=200;
[AO_amp_acc,AO_i_acc,HR_avg_acc,t]=Modified_PanTompkins(Signals,fs,1,1);
[AO_amp_gyr,AO_i_gyr,HR_avg_gyr,t]=Modified_PanTompkins(Signals,fs,2,1);


% T_HR = table(heartrate_acc,heartrate_gyr,...
%     'VariableNames',{'SCG_HR','GCG_HR'})
