close all;
clear all;
filename_scg = 'Acc_Wed Mar 11 140802 GMT+0000 2015.txt';
filename_gcg ='Gyro_Wed Mar 11 140802 GMT+0000 2015.txt';
% filename_scg ='Acc_2015-07-27 15_06_52.txt';
% filename_gcg ='Gyro_2015-07-27 15_06_52.txt';


delimiterIn = ' ';
headerlinesIn = 2;
A = importdata(filename_scg,delimiterIn,headerlinesIn);
B=importdata(filename_gcg,delimiterIn,headerlinesIn);
Fs=200

accZ=A.data(:,1);

GyrY=B.data(:,1);

if length(accZ)<length(GyrY)
    len=length(accZ);
else
    len=length(GyrY);
end
% accZ_b=accZ(1:len);
% accZ_a=detrend(accZ(1:len));
samples=1:numel(accZ(1:len));
t=samples/Fs;
% accZ=cumtrapz(t,accZ(1:len));
accZ=accZ(1:len);
% GyrY=cumtrapz(t,GyrY(1:len));
GyrY=GyrY(1:len);
figure
ax(1)=subplot(211),plot((HF_sphone(detrend(accZ(1:len)),1)));
ax(2)=subplot(212),plot(HF_sphone(GyrY,1));
linkaxes([ ax(2) ax(1)],'x');
%Signals=[accX,accY,accZ,GyrX,GyrY,GyrZ];
% Signals=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',GyrX,'gyroY',GyrY,'gyroZ',GyrZ);
% fs=200;
% [AO_amp_acc,AO_i_acc,HR_avg_acc,t]=Modified_PanTompkins(Signals,fs,1,1);
% [AO_amp_gyr,AO_i_gyr,HR_avg_gyr,t]=Modified_PanTompkins(Signals,fs,2,1);