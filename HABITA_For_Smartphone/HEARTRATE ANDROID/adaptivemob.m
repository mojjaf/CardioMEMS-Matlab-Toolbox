close all;
clear all;
filename_scg = 'Acc_Wed Mar 11 140802 GMT+0000 2015.txt';
filename_gcg ='Gyro_Wed Mar 11 140802 GMT+0000 2015.txt';
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
Signals=[accX,accY,accZ,GyrX,GyrY,GyrZ];
Fs=200

[AO_amp_acc,AO_i_acc,pk,ekg_locs]=APD(double(Signals),Fs,1,len,1,1);

[AO_amp_gyr,AO_i_gyr,pk,ekg_locs]=APD(double(Signals),Fs,1,len,2,1);


b2b_acc=diff(AO_i_acc)/Fs;
avg_HR_acc=60*(1/mean(b2b_acc));
b2b_gyr=diff(AO_i_gyr)/Fs;
avg_HR_gyr=60*(1/mean(b2b_gyr));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

T_HR = table(avg_HR_acc,avg_HR_gyr,...
    'VariableNames',{'SCG_HR','GCG_HR'})


