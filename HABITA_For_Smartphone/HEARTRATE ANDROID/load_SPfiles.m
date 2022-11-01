function [ Signals ] = load_SPfiles( ind )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% clear all
% close all

addpath('.\thtools2\')
% hfiles = findmeasfiles('./surgery/');
% hfiles = findmeasfiles('./rritest/');
% hfiles = findmeasfiles('./Sony Xperia Z1 Compact'); % Puhelimen valinta
% hfiles = findmeasfiles('./SUMM+proto'); % Puhelimen valinta
 hfiles = findmeasfiles('./Ischemia data');
%  hfiles = findmeasfiles('./HRV_Smartphone_Sep2015');
filnum = ind % TÄTÄ MUUTTAMALLA VAIHTAA MITTAUSTA! 

hfile = hfiles(filnum);
% fs = 200;
hfile = hfile{1}
tt = findstr(hfile,'\')

newmanpoints = 0;

infofilename = hfile(tt(end)+1:end)
pathname = hfile(1:tt(end))

gyrofilename = [ 'Gyro', infofilename(5:end)  ]
accfilename =  [ 'Acc', infofilename(5:end)  ]
rrifilename =  [ 'RRInt', infofilename(5:end)  ]

grdata = importdata([pathname,gyrofilename]);

accdata = importdata([pathname,accfilename]);

timegr = grdata(:,4);
timeacc = accdata(:,4);
tmin = min(min(timegr),min(timeacc));

timegr = timegr-tmin;
timeacc = timeacc-tmin;
fs = 10*round(1e8*1/mean(diff(timeacc)))

[b,a]=butter(4,[.01,.4]);  %0.01-0.4      % Bandpass digital filter design
%     h = fvtool(b,a);                 % Visualize filter
 
% b=1;
% a=1;

gyroX = grdata(:,1);
gyroY = grdata(:,2);
gyroZ = grdata(:,3);

accX = accdata(:,1);
accY = accdata(:,2);
accZ = accdata(:,3);


gyroX = filter(b,a,gyroX);
gyroY = filter(b,a,gyroY);
gyroZ = filter(b,a,gyroZ);
accX = filter(b,a,accX);
accY = filter(b,a,accY);
accZ = filter(b,a,accZ);

% gyroX = fft_filter(gyroX,fs,.5,40);
% gyroY = fft_filter(gyroY,fs,.5,40);
% gyroZ = fft_filter(gyroZ,fs,.5,40);
% accX = fft_filter(accX,fs,.5,40);
% accY = fft_filter(accY,fs,.5,40);
% accZ = fft_filter(accZ,fs,.5,40);

figure(11)
sb(1) = subplot(311);
plot(gyroX)
legend({'GyroX'})
sb(2) = subplot(312);
plot(gyroY,'r')
legend({'GyroY'})
sb(3) = subplot(313);
plot(gyroZ,'c')
legend({'GyroZ'})
linkaxes(sb,'x');

figure(12)
sb(1) = subplot(311);
plot(accX)
legend({'AccX'})
sb(2) = subplot(312);
plot(accY,'r')
legend({'AccY'})
sb(3) = subplot(313);
plot(accZ,'c')
legend({'AccZ'})
linkaxes(sb,'x');

header = load_header2([pathname,infofilename])

rriAX = xcorrsykelyhyt(accX,fs);
rriAY = xcorrsykelyhyt(accY,fs);
rriAZ = xcorrsykelyhyt(accZ,fs);
rriGX = xcorrsykelyhyt(accX,fs);
rriGY = xcorrsykelyhyt(accY,fs);
rriGZ = xcorrsykelyhyt(accZ,fs);

rric = combinerrint(rriAZ,rriGY);
rric = combinerrint(rric,rriGX);
rric = combinerrint(rric,rriGZ);
rric = combinerrint(rric,rriAX);
rric = combinerrint(rric,rriAY);

% figure(3)
% plot(timeacc(1:200))
% hold on
% plot(timegr(1:200),'r')
% hold off

% figure(4)
% subplot(211)
% plot(diff(timeacc))
% subplot(212)
% plot(diff(timegr))

figure(5)

subplot(321)
plot(rriAX)
title rriAX
subplot(322)
plot(rriAY)
title rriAY
subplot(323)
plot(rriAZ)
title rriAZ
subplot(324)
plot(rriGX)
title rriGX
subplot(325)
plot(rriGY)
title rriGY
subplot(326)
plot(rriGZ)
title rriGZ

figure(6)
% plot(abs(fft(accZ)))
plot(xcorr(accZ(1000:2000)))

figure(7) 
plot(rric)

try
rridata = importdata([pathname,rrifilename]);

rriandroid = rridata(:,2);
figure(8) 
plot(rriandroid)
hold on

plot(rric,'r--')
hold off

end

% figure(9)
% plot(rriandroid,'r')
% hold on
% plot(rriAZ)
% hold off
Signals=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyroX,'gyroY',gyroY,'gyroZ',gyroZ);

end

