function R = load_measurements(sensor,ind)
%sensor == 0 : accelerometer, sensor == 1: gyroscope
%ind > 0 : filenumber

addpath('.\thtools\')
% hfiles = findmeasfiles('./omat/');
% hfiles = findmeasfiles('./rritest/');
% hfiles = findmeasfiles('./Sony Xperia Z1 Compact'); % Puhelimen valinta
hfiles = findmeasfiles('./SUMM+proto'); % Puhelimen valinta
filnum = ind; % TÄTÄ MUUTTAMALLA VAIHTAA MITTAUSTA! 

hfile = hfiles(filnum);
% fs = 200;
hfile = hfile{1}
tt = findstr(hfile,'\');

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

gyroX = grdata(:,1);
gyroY = grdata(:,2);
gyroZ = grdata(:,3);

accX = accdata(:,1);
accY = accdata(:,2);
accZ = accdata(:,3);

%<filtering>
[b,a]=butter(4,[.01,.4]);        % Bandpass digital filter design
gyroX = filter(b,a,gyroX);
gyroY = filter(b,a,gyroY);
gyroZ = filter(b,a,gyroZ);
accX = filter(b,a,accX);
accY = filter(b,a,accY);
accZ = filter(b,a,accZ);
%</filtering>

if(sensor==0)
    R = zeros(4,length(timeacc));
    R(1,:) = timeacc;
    R(2,:) = accX';
    R(3,:) = accY';
    R(4,:) = accZ';
end

if(sensor==1)
    R = zeros(4,length(timegr));
    R(1,:) = timegr;
    R(2,:) = gyroX';
    R(3,:) = gyroY';
    R(4,:) = gyroZ';
end