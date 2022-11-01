function [ output_args ] = plotsmartphonedata( AccX,AccY,AccZ,GyrX,GyrY,GyrZ,Fs,gr )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%

fs=Fs;
samples=1:numel(AccX);
t=samples/fs;

if gr
figure
ax(1)=subplot(3,1,1);plot(t,AccX)
title('Accelerometer')
xlabel('X axis')
ylabel('Acceleration (G)')

ax(2)=subplot(3,1,2);plot(t,AccY)
xlabel('Y axis')
ylabel('Accelerations (G)')

ax(3)=subplot(3,1,3);plot(t,AccZ)
xlabel('Z axis')
ylabel('Accelerations (G)')

linkaxes([ax(3) ax(2) ax(1)],'x');
end
if gr
figure
ax(1)=subplot(2,3,1);plot(t,AccX)
title('Accelerometer')
xlabel('X axis')
ylabel('Acceleration (G)')
ax(2)=subplot(2,3,4);plot(t,GyrX,'r')
title('Gyroscope')
xlabel('X axis')
ylabel('Angular Velocity (rad/s)')

ax(3)=subplot(2,3,2);plot(t,AccY)
xlabel('Y axis')
ylabel('Accelerations (G)')
ax(4)=subplot(2,3,5);plot(t,GyrY,'r')
% title('Gyroscope')
xlabel('Y axis')
ylabel('Angular Velocity (rad/s)')
ax(5)=subplot(2,3,3);plot(t,AccZ)
xlabel('Z axis')
ylabel('Accelerations (G)')
ax(6)=subplot(2,3,6); plot(t,GyrZ,'r')
% title('Gyroscope')
xlabel('Z axis')
ylabel('Angular Velocity(rad/s)')
linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');

end
if gr
samples=1:numel(GyrY);
tt=samples/Fs;
% Disp_Y=cumtrapz(tt,GyrY(fs*10:fs*20));
% Disp_X=cumtrapz(tt,GyrX(fs*10:fs*20));
% Disp_Z=cumtrapz(tt,GyrZ(fs*10:fs*20));
% VelX=cumtrapz(tt,AccX(fs*10:fs*20));
% VelY=cumtrapz(tt,AccY(fs*10:fs*20));
% VelZ=cumtrapz(tt,AccZ(fs*10:fs*20));

Disp_Y=cumtrapz(tt,GyrY);
Disp_X=cumtrapz(tt,GyrX);
Disp_Z=cumtrapz(tt,GyrZ);
VelX=cumtrapz(tt,AccX);
VelY=cumtrapz(tt,AccY);
VelZ=cumtrapz(tt,AccZ);

opol = 15;
[p,s,mu] = polyfit(tt,Disp_Y',opol);
f_y = polyval(p,tt,[],mu);

dt_dispY = Disp_Y' - f_y;
% figure
% plot(dt_dispY)

opol = 15;
[p,s,mu] = polyfit(tt,Disp_X',opol);
f_x = polyval(p,tt,[],mu);

dt_dispX = Disp_X' - f_x;

opol = 15;
[p,s,mu] = polyfit(tt,Disp_Z',opol);
f_z = polyval(p,tt,[],mu);

dt_dispZ = Disp_Z' - f_z;
% figure
% plot(dt_dispZ)


figure
ax(1)=subplot(2,3,1);plot(tt,VelX)
title('Accelerometer')
xlabel('X axis')
ylabel('Velocity (m/s)')
ax(2)=subplot(2,3,4);plot(tt,dt_dispX,'r')
title('Gyroscope')
xlabel('X axis')
ylabel('Angular Displacement (rad)')

ax(3)=subplot(2,3,2);plot(tt,VelY)
xlabel('Y axis')
ylabel('Velocity (m/s)')
ax(4)=subplot(2,3,5);plot(tt,dt_dispY,'r')
% title('Gyroscope')
xlabel('Y axis')
ylabel('Angular Angular Displacement (rad)')
ax(5)=subplot(2,3,3);plot(tt,VelZ)
xlabel('Z axis')
ylabel('Velocity (m/s)')
ax(6)=subplot(2,3,6); plot(tt,dt_dispZ,'r')
% title('Gyroscope')
xlabel('Z axis')
ylabel('Angular Angular Displacement (rad)')
linkaxes([ax(6) ax(5) ax(4) ax(3) ax(2) ax(1)],'x');
end

total_disp_mag=sqrt(dt_dispX.^2+dt_dispY.^2+dt_dispZ.^2);
total_velocity_mag=sqrt(VelX.^2+VelY.^2+VelZ.^2);

total_rot=(GyrX.^2+GyrY.^2+GyrZ.^2);
total_trans=(AccX.^2+AccY.^2+AccZ.^2);

tot_rot_env= envelope( total_rot );
tot_trans_env= envelope( total_trans );

figure
ax(1)=subplot(2,1,1);plot(tt,total_disp_mag);
% hold on
% plot(t,tot_trans_env,'r')
title('Total RotationalDisplacement Magnitude')
xlabel('time (sec)')
ylabel('Arbitrary Units')
ax(2)=subplot(2,1,2);plot(tt,total_velocity_mag);
% hold on
% plot(t,tot_rot_env,'r');
title('Total Translational Velocity Magnitude')
xlabel('time (sec)')
ylabel('Arbitrary Units')
linkaxes([ax(2) ax(1)],'x');

figure
ax(1)=subplot(2,1,1);plot(t,total_trans);
hold on
plot(t,tot_trans_env,'r')
title('Total Translational Energy')
xlabel('time (sec)')
ylabel('Arbitrary Units')
ax(2)=subplot(2,1,2);plot(t,total_rot);
hold on
plot(t,tot_rot_env,'r');
title('Total Rotational Energy')
xlabel('time (sec)')
ylabel('Arbitrary Units')
linkaxes([ax(2) ax(1)],'x');

figure
ax(1)=subplot(3,1,1);plot(t,AccZ);
title('SCG Z')
xlabel('time (sec)')
ylabel('Acceleration (mG)')
ax(2)=subplot(3,1,2);plot(t,GyrY);
title('GCG Y')
xlabel('time (sec)')
ylabel('Angular Velocity (dps)')
ax(3)=subplot(3,1,3);plot(tt,dt_dispY,'r');
xlabel('Y axis')
ylabel('Angular Angular Displacement (rad)')
linkaxes([ ax(3) ax(2) ax(1)],'x');

end

