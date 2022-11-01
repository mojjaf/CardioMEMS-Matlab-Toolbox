function [ rricomb,rriekg,HR_6x,HR_ekg ] = sixaxisxcorr_sp( data,gr,Fs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

if ~isstruct(data)
  error('data must be a structure based input');
end
Fs=800;

if nargin < 4
    gr = 1;   % on default the function always plots
end
if isfield(data,'EKG')
    ECG=double(data.EKG);
  elseif  isfield(data,'ECG2')

  ECG=double(data.ECG2);

elseif isfield(data,'ECG')
ECG=double(data.ECG);
else
    ECG=double(data.accX);
end


GyrY=double(data.gyroY)*(2000/32767)/256;
GyrX=double(data.gyroX)*(2000/32767)/256;
GyrZ=double(data.gyroZ)*(2000/32767)/256;
AccX=double(data.accX)*(2000/32760);
AccY=double(data.accY)*(2000/32760);
AccZ=double(data.accZ)*(2000/32760);
ECG=fft_filter(ECG,Fs,5,40);
ECG = detrend((ECG-mean(ECG))/std(ECG));

indvec=validseismodata(GyrY, Fs);
start=indvec(100);
stop=indvec(end-100);
samples=1:numel(GyrY(start:stop));
t=samples/Fs;


GyrY=GyrY(start:stop);
GyrX=GyrX(start:stop);
GyrZ=GyrZ(start:stop);

AccZ=AccZ(start:stop);
AccY=AccY(start:stop);
AccX=AccX(start:stop);

ECG=ECG(start:stop);
fieldnames = {'accZ','accY','accX','gyroY','gyroZ','gyroX','ekg'};
fieldnum = 7;
s = struct(fieldnames{1},AccZ,fieldnames{2},AccY,fieldnames{3},AccX,fieldnames{4},...
    GyrY,fieldnames{5},GyrZ,fieldnames{6},GyrX,fieldnames{7},ECG)




rriekgvec = [];
rricombvec = [];
grriall     = [  ];
epeakrriall = [ ];
ekgacall    = [   ];
esn = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             BAND PASS FILTER
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nfir = 16;
nfir2 = 32;
bhp = -ones(1,nfir)/nfir;
bhp(1) = bhp(1)+1;
blp = triang(nfir2);


nfir = 16;
nfir2 = 32;
bhp = -ones(1,nfir)/nfir;
bhp(1) = bhp(1)+1;
blp = triang(nfir2);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data.measurement_name
% data.sensors
% data.test_subject
% % indvec = 1.5e5:4.5e5;
% afib = data.AFib;
for ss = 1 : 6;
    fieldname = fieldnames{ss};
    tmpdat = s.(fieldname);
    tmpdat = filter(bhp,1,tmpdat);
    tmpdat = filter(blp,1,tmpdat);
    tmpdat = [tmpdat(16:end); zeros(15,1)];
    rrivec{ss} = xcorrsykelyhyt(tmpdat,Fs);
    if ss == 1
        rricomb = rrivec{ss};
    else
        rricomb = combinerrint(rricomb,rrivec{ss});
    end
end


rricomb = rricomb(2:end) ;
ekg = fft_filter(double(s.ekg),Fs,1,45);
ekg(1:400)=0;
rriekg = ekgrrint(ekg);

HR_6x=60*(1/(mean(rricomb)/Fs));
HR_ekg=60*(1/(mean(rriekg)/Fs));

if gr
figure
plot(rriekg)
hold on
plot(rricomb,'r--')
hold off
legend({'RRI ekg','RRI 6axAC'})
title('ekgrri vs acrri')


figure
plen = min(length(rriekg),length(rricomb));
pindvec  =1:plen;
plot(rriekg(pindvec)-rricomb(pindvec))
title rrierror
end
end




