function [AO_amp_raw,AO_i_raw,heartrate_avg,t]=Modified_PanTompkins(data,fs,flag,motionArtifact,gr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%% Author : Mojtaba Jafari Tadi
% Technology Research Center, University of Turku
% email : mojjaf@utu.fi
% 

% Any direct or indirect use of this code should be referenced 
% Copyright JUNE 2015
%%
if ~isstruct(data)
  error('data must be a structure based input');
end
Fs=fs;

if nargin < 4
    gr = 1;   % on default the function always plots
end
% if nargin<5
%     coef=[1 1];
% end
% if isfield(data,'EKG')
%   ECG=double(data.EKG);
%     
%  elseif  isfield(data,'ECG2')
% 
%  ECG=double(data.ECG2);
% 
% elseif isfield(data,'ECG')
% ECG=double(data.ECG);
% else
%     ECG=double(data.accX);
% end
% ECG=medfilt1(ECG,3);

gyrY=double(data.gyroY);
gyrX=double(data.gyroX);
gyrZ=double(data.gyroZ);
accX=double(data.accX);
accY=double(data.accY);
accZ=double(data.accZ);
% ECG=fft_filter(ECG,800,5,40);
% ECG = detrend((ECG-mean(ECG))/std(ECG));

if motionArtifact

indvec=validseismodata(accZ, Fs);



readingJitter = Fs*2;
start=indvec(readingJitter);
stop=indvec(end-readingJitter);
samples=1:numel(gyrY(start:stop));
t=samples/Fs;

GyrY=fft_filter(gyrY(start:stop),fs,1,20,0);
GyrX=fft_filter(gyrX(start:stop),fs,1,20,0);
GyrZ=fft_filter(gyrZ(start:stop),fs,1,20,0);

AccZ=fft_filter(accZ(start:stop),fs,4,40,0);
AccY=fft_filter(accY(start:stop),fs,4,40,0);
AccX=fft_filter(accX(start:stop),fs,4,40,0);

% ECG=ECG(start:stop);
% sig_acc_z=fft_filter(AccZ,800,4,40);
% sig_gyr_y=fft_filter(GyrY,800,1,20);
SCG_Z=accZ(start:stop);
GCG_Y=gyrY(start:stop);
else
samples=1:numel(gyrY);
t=samples/Fs;
% 
% GyrY=fft_filter(gyrY,fs,1,20);
% GyrX=fft_filter(gyrX,fs,1,20);
% GyrZ=fft_filter(gyrZ,fs,1,20);
% 
% AccZ=fft_filter(accZ,fs,4,40);
% AccY=fft_filter(accY,fs,4,40);
% AccX=fft_filter(accX,fs,4,40); 
data_maf=struct('accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyrX,'gyroY',gyrY,'gyroZ',gyrZ);
  data_ssa=ssa_Analysis( data_maf,1,12 );
SCG_Z=data_ssa.accZ;
GCG_Y=data_ssa.gyroY;
end 

sig_acc_z=SCG_Z;
sig_gyr_y=GCG_Y;


% ECG=filter([1/15, 1/15, 1/15],1,ECG);
%% Initialize
AO_c =[]; %amplitude of R
AO_i =[]; %index
SIG_LEV = 0; 
nois_c =[];
nois_i =[];
delay = 0;
skip = 0; % becomes one when a T wave is detected
not_nois = 0; % it is not noise when not_nois = 1
selected_RR =[]; % Selected RR intervals
m_selected_AOAO = 0;
mean_AO = 0;
AO_i_raw =[];
AO_amp_raw=[];
ser_back = 0; 
test_m = 0;
SIGL_buf = [];
NOISL_buf = [];
THRS_buf = [];
SIGL_buf1 = [];
NOISL_buf1 = [];
THRS_buf1 = [];

%% 
%% Filtering  motion signals 


if flag==1

% sig_single=sig_acc_z(fs:end,1);
sig_single=sig_acc_z;

[ tmptrigy ] = HF_sphone(SCG_Z,flag );
%   tmptrigy= tmptrigy(fs:end,1);

if gr
    figure
ax(1)=subplot(211);plot(sig_single);axis tight;title('Original SCG Signal');
ax(2)=subplot(212);plot(tmptrigy);axis tight;title('Original SCG Signal');

end

elseif flag==2
%    sig_single= sig_gyr_y(fs:end,1);
       sig_single= diff(sig_gyr_y);

   [ tmptrigy ] = HF_sphone(GCG_Y,flag );
  tmptrigy= tmptrigy(fs:end,1);
   if gr
       figure
   ax(1)=subplot(211);plot(sig_single);axis tight;title('Original GCG signal');
   ax(2)=subplot(212);plot(tmptrigy);axis tight;title('Original GCG Signal');

   end
end
%     
% if flag==1
% 
% sig_single=sig_acc_z(fs:end,1);
% [ tmptrigy1 ] = HF_sphone(SCG_Z,flag );
% [ tmptrigy2 ] = HF_sphone(AccX,flag );
% tmptrigy= tmptrigy1+tmptrigy2;
%   tmptrigy= tmptrigy(fs:end,1);
% %  ECG=ECG(fs:end,1);
% if gr
% ax(2)=subplot(312);plot(sig_single);axis tight;title('Original SCG Signal');
% ax(3)=subplot(313);plot(tmptrigy);axis tight;title('Original SCG Signal');
% 
% end
% 
% elseif flag==2
%    sig_single= sig_gyr_y(fs:end,1);
% %     ECG=ECG(fs:end,1);
%    [ tmptrigy1 ] = HF_sphone(GCG_Y,flag );
%    [ tmptrigy2 ] = HF_sphone(GyrX,flag );
%    tmptrigy= tmptrigy1+tmptrigy2;
%   tmptrigy= tmptrigy(fs:end,1);
%    if gr
%    ax(2)=subplot(312);plot(sig_single);axis tight;title('Original SCG signal');
%    ax(3)=subplot(313);plot(tmptrigy);axis tight;title('Original SCG Signal');
% 
%    end
% end








%% Navigator signal for Peak Detection, i.e. Filtered version of SCG and GCG 
% Note : a minimum distance of 350 samples is considered between each heart beat wave
% since in physiological point of view no RR wave can occur in less than
% 200-300 msec distance
sig_comb_m=tmptrigy;
% [~,locs_comb]=findpeaks(sig_comb,'MinPeakHeight',(max(sig_comb)*0.3), 'MinPeakDistance',300);
[pks,locs]=findpeaks(sig_comb_m,'MinPeakHeight',(max(sig_comb_m(1:fs*10))*0.35), 'MinPeakDistance',fs*0.5);



%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(sig_comb_m(1:2*fs))*1/3; % 0.25 of the max amplitude
THR_NOISE = mean(sig_comb_m(1:2*fs))*1/2; % 0.5 of the mean signal is considered to be noise**2
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;


%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
THR_SIG1 = max(sig_single(1:2*fs))*1/3; % 0.25 of the max amplitude
THR_NOISE1 = mean(sig_single(1:2*fs))*1/2; %
SIG_LEV1 = THR_SIG1; % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1; % Noise level in Bandpassed filter
%% Thresholding and online desicion rule

for i = 1 : length(pks)
    
    %% locate the corresponding peak in the filtered signal
    if locs(i)-round(0.70*fs)>= 1 && locs(i)<= length(sig_single)
        [y_i x_i] = max(sig_single(locs(i)-round(0.7*fs):locs(i)));
    else
        if i == 1
            [y_i x_i] = max(sig_single(1:locs(i)));
            ser_back = 1;
        elseif locs(i)>= length(sig_single)
            [y_i x_i] = max(sig_single(locs(i)-round(0.70*fs):end));
        end
        
    end
    
    
    %% update the heart_rate (Two heart rate means one the moste recent and the other selected)
    if length(AO_c) >= 9
        
        diffRR = diff(AO_i(end-8:end)); %calculate AO-AO interval
        mean_AO = mean(diffRR); % calculate the mean of 8 previous AO waves interval
        comp =AO_i(end)-AO_i(end-1); %latest AO-AO
        if comp <= 0.92*mean_AO || comp >= 1.16*mean_AO
            % lower down thresholds to detect better in MVI
            THR_SIG = 0.5*(THR_SIG);
            %THR_NOISE = 0.5*(THR_SIG);
            % lower down thresholds to detect better in Bandpass filtered
            THR_SIG1 = 0.5*(THR_SIG1);
            %THR_NOISE1 = 0.5*(THR_SIG1);
            
        else
            m_selected_AOAO = mean_AO; %the latest regular beats mean
        end
        
    end
    
    %% calculate the mean of the last 8 AO waves to make sure that AO is not
    % missing(If no AO detected , trigger a search back) 1.66*mean
    
    if m_selected_AOAO
        test_m = m_selected_AOAO; %if the regular AO-AO availabe use it
    elseif mean_AO && m_selected_AOAO == 0
        test_m = mean_AO;
    else
        test_m = 0;
    end
    
    if test_m
        if (locs(i) - AO_i(end)) >= round(1.66*test_m)% it shows a AO is missed
            [pks_temp,locs_temp] = max(sig_comb_m(AO_i(end)+ round(0.250*fs):locs(i)-round(0.250*fs))); % search back and locate the max in this interval
            locs_temp = AO_i(end)+ round(0.250*fs) + locs_temp -1; %location
            
            if pks_temp > THR_NOISE
                AO_c = [AO_c pks_temp];
                AO_i = [AO_i locs_temp];
                
                % find the location in filtered sig
                if locs_temp <= length(sig_single)
                    [y_i_t x_i_t] = max(sig_single(locs_temp-round(0.70*fs):locs_temp));
                else
                    [y_i_t x_i_t] = max(sig_single(locs_temp-round(0.70*fs):end));
                end
                % take care of bandpass signal threshold
                if y_i_t > THR_NOISE1
                    
                    AO_i_raw = [AO_i_raw locs_temp-round(0.30*fs)+ (x_i_t - 1)];% save index of bandpass
                    AO_amp_raw =[AO_amp_raw y_i_t]; %save amplitude of bandpass
                    SIG_LEV1 = 0.25*y_i_t + 0.75*SIG_LEV1; %when found with the second thres
                end
                
                not_nois = 1;
                SIG_LEV = 0.25*pks_temp + 0.75*SIG_LEV ;  %when found with the second threshold
            end
            
        else
            not_nois = 0;
            
        end
    end
    
    
    
    
    %%  find noise and QRS peaks
    if pks(i) >= THR_SIG
        
        % if a AO candidate occurs within 360ms of the previous QRS
        % ,the algorithm determines if its AC wave or AO
        if length(AO_c) >= 3
            if (locs(i)-AO_i(end)) <= round(0.3600*fs)
                Slope1 = mean(diff(sig_comb_m(locs(i)-round(0.075*fs):locs(i)))); %mean slope of the waveform at that position
                Slope2 = mean(diff(sig_comb_m(AO_i(end)-round(0.075*fs):AO_i(end)))); %mean slope of previous AO wave
                if abs(Slope1) <= abs(0.5*(Slope2)) || (locs(i)-AO_i(end)) <= round(0.4*test_m)  % slope less then 0.5 of previous AO
                    nois_c = [nois_c pks(i)];
                    nois_i = [nois_i locs(i)];
                    skip = 1; % T wave identification
                    % adjust noise level in both filtered and
                    % MVI
                    NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
                    NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
                else
                    skip = 0;
                end
                
            end
        end
        
        if skip == 0  % skip is 1 when a AC wave is detected
            AO_c = [AO_c pks(i)];
            AO_i = [AO_i locs(i)];
            
            % bandpass filter check threshold
            if y_i >= THR_SIG1
                if ser_back
                    AO_i_raw = [AO_i_raw x_i];  % save index of bandpass
                else
                    AO_i_raw = [AO_i_raw locs(i)-round(0.70*fs)+ (x_i - 1)];% save index of bandpass
                end
                AO_amp_raw =[AO_amp_raw y_i];% save amplitude of bandpass
                SIG_LEV1 = 0.125*y_i + 0.875*SIG_LEV1;% adjust threshold for bandpass filtered sig
            end
            
            % adjust Signal level
            SIG_LEV = 0.125*pks(i) + 0.875*SIG_LEV ;
        end
        
        
    elseif THR_NOISE <= pks(i) && pks(i)<THR_SIG
        
        %adjust Noise level in filtered sig
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        %adjust Noise level in MVI
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
        
        
        
    elseif pks(i) < THR_NOISE
        nois_c = [nois_c pks(i)];
        nois_i = [nois_i locs(i)];
        
        % noise level in filtered signal
        NOISE_LEV1 = 0.125*y_i + 0.875*NOISE_LEV1;
        %end
        
        %adjust Noise level in MVI
        NOISE_LEV = 0.125*pks(i) + 0.875*NOISE_LEV;
        
        
    end
    
    
    
    
    
    %% adjust the threshold with SNR
    if NOISE_LEV ~= 0 || SIG_LEV ~= 0
        THR_SIG = NOISE_LEV + 0.25*(abs(SIG_LEV - NOISE_LEV));
        THR_NOISE = 0.5*(THR_SIG);
    end
    
    % adjust the threshold with SNR for bandpassed signal
    if NOISE_LEV1 ~= 0 || SIG_LEV1 ~= 0
        THR_SIG1 = NOISE_LEV1 + 0.25*(abs(SIG_LEV1 - NOISE_LEV1));
        THR_NOISE1 = 0.5*(THR_SIG1);
    end
    
    
    % take a track of thresholds of smoothed signal
    SIGL_buf = [SIGL_buf SIG_LEV];
    NOISL_buf = [NOISL_buf NOISE_LEV];
    THRS_buf = [THRS_buf THR_SIG];
    
    % take a track of thresholds of filtered signal
    SIGL_buf1 = [SIGL_buf1 SIG_LEV1];
    NOISL_buf1 = [NOISL_buf1 NOISE_LEV1];
    THRS_buf1 = [THRS_buf1 THR_SIG1];
    
    
    
    
    skip = 0; %reset parameters
    not_nois = 0; %reset parameters
    ser_back = 0;  %reset bandpass param
end

if gr
    hold on,scatter(AO_i,AO_c,'m');
    hold on,plot(locs,NOISL_buf,'--k','LineWidth',2);
    hold on,plot(locs,SIGL_buf,'--r','LineWidth',2);
    hold on,plot(locs,THRS_buf,'--g','LineWidth',2);
%     if ax(:)
%         linkaxes(ax,'x');
%         zoom on;
%     end
end

%%%%%%%%%%%%%%%%%%%
% [pk,ekg_locs] = findpeaks(ECG,'MinPeakHeight',(max(ECG)*0.3), 'MinPeakDistance',Fs*0.3);
% [pk,ekg_locs,delay]=pan_tompkin(ECG,Fs,1)
% ekg_locs=(ekg_locs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% kick out the outlier picks
% [a,b]=find(AO_amp_raw<(mean(AO_amp_raw)*0.35));
% % [c,d]=find(pk<(mean(pk)*0.35));
% 
% AO_amp_raw(b)=[];
% AO_i_raw(b)=[];
% 
% [aa,bb]=find(diff(AO_i_raw)<=300);
% AO_amp_raw(bb)=[];
% AO_i_raw(bb)=[];

%% overlay on the signals
if gr
    figure
    az(1)=subplot(311);plot(sig_single);
if flag==1
title('AO on Filtered Acc Signal');
elseif flag==2
title('AO on Filtered gyro Signal');
end
    
    axis tight;
    hold on,scatter(AO_i_raw,AO_amp_raw,'r');
    hold on,plot(locs,NOISL_buf1,'LineWidth',2,'Linestyle','--','color','k');
    hold on,plot(locs,SIGL_buf1,'LineWidth',2,'Linestyle','-.','color','g');
    hold on,plot(locs,THRS_buf1,'LineWidth',2,'Linestyle','-.','color','c');
    az(2)=subplot(312);plot(sig_comb_m);title('AO on MVI signal and Noise level(black),Signal Level (g) and Adaptive Threshold(cyan)');axis tight;
    hold on,scatter(AO_i,AO_c,'r');
    hold on,plot(locs,NOISL_buf,'LineWidth',2,'Linestyle','--','color','k');
    hold on,plot(locs,SIGL_buf,'LineWidth',2,'Linestyle','-.','color','g');
    hold on,plot(locs,THRS_buf,'LineWidth',2,'Linestyle','-.','color','c');
    az(3)=subplot(313);plot(sig_single-mean(sig_single));title('Pulse train of the found AO on MEMS signal');axis tight;
    line(repmat(AO_i_raw,[2 1]),repmat([min(sig_single-mean(sig_single))/2; max(sig_single-mean(sig_single))/2],size(AO_i_raw)),'LineWidth',2.5,'LineStyle','-.','Color','r');
    linkaxes(az,'x');
    zoom on;
end
% [CTI_values] = STIPD2( sig_single,AO_i_raw,fs,flag,gr)
[ heartrate_avg ] = heartrate( sig_single,AO_i_raw,fs)



end

