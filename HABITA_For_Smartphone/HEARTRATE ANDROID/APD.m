function [AO_amp_raw,AO_i_raw,pk,ekg_locs]=APD(data,fs,start,stop,flag,gr)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
% Complete implementation of Seismocardiogram and Gyrocardiogram peak detection algorithm

%% Inputs
% data : raw measured data including ECG, 3 axis Accelerometer(AccX, AccY and AccZ) and 3 Axes
% Gyroscope signals(GyrX,GyrY and GyrZ)
% fs : sampling frequency e.g. 800 Hz
% gr : flag to plot or not plot (set it 1 to have a plot or set it zero not
% to see any plots
% flag: 1 to have peak detection in Gyroscope signal(Y axis) and 2 to have peak
% detection in Accelerometer (Z axis)
%% Outputs
% AO_amp_raw : amplitude of AO waves amplitudes
% qrs_i_raw : index/ location of AO waves
% delay : number of samples which the signal is delayed due to the
% filtering
%% Method :

%% PreProcessing
% % 1) Signal is preprocessed, first the raw motion signals obtained from each axes (
% i.e.  X,Y and Z)of motion sensors were combined together by summing up the acceleration and rotational informa-
% tion  to  form  a  new  signal.   The  new  summed  signals  called  3D-SCG  and
% 3D-GCG  are  pre-processed  by  applying  a  moving  average  ?lter  as  well  as
% band-pass ?ltering to reduce noise,  baseline wander,  muscle noise and etc 
% (a combination of low pass and high pass filter 5-15 Hz)
% % to get rid of the baseline wander and muscle noise. 

% 2) The filtered signal
% is derivated using a derivating filter to high light the QRS complex.

% 3) Signal is squared.4)Signal is averaged with a moving window to get rid
% of noise (0.150 seconds length).

% 5) depending on the sampling frequency of your signal the filtering
% options are changed to best match the characteristics of your ecg signal

% 6) Unlike the other implementations in this implementation the desicion
% rule of the Pan tompkins is implemented completely.

%% Decision Rule 
% At this point in the algorithm, the preceding stages have produced a roughly pulse-shaped
% waveform at the output of the MWI . The determination as to whether this pulse
% corresponds to a QRS complex (as opposed to a high-sloped T-wave or a noise artefact) is
% performed with an adaptive thresholding operation and other decision
% rules outlined below;

% a) FIDUCIAL MARK - The waveform is first processed to produce a set of weighted unit
% samples at the location of the MWI maxima. This is done in order to localize the QRS
% complex to a single instant of time. The w[k] weighting is the maxima value.

% b) THRESHOLDING - When analyzing the amplitude of the MWI output, the algorithm uses
% two threshold values (THR_SIG and THR_NOISE, appropriately initialized during a brief
% 2 second training phase) that continuously adapt to changing ECG signal quality. The
% first pass through y[n] uses these thresholds to classify the each non-zero sample
% (CURRENTPEAK) as either signal or noise:
% If CURRENTPEAK > THR_SIG, that location is identified as a “QRS complex
% candidate” and the signal level (SIG_LEV) is updated:
% SIG _ LEV = 0.125 ×CURRENTPEAK + 0.875× SIG _ LEV

% If THR_NOISE < CURRENTPEAK < THR_SIG, then that location is identified as a
% “noise peak” and the noise level (NOISE_LEV) is updated:
% NOISE _ LEV = 0.125×CURRENTPEAK + 0.875× NOISE _ LEV
% Based on new estimates of the signal and noise levels (SIG_LEV and NOISE_LEV,
% respectively) at that point in the ECG, the thresholds are adjusted as follows:
% THR _ SIG = NOISE _ LEV + 0.25 × (SIG _ LEV ? NOISE _ LEV )
% THR _ NOISE = 0.5× (THR _ SIG)
% These adjustments lower the threshold gradually in signal segments that are deemed to
% be of poorer quality.


% c) SEARCHBACK FOR MISSED QRS COMPLEXES - In the thresholding step above, if
% CURRENTPEAK < THR_SIG, the peak is deemed not to have resulted from a QRS
% complex. If however, an unreasonably long period has expired without an abovethreshold
% peak, the algorithm will assume a QRS has been missed and perform a
% searchback. This limits the number of false negatives. The minimum time used to trigger
% a searchback is 1.66 times the current R peak to R peak time period (called the RR
% interval). This value has a physiological origin - the time value between adjacent
% heartbeats cannot change more quickly than this. The missed QRS complex is assumed
% to occur at the location of the highest peak in the interval that lies between THR_SIG and
% THR_NOISE. In this algorithm, two average RR intervals are stored,the first RR interval is 
% calculated as an average of the last eight QRS locations in order to adapt to changing heart 
% rate and the second RR interval mean is the mean 
% of the most regular RR intervals . The threshold is lowered if the heart rate is not regular 
% to improve detection.

% d) ELIMINATION OF MULTIPLE DETECTIONS WITHIN REFRACTORY PERIOD - It is
% impossible for a legitimate QRS complex to occur if it lies within 200ms after a previously
% detected one. This constraint is a physiological one – due to the refractory period during
% which ventricular depolarization cannot occur despite a stimulus[1]. As QRS complex
% candidates are generated, the algorithm eliminates such physically impossible events,
% thereby reducing false positives.

% e) T WAVE DISCRIMINATION - Finally, if a QRS candidate occurs after the 200ms
% refractory period but within 360ms of the previous QRS, the algorithm determines
% whether this is a genuine QRS complex of the next heartbeat or an abnormally prominent
% T wave. This decision is based on the mean slope of the waveform at that position. A slope of
% less than one half that of the previous QRS complex is consistent with the slower
% changing behaviour of a T wave – otherwise, it becomes a QRS detection.
% Extra concept : beside the points mentioned in the paper, this code also
% checks if the occured peak which is less than 360 msec latency has also a
% latency less than 0,5*mean_RR if yes this is counted as noise

% f) In the final stage , the output of R waves detected in smoothed signal is analyzed and double
% checked with the help of the output of the bandpass signal to improve
% detection and find the original index of the real R waves on the raw ecg
% signal



%% Author : Mojtaba Jafari Tadu
% Technology Research Center, University of Turku
% email : mojjaf@utu.fi
% 

% Any direct or indirect use of this code should be referenced 
% Copyright MAY 2015
%%
% if ~isstruct(data)
%   error('data must be a structure based input');
% end
Fs=fs;

if nargin < 4
    gr = 1;   % on default the function always plots
end
% if isnumeric(data.EKG)
%   ECG=double(data.EKG);
% elseif isnumeric(data.ECG)
% end
if isstruct(data)
GyrY=double(data.gyroY)*(2000/32767)/256;
GyrX=double(data.gyroX)*(2000/32767)/256;
GyrZ=double(data.gyroZ)*(2000/32767)/256;
AccX=double(data.accX)*(2000/32760);
AccY=double(data.accY)*(2000/32760);
AccZ=double(data.accZ)*(2000/32760);
samples=1:numel(GyrY(start:stop));
t=samples/Fs;
ECG=double(data.ECG);


GyrY=GyrY(start:stop);
GyrX=GyrX(start:stop);
GyrZ=GyrZ(start:stop);

AccZ=AccZ(start:stop);
AccY=AccY(start:stop);
AccX=AccX(start:stop);

ECG=ECG(start:stop);
sig_acc_z=fft_filter(AccZ,Fs,4,40);
sig_gyr_y=fft_filter(GyrY,Fs,1,20);
ECG=fft_filter(ECG,Fs,1,40);
ECG = (ECG-mean(ECG))/std(ECG);
else
AccX=data(:,1);
AccY=data(:,2);
AccZ=data(:,3);
GyrX=data(:,4);
GyrY=data(:,5);
GyrZ=data(:,6);
GyrY=GyrY(start:stop);
GyrX=GyrX(start:stop);
GyrZ=GyrZ(start:stop);
AccZ=AccZ(start:stop);
AccY=AccY(start:stop);
AccX=AccX(start:stop);
samples=1:numel(GyrY);
t=samples/Fs;
sig_acc_z=fft_filter(AccZ,Fs,4,40);
sig_gyr_y=fft_filter(GyrY,Fs,1,20);
end
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

%% Plot ECG 
if gr 
    figure
       if flag==1
         
        annotation('textbox', [0, 0.5, 0, 0], 'string', 'My Gyroscope')
       else 
      annotation('textbox', [0, 0.5, 0, 0], 'string', 'My Acceleromter')
     
    end
 if fs == 800 && exist('ECG','var')
    ax(1)=subplot(321);plot(ECG);axis tight;title('Raw ECG Signal');
 else
     
end    
%% Summing up motion signals
sig_acc=detrend((AccZ+AccX+AccY));
sig_gyr=(GyrY+GyrX+GyrZ);


if gr 
    if ~exist('ECG','var')
        ax(1)=subplot(3,2,1);plot(sig_acc);axis tight;title('Raw combined 3D SCG');
        ax(2)=subplot(322);plot(sig_gyr);axis tight;title('Raw combined 3D GCG');

    else
        ax(2)=subplot(322);plot(sig_acc);axis tight;title('Raw combined 3D SCG');
        ax(3)=subplot(323);plot(sig_gyr);axis tight;title('Raw combined 3D GCG');
    end
end




%% Hilbert Transform and Noise cancelation(Filtering) % Moving average and Bandpass Filters 
% (MV= 5 samples, cutt of freq for BP Filter 4-40 Hz for SCG and 1-15 Hz for GCG)

sig_acc_filt= Transform( sig_acc,Fs,1 );
sig_gyr_filt= Transform( sig_gyr,Fs ,2);
if gr
if ~exist('ECG','var')
    ax(3)=subplot(323);plot(sig_acc_filt);axis tight;title('Transformed 3D SCG');
    ax(4)=subplot(324);plot(sig_gyr_filt);axis tight;title('Transformed 3D GCG');
    
else
    ax(4)=subplot(324);plot(sig_acc_filt);axis tight;title('Transformed 3D SCG');
    ax(5)=subplot(325);plot(sig_gyr_filt);axis tight;title('Transformed 3D GCG');
end
end

%% Combining transfomed version of SCG and GCG by summing up 3D-SCG and 3D-GCG
acc_int=cumtrapz(t,sig_acc_filt);
gyr_int=cumtrapz(t,sig_gyr_filt);
% sig_comb=sig_acc_filt+sig_gyr_filt*10;
% [ acc_int ] = normalize( acc_int );
% [ gyr_int ] = normalize( sig_gyr_filt );


sig_comb_int=acc_int+gyr_int;

% sig_comb_int=cumtrapz(t,sig_comb);
sig_comb_int=filter(0.3,1,sig_comb_int);

sig_comb_int = conv(sig_comb_int ,ones(1 ,round(0.150*fs))/round(0.150*fs));
delay = delay + 15;
sig_comb_m=sig_comb_int;
if ~exist('ECG','var')
    ax(6)=subplot(3,2,[5 6]);plot(sig_comb_int);axis tight;title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
    axis tight;
else
    ax(6)=subplot(326);plot(sig_comb_int);axis tight;title('Averaged with 30 samples length,Black noise,Green Adaptive Threshold,RED Sig Level,Red circles QRS adaptive threshold');
    axis tight;
end

%% Navigator signal for Peak Detection, i.e. Integrated version of Combined 3D SCG and GCG 
% Note : a minimum distance of 80 samples is considered between each heart beat wave
% since in physiological point of view no RR wave can occur in less than
% 200-300 msec distance
wander=abs(max(sig_comb_int)-min(sig_comb_int))/2;
sig_comb_m=sig_comb_int+wander;
% [~,locs_comb]=findpeaks(sig_comb,'MinPeakHeight',(max(sig_comb)*0.3), 'MinPeakDistance',300);
[pks,locs]=findpeaks(sig_comb_m,'MinPeakHeight',(max(sig_comb_m)*0.25), 'MinPeakDistance',80);

if flag==1

sig_comb=sig_acc_z;
elseif flag==2
   sig_comb= sig_gyr_y;
end
    

%% initialize the training phase (2 seconds of the signal) to determine the THR_SIG and THR_NOISE
THR_SIG = max(sig_comb_m(1:2*fs))*1/3; % 0.25 of the max amplitude 
THR_NOISE = mean(sig_comb_m(1:2*fs))*1/2; % 0.5 of the mean signal is considered to be noise
SIG_LEV= THR_SIG;
NOISE_LEV = THR_NOISE;


%% Initialize bandpath filter threshold(2 seconds of the bandpass signal)
THR_SIG1 = max(sig_comb(1:2*fs))*1/3; % 0.25 of the max amplitude 
THR_NOISE1 = mean(sig_comb(1:2*fs))*1/2; %
SIG_LEV1 = THR_SIG1; % Signal level in Bandpassed filter
NOISE_LEV1 = THR_NOISE1; % Noise level in Bandpassed filter
%% Thresholding and online desicion rule

for i = 1 : length(pks)
    
   %% locate the corresponding peak in the filtered signal 
    if locs(i)-round(0.70*fs)>= 1 && locs(i)<= length(sig_comb)
          [y_i x_i] = max(sig_comb(locs(i)-round(0.7*fs):locs(i)));
       else
          if i == 1
            [y_i x_i] = max(sig_comb(1:locs(i)));
            ser_back = 1;
          elseif locs(i)>= length(sig_comb)
            [y_i x_i] = max(sig_comb(locs(i)-round(0.70*fs):end));
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
               if locs_temp <= length(sig_comb)
                [y_i_t x_i_t] = max(sig_comb(locs_temp-round(0.70*fs):locs_temp));
               else
                [y_i_t x_i_t] = max(sig_comb(locs_temp-round(0.70*fs):end));
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
if ax(:)
linkaxes(ax,'x');
zoom on;
end
end

%%%%%%%%%%%%%%%%%%%
if isstruct(data)
[pk,ekg_locs] = findpeaks(ECG,'MinPeakHeight',(max(ECG)*0.3), 'MinPeakDistance',Fs*0.5);
else
    pk=0;
    ekg_locs=0;
end




%% overlay on the signals
if gr 
    figure
    if flag==1
         
        annotation('textbox', [0, 0.5, 0, 0], 'string', 'My Gyroscope')
    else 
      annotation('textbox', [0, 0.5, 0, 0], 'string', 'My Acceleromter')
     
    end
    if exist('ECG','var')
        az(1)=subplot(411);plot(ECG)
        hold on
        plot(ekg_locs,ECG(ekg_locs),'g+');
        title('Reference ECG')
        axis tight;
        xlabel('Time (s)')
        ylabel('ECG (v)')
         az(2)=subplot(412);plot(sig_comb);title('AO on Filtered Signal');axis tight;
    else
         az(2)=subplot(4,1,[1 2]);plot(sig_comb);title('AO on Filtered Signal');axis tight;
    end
    
%     az(2)=subplot(412);plot(sig_comb);title('AO on Filtered Signal');axis tight;
    hold on,scatter(AO_i_raw,AO_amp_raw,'r');
    hold on,plot(locs,NOISL_buf1,'LineWidth',2,'Linestyle','--','color','k');
    hold on,plot(locs,SIGL_buf1,'LineWidth',2,'Linestyle','-.','color','g');
    hold on,plot(locs,THRS_buf1,'LineWidth',2,'Linestyle','-.','color','c');
    az(3)=subplot(413);plot(sig_comb_m);title('AO on MVI signal and Noise level(black),Signal Level (g) and Adaptive Threshold(cyan)');axis tight;
    hold on,scatter(AO_i,AO_c,'r');
    hold on,plot(locs,NOISL_buf,'LineWidth',2,'Linestyle','--','color','k');
    hold on,plot(locs,SIGL_buf,'LineWidth',2,'Linestyle','-.','color','g');
    hold on,plot(locs,THRS_buf,'LineWidth',2,'Linestyle','-.','color','c');
    az(4)=subplot(414);plot(sig_comb-mean(sig_comb));title('Pulse train of the found AO on MEMS signal');axis tight;
    line(repmat(AO_i_raw,[2 1]),repmat([min(sig_comb-mean(sig_comb))/2; max(sig_comb-mean(sig_comb))/2],size(AO_i_raw)),'LineWidth',2.5,'LineStyle','-.','Color','r');
    linkaxes(az,'x');
    zoom on;
end
end
 













