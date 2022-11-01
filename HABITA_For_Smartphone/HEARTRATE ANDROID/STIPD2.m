function [CTI_values] = STIPD2( sig_single,AO_i_raw,fs,flag,gr)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Fs=fs;
MC_location=[];
MI=[];
% MA=ones(1,numel(AO_i_raw));
MA=[];
V1=[];
% V2=zeros(1,numel(AO_i_raw));
% V2=ones(1,numel(AO_i_raw));
V2=[];
MO=[];
AC=[];
R_MC=[];% electromechanical systole (QS2)
R_AO=[];% pre-ejection period (PEP)
LVET=[]; %left ventricular ejection time (LVET)
MC_AO=[]; %Isovolumetric contraction time (IVCT)
dias=[];

for i=2:(numel(AO_i_raw)-2)
        peak_sample_n = AO_i_raw(i); %in samples
        start = peak_sample_n-0.125*fs;
        stop = peak_sample_n-0.019*fs ;
        sample1 = sig_single(start:stop);
        
     for j=1:numel(sample1)
         if sample1(j)==min(sig_single(start:stop))
           MC_location(i)=start+j-1;
         end
     end
       
end
for i=2:(numel(AO_i_raw)-2)
%         peak_sample_n = locs_AOwave(i); %in samples
        start = MC_location(i);
        stop = AO_i_raw(i) ;
        sample2 = sig_single(start:stop);
        
     for j=1:numel(sample2)
         if sample2(j)==min(sig_single(start:stop))
           MI(i)=start+j-1;
         end
     end
       
end
for i=2:(numel(AO_i_raw)-2)
        peak_sample_n = AO_i_raw(i); %in samples
        start = peak_sample_n-1;
        stop = peak_sample_n+0.044*fs ;
        sample3 = sig_single(start:stop);
        
     for j=1:numel(sample3)
         if sample3(j)==min(sig_single(start:stop))
           MA(i)=start+j-1;
         end
     end
       
end



for i=2:(numel(AO_i_raw)-2)
    peak_sample_n = AO_i_raw(i); %in samples
    start = peak_sample_n+0.25*fs ;
    if flag==1
        stop = start+0.25*fs;
    else
        stop = start+0.375*fs;
    end
    sample4 = sig_single(start:stop);
    
    for j=1:numel(sample4)
        if sample4(j)==min(sig_single(start:stop))
            V2(i)=start+j-1;
        end
    end
    
end

for i=2:(numel(V2)-2)
        peak_sample_n = V2(i); %in samples
        start = peak_sample_n+1 ;
        stop = peak_sample_n+0.0375*fs;
        sample5 = sig_single(start:stop);
        
     for j=1:numel(sample5)
         if sample5(j)==max(sig_single(start:stop))
           MO(i)=start+j-1;
         end
     end
       
end



for i=2:(numel(V2)-2)
        peak_sample_n = V2(i); %in samples
        start = peak_sample_n-0.0625*fs ;
        stop = peak_sample_n-0.006;
        sample6 = sig_single(start:stop);
        
     for j=1:numel(sample6)
         if sample6(j)==max(sig_single(start:stop))
           V1(i)=start+j-1;
         end
         
     end
       
end
% for i=1:(numel(AO_i_raw)-1)
%         peak_sample_n = V1(i); %in samples
%         start = peak_sample_n-45 ;
%         stop = peak_sample_n-10;
%         sample = sig_single(start:stop);
%         
%      for j=1:numel(sample)
%          if sample(j)==max(sig_single(start:stop))
%            AC(i)=start+j-1;
%          end
%          
%      end
%        
% end
V1=V1(2:end);
V2=V2(2:end);
AC=V2(2:end);
MC_location=MC_location(2:end);
MI=MI(2:end);
% MA=ones(1,numel(AO_i_raw));
MA=MA(2:end);

MO=MO(2:end);
AO_i_raw=AO_i_raw(2:end-1);
for m=1:length(AO_i_raw)-1
%     R_MC(m-1)=MC_location(m)-ekg_locs(m);
%     R_AO(m-1)=AO_i_raw(m)-ekg_locs(m);
% MC_AO(m-1)=AO_i_raw(m)-MC_location(m);
LVET(m)=V2(m)-AO_i_raw(m);
dias(m)=AO_i_raw(m+1)-V2(m);
end
% end

QS2=(mean(R_MC)/Fs)*1000;
PEP=(mean(R_AO)/Fs)*1000;
IVCT=(mean(MC_AO)/Fs)*1000;
LVET=(mean(LVET(3:end-3))/Fs)*1000;
LVRT=(mean(dias(3:end-3))/Fs)*1000;
% electromechanical systole (QS2),  and 
CTI_values=[LVET,LVRT];
if gr
figure
plot(sig_single);
hold on
plot(AO_i_raw,sig_single(AO_i_raw),'om');
hold on
% plot(MA,sig_single(MA),'ks','MarkerFaceColor','c');
% hold on
plot(V2,sig_single(V2),'oc');
hold on

if flag==1
plot(AC,sig_single(AC),'*r')
hold on
plot(MC_location,sig_single(MC_location),'*c');
hold on
plot( MI,sig_single( MI),'*g');
hold on
plot(MO,sig_single(MO),'*b')
hold on
end
% line(repmat(V2,[2 1]),repmat([min(sig_single-mean(sig_single))/2; max(sig_single-mean(sig_single))/2],size(V2)),'LineWidth',2.5,'LineStyle','-.','Color','r');
if flag==1
xlabel('Samples'); ylabel('Accelerations')
title('SYTOLIC TIME INTERVAL PEAK DETECTION- SCG')
legend('SCG Signal','AO','MC','MI','MA','V2','AC','MO');
elseif flag==2
  xlabel('Samples'); ylabel('Angular velocity')
title('SYTOLIC TIME INTERVAL PEAK DETECTION- GCG')
legend('GCG Signal','AO','AC');
end
line(repmat(V2,[2 1]),repmat([min(sig_single-mean(sig_single))/2; max(sig_single-mean(sig_single))/2],size(V2)),'LineWidth',2.5,'LineStyle','-.','Color','r');

end
n=numel(AO_i_raw)-2;
mat_MEMS=zeros(n,Fs*1.5);
Jittar=0.125*fs;
for i=2:n
        peak_sample_n = AO_i_raw(i); %in samples
        start = peak_sample_n-Jittar ;
        stop = AO_i_raw(i+1)-Jittar;
        sample = sig_single(start:stop);
        
%       j=numel(sample);
          mat_MEMS(i,1:numel(sample))=sample;
          
    
end
     if gr
         figure
plot((mat_MEMS'))  
hold on
plot(mean((mat_MEMS)),'LineWidth',4.5,'LineStyle','-','Color','r')  

     end

end

