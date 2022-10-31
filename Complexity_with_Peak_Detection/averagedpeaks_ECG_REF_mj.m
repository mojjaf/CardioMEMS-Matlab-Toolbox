function [all_traces,avpeak,FS_new] = averagedpeaks_ECG_REF_mj(datain,fs, ploc)
% [avpeak,ploc] = averagedpeaks(signalin,fs,fignum)

if nargin<2
    fs = 200;
end
%UUOLLI
%   figure(fignum)
%UUOLLI
strech=1;
if strech
outputlen = 1000;

olli_avvi_median(1,1:outputlen)=0;
FS_new=[];
%ploc = findheartbeats(datain,fs);

num = numel(ploc)-1;
avpeak = zeros(1,outputlen);

else
    if fs<400
  outputlen = 300;
    else
          outputlen = 500;
    end
num = numel(ploc)-1; 
olli_avvi_median(num,1:outputlen)=0;
FS_new=[];
%ploc = findheartbeats(datain,fs);

 
    avpeak = zeros(num,outputlen);
end
% figure(1)


% olli_avvi_median={};
for ii = 1 : num
    
    
    if strech
        tmp = datain(ploc(ii):ploc(ii+1)-20);
        len = numel(tmp);
      
        %%%%% Mojtaba edited
        dx=1/fs;
        dy=dx*(len/outputlen);
        tmp2 = interpft(tmp,outputlen);
        
        fs_nw=round(1/dy);
%         
%         newFs = round(1/dy);
%         desiredFs = 1000;
% 
%         [p,q] = rat(desiredFs / newFs)
%         y = resample(tmp2,p,q);
        %%%%%
        
        avpeak = [avpeak ; tmp2(:)'];
        
        %plot(tmp2(:)'); hold on;
        FS_new(ii)=fs_nw;
        olli_avvi_median(ii,1:outputlen)=tmp2(:)';
        
    else
        tmp2 = datain(ploc(ii):ploc(ii+1));
%         avpeak = avpeak + tmp2(:)';
%         avpeak(ii,1:numel(tmp2)) = tmp2;
        FS_new(ii)=fs;
        olli_avvi_median(ii,1:numel(tmp2))=tmp2(:)';
%         olli_avvi_median{ii}=tmp2(:)';
    end
    
end


% if nargin>2
% %    figure(fignum)
% %    plot(avpeak)
%
% %plot(median(olli_avvi_median,1),'Color','b','LineWidth',3);
%
% % hold off;
% end
all_traces=olli_avvi_median;
avpeak=median(olli_avvi_median,1);
% avpeak=(olli_avvi_median);
