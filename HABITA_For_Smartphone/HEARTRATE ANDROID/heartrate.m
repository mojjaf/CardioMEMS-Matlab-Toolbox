function [ heartrate_median,HR_mean ] = heartrate( sig_in,locs,Fs)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
sig_len=numel(sig_in);

num=floor((sig_len-Fs)/Fs);
% rr_int=diff(locs);
rr_intval=Fs*10;
for ii=1:num-1    
    indvec =( (1:rr_intval) + (ii-1)*Fs);
    if max(indvec) < sig_len
        HR(ii)=60*(1/(mean(diff(locs(locs>min(indvec)& locs<max(indvec))))/Fs));
    end
end
heartrate_median=median(HR);
HR_mean=mean(HR);
end

