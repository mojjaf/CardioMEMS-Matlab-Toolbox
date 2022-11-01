function [ sig_out ] = Transform( sig_in,Fs,flag )
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
if nargin < 3
    Fs = 800;   % on default the function always plots
end
sig_mv=filter([1/3, 1/3, 1/3],1,sig_in);
if flag==1
sig_mv_fl=fft_filter(sig_mv,Fs,4,30);
elseif flag==2
    sig_mv_fl=fft_filter(sig_mv,Fs,1,15);
end
 x=sig_mv_fl;
signal_ht=sqrt(x.^2+abs(hilbert(x)).^2);
sig_ht_filt=fft_filter(signal_ht,Fs,0.5,2.5);
sig_out=filter([1/3, 1/3, 1/3],1,sig_ht_filt);

end

