function [ sig_out ] = hilbert_transform( sig_in,Fs )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

x=sig_in;

signal_ht=sqrt(x.^2+abs(hilbert(x)).^2);
sig_ht_filt=fft_filter(signal_ht,Fs,0.5,2.5,0);
% sig_out=filter([1/3, 1/3, 1/3],1,sig_ht_filt);
% sig_out=filter([1/5, 1/5, 1/5],1,sig_ht_filt); Changed 170715
sig_out=conv(sig_ht_filt,ones(5,1)/5,'same');
end

