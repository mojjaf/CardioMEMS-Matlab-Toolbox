function [slope] = wave_slope(frame,fs)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
sample=1:numel(frame);
tt=sample/fs;
coefficients = polyfit(tt, frame, 1);
slope = coefficients(1);
end

