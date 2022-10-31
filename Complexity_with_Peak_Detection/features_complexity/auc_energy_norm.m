function [AUC] = auc_energy_norm(frame)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
AUC=trapz(frame.^2)/max(frame.^2);
end

