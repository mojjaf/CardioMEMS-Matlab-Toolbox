function [total_vib_comlx,total_vib_mob] = complxity_and_mobility_segmented(vibtraces)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

n=size(vibtraces,1);
MC_temp=[];
MM_temp=[];
for i=1:n

[mobility(i),complexity(i)] = HjorthParameters(vibtraces(i,:)');
% MC_temp=[MC_temp;complexity];
% MM_temp(i)=[MM_temp;mobility];
end
total_vib_comlx=(mobility);
total_vib_mob=(complexity);

end 
