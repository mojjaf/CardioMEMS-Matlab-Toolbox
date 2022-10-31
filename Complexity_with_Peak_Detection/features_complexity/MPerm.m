function MPE = MPerm(X,m,t,Scale)
%  Calculate the Multiscale Permutation Entropy (MPE)
%  Input:   X: time series;
%           m: order of permuation entropy
%           t: delay time of permuation entropy, 
%           Scale: the scale factor
% Output: 
%           MPE: multiscale permuation entropy
%Ref: G Ouyang, J Li, X Liu, X Li, Dynamic Characteristics of Absence EEG Recordings with Multiscale Permutation %     %                             Entropy Analysis, Epilepsy Research, doi: 10.1016/j.eplepsyres.2012.11.003
%     G Ouyang, C Dang, X Li, Complexity Analysis of EEG Data with Multiscale Permutation Entropy, Advances in %       %                      Cognitive Neurodynamics (II), 2011, pp 741-745 
MPE=[];
for j=1:Scale
    Xs = Multi(X,j);
    PE = pec(Xs,m,t);
    MPE=[MPE PE];
end
function M_Data = Multi(Data,S)
%  generate the consecutive coarse-grained time series
%  Input:   Data: time series;
%           S: the scale factor
% Output: 
%           M_Data: the coarse-grained time series at the scale factor S
L = length(Data);
J = fix(L/S);
for i=1:J
M_Data(i) = mean(Data((i-1)*S+1:i*S));
end
