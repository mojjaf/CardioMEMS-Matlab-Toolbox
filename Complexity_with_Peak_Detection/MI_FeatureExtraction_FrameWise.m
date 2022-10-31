function [Features_Median] = MI_FeatureExtraction_FrameWise(frame, fs)
%%

% This function computes basic complexity feature sequencies for an
% accelerometer or gyroscope signal frame, on a short-term basis.
%
% ARGUMENTS:
%  - signal frame:    the Acc/Gyro signal
%  - fs:        the sampling frequency
%
% RETURNS:
%  - Features: a [MxN] matrix, where M is the number of short-term windows and N is
%  the total number of features. Each row of the matrix
%  corresponds to a seperate feature sequence
%

% compute the total number of frames:
% i=numOfFrames 
% number of features to be computed:
numOfFeatures = 15;

% Ham = window(@hamming, windowLength);

[row, col]=size(frame);
   Features = zeros(row,numOfFeatures);  
%   frame  = signal(curPos:curPos+windowLength-1);
   
  for i=1:row
        % compute time-domain features:
%          Features(i,1) = feature_zcr(frame(i,:));
         Features(i,1) = feature_energy(frame(i,:));
         Features(i,2) = feature_energy_entropy(frame(i,:), 10);
         Features(i,3) = auc_energy_norm(frame(i,:));
         Features(i,4) =  ApEn( 1, 0.2 * std(frame(i,:)), frame(i,:), 1 );
         [Features(i,5),Features(i,6)] = HjorthParameters(frame(i,:)');
         Features(i,7) = ShanEn( frame(i,:) );
         Features(i,8) = MAS(frame(i,:)); 
         Features(i,9) = tpr(frame(i,:));
         Features(i,10)=Tsallis_entro(frame(i,:)/max(frame(i,:)),2);
         Features(i,11)=renyi_entro(frame(i,:),2);
%          
         Features(i,12) = Katz_FD(frame(i,:));
         Features(i,13) = Higuchi_FD(frame(i,:), 3);
         Features(i,14) = genhurst(frame(i,:));
%         [Features(i,17),Features(i,18)] = DFA_fun(frame(i,:),20,1);
         %Features(i,15) = MPerm(frame(i,:),3,100,1); %%multiscale permutation entropy
         Features(i,15)=iqr(frame(i,:)/max(frame(i,:)));
         Features(i,16)=max(frame(i,:))-min(frame(i,:));
        
  end
  Features_Median=median(Features);
end

