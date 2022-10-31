function [ features,features_avg,features_fusion,features_acc,features_gyro] = Feature_Extraction_ADXL( data,fs,win, step,lab, id )
%UNTITLED17 Summary of this function goes here
%   Detailed explanation goes here


 Features_Ax_hf = FeatureExtraction_adxl(data.accX,fs , win, step );
  Features_Ay_hf = FeatureExtraction_adxl(data.accY,fs , win, step );
   Features_Az_hf = FeatureExtraction_adxl(data.accZ,fs , win, step );
    Features_Ax_lf = FeatureExtraction_adxl(data.accXLF,fs , win, step );
     Features_Ay_lf = FeatureExtraction_adxl(data.accYLF,fs , win, step);
      Features_Az_lf = FeatureExtraction_adxl(data.accZLF,fs , win, step );
 if lab
labels=ones(size(Features_Ax_lf,2),1);
% labels_x=ones(6,1);
labels_x=ones(size(Features_Ax_lf,2),1);
 else
   labels=zeros(size(Features_Ax_lf,2),1);
   labels_x=zeros(size(Features_Ax_lf,2),1);
 end
idx=zeros(size(labels));
idx(:,:)=id;

ID_X=zeros(size(labels_x));
ID_X(:,:)=id;
features=[Features_Ax_hf',Features_Ay_hf',Features_Az_hf',Features_Ax_lf',Features_Ay_lf',Features_Az_lf',labels,idx];
% 
% F_ax=median(Features_Ax');
% F_ay=median(Features_Ay');
% F_az=median(Features_Az');
% F_gx=median(Features_Gx');
% F_gy=median(Features_Gy');
% F_gz=median(Features_Gz');
Features_3D(1,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Ax_hf;
Features_3D(2,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Ay_hf;
Features_3D(3,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Az_hf;
Features_3D(4,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Ax_lf;
Features_3D(5,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Ay_lf;
Features_3D(6,1:size(Features_Ax_hf,1),1:size(Features_Ax_hf,2))=Features_Az_lf;
features_avg=median(Features_3D,1);
% features_avg=[F_gx,F_gy,F_gz,F_ax,F_ay,F_az,labels_x,ID_X];
feat_temp = reshape(features_avg,size(features_avg,1)*size(features_avg,2),size(features_avg,3));
features_avg=[feat_temp',labels_x,ID_X];

features_fusion=[Features_Ay_hf',Features_Az_lf',labels,idx];

features_acc=[Features_Ax_lf',Features_Ay_lf',Features_Az_lf',labels,idx];
features_gyro=[Features_Ax_hf',Features_Ay_hf',Features_Az_hf',labels,idx ];
end

