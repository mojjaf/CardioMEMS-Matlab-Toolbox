%  close all;
 clear all;
Features_energy_complexity_Post=[];Features_energy_complexity_pre=[];
Features_Frame_wise_Pre=[];Features_sensorFusion_Pre=[];Features_medianFusion_Pre=[];Features_Accelerometer_Pre=[];Features_Gyroscope_Pre=[];
Features_Frame_wise_Post=[];Features_sensorFusion_Post=[];Features_medianFusion_Post=[];Features_Accelerometer_Post=[];Features_Gyroscope_Post_DIAS=[];

Features_Frame_wise_Pre_DIAS=[];Features_sensorFusion_Pre_DIAS=[];Features_medianFusion_Pre_DIAS=[];Features_Accelerometer_Pre_DIAS=[];Features_Gyroscope_Pre_DIAS=[];
Features_Frame_wise_Post_DIAS=[];Features_sensorFusion_Post_DIAS=[];Features_medianFusion_Post_DIAS=[];Features_Accelerometer_Post_DIAS=[];Features_Gyroscope_Post=[];
Features_similarity_framwise_pre=[]; Features_similarity_framwise_Post=[];
Features_angle_framwise_post=[];
Features_angle_framwise_pre=[];
Features_envhflf_framwise_post=[];
Features_envhflf_framwise_pre=[];

Features_Frame_wise_Pre_adxl=[]; Features_Frame_wise_Post_adxl=[];
 Features_Frame_wise_PostAVG_adxl=[]; Features_Frame_wise_PostFUSIONmedian_adxl=[];
 Features_envhflf_framwise_post_AVGhflf=[]; Features_envhflf_framwise_post_FUSIONMedianhflf=[];
  Features_Frame_wise_preAVG_adxl=[]; Features_Frame_wise_preFUSIONmedian_adxl=[];
 Features_envhflf_framwise_pre_AVGhflf=[]; Features_envhflf_framwise_pre_FUSIONMedianhflf=[];

ID=0;
labels=[];
subID=[];
for i=1:2
    switch i
        case 1
            
            % %             STEMI
            hfiles = 'D:\BIOSIGNAL_PROCESSING\STENOSIS STUDY\SAFE STEMI\New CLEANED_ byOLLI\STEMI/'; %% ALL STEMIs with HQ ECG
            %             NSTEMI data
%              hfiles = 'D:\Seafile\Seafile\Saeed_And_Mojtaba_Underwork\MI_ASRO\NEW CLEANED _ 15012019\STEMI_HQekg_noAF/';
            %             STEMI  data only
            %                                     hfiles = 'D:\BIOSIGNAL_PROCESSING\Physiological Recordings\STEMI_Cleaned\MI_group_pairmatched - HQekg\MI_prePCI/';
            %                                     hfiles = 'D:\BIOSIGNAL_PROCESSING\Physiological Recordings\STEMI_Cleaned\MI_group_pairmatched - HQekg\POST_POST_study\post_post/';
            
            Is_MI=1;
            filepattern=fullfile(hfiles,'*.mat')
            matfiles=dir(filepattern)
            %             Features_pre=zeros(length(matfiles),66);
            
        case 2
            %             all controls (postpost)
%             hfiles='D:\BIOSIGNAL_PROCESSING\MI Detect project\STEMI_Mechanical Complexity\DATA_CLEANED_04092018\NEW CLEANED _ 15012019\CONTROL_HQekg_noAF/'
%                          hfiles = 'D:\BIOSIGNAL_PROCESSING\MI Detect project\STEMI_Mechanical Complexity\DATA_CLEANED_04092018\NEW CLEANED _ 15012019\CONTROLS_HIGHQUALITY_SELECTED/';
            %             STEMI Control
                        hfiles = 'D:\BIOSIGNAL_PROCESSING\STENOSIS STUDY\SAFE STEMI\New CLEANED_ byOLLI\PostSTEMI/';
            %             STEMI post
            %                         hfiles = 'D:\BIOSIGNAL_PROCESSING\MI Detect project\STEMI_Mechanical Complexity\DATA_CLEANED_04092018\POST STEMI - without AF/';
            Is_MI=0;
            filepattern=fullfile(hfiles,'*.mat')
            matfiles=dir(filepattern)
            %             Features_post=zeros(length(matfiles),66);
        case 3
            Is_MI=-1;
            
        otherwise
            disp('Unknown path.')
            
    end
    
    
    for k=1:length(matfiles)
        %             try
%      for k=1:1
        fs=200; % sampling freqeuncy
        gr=0; %graphic
        filt=1;  % 1= normal BPF 2= Envelope 3= SSA
        %           data=import_smartphone_data(hfiles, k,filt,fs, gr );% load smartphone data
        baseFileName = matfiles(k).name;
        load([fullfile(hfiles, baseFileName)]);
        %         tt=(1:numel(data.accX))/fs;
        
        
        fs_lsm=200;
        fs_adxl=400;
        jitter=floor(min(length(data.accX)/fs_lsm,length(data.accX_adxl)/fs_adxl));
        
        gyroX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,3,30);
        gyroY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,3,30);
        gyroZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,3,30);
        accX=fft_filter(data.accX(1:jitter*fs_lsm),fs_lsm,5,30);
        accY=fft_filter(data.accY(1:jitter*fs_lsm),fs_lsm,5,30);
        accZ=fft_filter(data.accZ(1:jitter*fs_lsm),fs_lsm,5,30);
        ekg_lf=fft_filter(data.ekg(1:jitter*fs_lsm),fs_lsm,1,49);
        
        data_lsm=struct('ekg',ekg_lf,'accX',accX,'accY',accY,'accZ',accZ,'gyroX',gyroX,'gyroY',gyroY,'gyroZ',gyroZ);
        
        
        
        accX_hf=fft_filter(data.accX_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        accY_hf=fft_filter(data.accY_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        accZ_hf=fft_filter(data.accZ_adxl(1:jitter*fs_adxl),fs_adxl,10,80);
        accX_lf=fft_filter(data.accX_adxl(1:jitter*fs_adxl),fs_adxl,1,10);
        accY_lf=fft_filter(data.accY_adxl(1:jitter*fs_adxl),fs_adxl,1,10);
        accZ_lf=fft_filter(data.accZ_adxl(1:jitter*fs_adxl),fs_adxl,1,10);
        ekg_hf=fft_filter(data.ekg_hf(1:jitter*fs_adxl),400,1,49);
                
        data_adxl=struct('ekg',ekg_hf,'accX',accX_hf,'accY',accY_hf,'accZ',accZ_hf,'accXLF',accX_lf,'accYLF',accY_lf,'accZLF',accZ_lf);
        
%        adxl_energy(data_adxl,fs_adxl,1,k,Is_MI)
% %        close all
%          pause
%          continue
        %%
        simscore=0;
        if simscore
            if Is_MI==1
                [averaged_signal] = ensemble_morph_sim(data,fs);
                all_ensembles_pre{k}=struct2cell(averaged_signal)
            else
                [averaged_signal] = ensemble_morph_sim(data,fs);
                all_ensembles_post{k}=struct2cell(averaged_signal)
            end
            
            continue
            [similarity_Scores,pvals] = similarity_scores(all_ensembles_pre,all_ensembles_post);
            col=@(x)reshape(x,numel(x),1);
            boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
                cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
                'Labels',{'Acc_x','Acc_y','Acc_z','Gyro_x','Gyro_y','Gyro_z'},'Whisker',1);
            figure
            boxplot2({similarity_Scores.SimScor_accX,similarity_Scores.SimScor_accY,similarity_Scores.SimScor_accZ,similarity_Scores.SimScor_gyroX,similarity_Scores.SimScor_gyroY, similarity_Scores.SimScor_gyroZ});
            % title('Rotational Complexity Ratio ')
            %  title(' Similarity Score Between Healthy Subjects  (Serial Measurement)')
            %   title(' Similarity Score Between STEMI Subjects  (Pre PCI vs Post PCI)')
            xlabel('Sensing Modality and Motion Axis');
            ylabel('Pearsons Correlation Coefficient (r)');
            %                  close all
            pval_ax= mean(pvals.Pval_accX)
            pval_ay=mean(pvals.Pval_accY)
            pval_az=mean(pvals.Pval_accZ)
            pval_gx=mean(pvals.Pval_gyroX)
            pval_gy=mean(pvals.Pval_gyroY)
            pval_gz=mean(pvals.Pval_gyroZ)
            
        end
        
        
        %% Feature Extraction
        win=10;
        step=5;
        [ features,features_avgerage_sensors,features_fusion_median,features_acc,features_gyro,labID_ext] =...
            Feature_Extraction_MI_Complexity(data_lsm,fs_lsm,win,step,Is_MI, k,'systole');
        
           [ features_diastole,features_avgerage_sensors_diastole,features_fusion_median_diastole,...
            features_acc_diastole,features_gyro_diastole,labID_ext_diastole] =...
            Feature_Extraction_MI_Complexity(data_lsm,fs_lsm,win,step,Is_MI, k,'diastole');
        
        
        features_angle = angle_features(data_lsm, fs_lsm, win, step,'fullcycle');
        [features_compMob] = energy_complexity_mobility_segmented(data_lsm, fs_lsm, win, step,'fullcycle');
     
        [ features_adxl,features_avg_adxl,features_fusion, feat_accHF, feat_accLF ] =...
            Feature_Extraction_MI_Complexity_ADXL(data_adxl,fs_adxl,win, step,Is_MI,k,'fullcycle');
%         [ features_adxl,features_avg_adxl,features_fusion, feat_accHF, feat_accLF ] = Feature_Extraction_ADXL( data_adxl,fs_adxl,win, step,Is_MI,k );
        [ features_env_hflf,features_avg_envhflf,features_fusion_median_envhflf,labID_ext_envhflf] = Feature_Extraction_MI_Complexity_HIGHFREQ_ENVELOPE(data_adxl,fs_adxl,win,step,Is_MI, k);
        
        if Is_MI==1
            
            %%%%%%%%%SYSTOLIC FEATURES
            Features_Frame_wise_Pre=[Features_Frame_wise_Pre;features];
            Features_sensorFusion_Pre=[Features_sensorFusion_Pre;features_avgerage_sensors];
            Features_medianFusion_Pre=[Features_medianFusion_Pre;features_fusion_median];
            Features_Accelerometer_Pre=[Features_Accelerometer_Pre;features_acc];
            Features_Gyroscope_Pre=[Features_Gyroscope_Pre;features_gyro];
            
            %%%%%% DIASTOLIC FEATURES
            Features_Frame_wise_Pre_DIAS=[Features_Frame_wise_Pre_DIAS;features_diastole];
            Features_sensorFusion_Pre_DIAS=[Features_sensorFusion_Pre_DIAS;features_avgerage_sensors_diastole];
            Features_medianFusion_Pre_DIAS=[Features_medianFusion_Pre_DIAS;features_fusion_median_diastole];
            Features_Accelerometer_Pre_DIAS=[Features_Accelerometer_Pre_DIAS;features_acc_diastole];
            Features_Gyroscope_Pre_DIAS=[Features_Gyroscope_Pre_DIAS;features_gyro_diastole];
            
            
            Features_energy_complexity_pre=[Features_energy_complexity_pre;features_compMob'];
%             Features_similarity_framwise_pre=[Features_similarity_framwise_pre;similarity_score_framwise];
            Features_angle_framwise_pre=[Features_angle_framwise_pre;features_angle];
            
        
            
                
            %%%%%% ADXL sensor FEATURES
            Features_Frame_wise_Pre_adxl=[Features_Frame_wise_Pre_adxl;features_adxl];
            Features_envhflf_framwise_pre=[Features_envhflf_framwise_pre;features_env_hflf];
             Features_Frame_wise_PreAVG_adxl=[Features_Frame_wise_preAVG_adxl;features_avg_adxl];
            Features_Frame_wise_preFUSIONmedian_adxl=[Features_Frame_wise_preFUSIONmedian_adxl;features_fusion];
             Features_envhflf_framwise_pre_AVGhflf=[Features_envhflf_framwise_pre_AVGhflf;features_avg_envhflf];
            Features_envhflf_framwise_pre_FUSIONMedianhflf=[Features_envhflf_framwise_pre_FUSIONMedianhflf;features_fusion_median_envhflf];
            ID=ID+1;
            subID=[subID,ID];
            labels=[labels,1];
            
        elseif Is_MI==0
            
            
            %%%%%LSM sensor systolic features POST/control
            Features_Frame_wise_Post=[Features_Frame_wise_Post;features];
            Features_sensorFusion_Post=[Features_sensorFusion_Post;features_avgerage_sensors];
            Features_medianFusion_Post=[Features_medianFusion_Post;features_fusion_median];
            Features_Accelerometer_Post=[Features_Accelerometer_Post;features_acc];
            Features_Gyroscope_Post=[Features_Gyroscope_Post;features_gyro];
   

            
            %%%%%% LSM sensor diastolic features POST/control
            Features_Frame_wise_Post_DIAS=[Features_Frame_wise_Post_DIAS;features_diastole];
            Features_sensorFusion_Post_DIAS=[Features_sensorFusion_Post_DIAS;features_avgerage_sensors_diastole];
            Features_medianFusion_Post_DIAS=[Features_medianFusion_Post_DIAS;features_fusion_median_diastole];
            Features_Accelerometer_Post_DIAS=[Features_Accelerometer_Post_DIAS;features_acc_diastole];
            Features_Gyroscope_Post_DIAS=[Features_Gyroscope_Post_DIAS;features_gyro_diastole];
            
            
            %%%%%% LSM Angle and energy features
            Features_energy_complexity_Post=[Features_energy_complexity_Post;features_compMob'];
            Features_angle_framwise_post=[Features_angle_framwise_post;features_angle];
            
            
            %%%%%%% ADXL Sensor features POST
            Features_Frame_wise_Post_adxl=[Features_Frame_wise_Post_adxl;features_adxl];
            Features_Frame_wise_PostAVG_adxl=[Features_Frame_wise_PostAVG_adxl;features_avg_adxl];
            Features_Frame_wise_PostFUSIONmedian_adxl=[Features_Frame_wise_PostFUSIONmedian_adxl;features_fusion];
            Features_envhflf_framwise_post=[Features_envhflf_framwise_post;features_env_hflf];
            Features_envhflf_framwise_post_AVGhflf=[Features_envhflf_framwise_post_AVGhflf;features_avg_envhflf];
            Features_envhflf_framwise_post_FUSIONMedianhflf=[Features_envhflf_framwise_post_FUSIONMedianhflf;features_fusion_median_envhflf];

            
            
            ID=ID+1;
            subID=[subID,ID];
            labels=[labels,0];
        end
        
        
    end
    
    
end

Features_Frame_wise_Post(:,end)=Features_Frame_wise_Post(:,end)+max(Features_Frame_wise_Pre(:,end));

Features_pre=[Features_angle_framwise_pre,Features_energy_complexity_pre,Features_envhflf_framwise_pre(:,1:end-2),Features_Frame_wise_Pre_adxl(:,1:end-2),Features_Frame_wise_Pre];
Features_post=[Features_angle_framwise_post,Features_energy_complexity_Post,Features_envhflf_framwise_post(:,1:end-2),Features_Frame_wise_Post_adxl(:,1:end-2),Features_Frame_wise_Post];
% Features_pre=[Features_Frame_wise_Pre];
% Features_post=[Features_Frame_wise_Post];

Feature_matrix_MI_pre_vs_post_SYS=[Features_pre;Features_post];%%Systolic feature
Feature_matrix_MI_pre_vs_post_DIASTOLE=[Features_Frame_wise_Pre_DIAS(:,1:end-2);Features_Frame_wise_Post_DIAS(:,1:end-2)];%%distolic features
Features_SYS_DIAS=[Feature_matrix_MI_pre_vs_post_DIASTOLE,Feature_matrix_MI_pre_vs_post_SYS]; %% Systolic and diastolic data fusion


Feature_matrix_MI_SensorFUSION_sys=[Features_sensorFusion_Pre;Features_sensorFusion_Post];%%Systolic feature fusion
Feature_matrix_MI_SensorFUSION_dias=[Features_sensorFusion_Pre_DIAS;Features_sensorFusion_Post_DIAS];%%Diastolic feature fusion
Features_SYS_DIAS_SenFused=[Feature_matrix_MI_SensorFUSION_dias,Feature_matrix_MI_SensorFUSION_sys];

Features_dualsen_fused_average=[Features_envhflf_framwise_pre_AVGhflf;Features_envhflf_framwise_post_AVGhflf];

Features_lsm_adxl_avg=[Features_dualsen_fused_average,Features_SYS_DIAS_SenFused];

Features_Median_MI_sys=[Features_medianFusion_Pre;Features_medianFusion_Post];
Features_Median_MI_dias=[Features_medianFusion_Pre_DIAS;Features_medianFusion_Post_DIAS];
Features_Median_MI_sysdias=[Features_Median_MI_sys,Features_Median_MI_dias];
Features_dualsen_fused_median=[Features_envhflf_framwise_pre_FUSIONMedianhflf;Features_envhflf_framwise_post_FUSIONMedianhflf];
Features_Median_MI_dualsen=[Features_dualsen_fused_median,Features_Median_MI_sysdias];

%%%%% NORMALIZATION
Xtrain_SYS_DISA=zscore(Features_SYS_DIAS(:,1:end-2));
Xtrain_sys=zscore(Feature_matrix_MI_pre_vs_post_SYS(:,1:end-2));
Xtrain_SenFused_Sys=zscore(Feature_matrix_MI_SensorFUSION_sys);
Xtrain_SenFused_SysDias=zscore(Features_SYS_DIAS_SenFused);
Xtrain_dualsen_lsm_adxl=zscore(Features_lsm_adxl_avg);


Xtrain_Median=zscore(Features_Median_MI_sys);
Xtrain_Dias=zscore(Feature_matrix_MI_pre_vs_post_DIASTOLE);
Xtrain_Median_SysDias=zscore(Features_Median_MI_sysdias);
Xtrain_Median_dualsen=zscore(Features_Median_MI_dualsen);
%%%% LABELLING AND IDs
all_labels=Feature_matrix_MI_pre_vs_post_SYS(:,end-1);
Subjects_IDs=Feature_matrix_MI_pre_vs_post_SYS(:,end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%
% [ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_new2(Subjects_IDs, all_labels', Xtrain_SenFused_Sys , 0, 0);
[ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(Subjects_IDs, all_labels', Xtrain_sys , 0, 0);

[ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(Subjects_IDs, all_labels', Xtrain_SYS_DISA , 0, 0);
[ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(subID, labels, Xtrain_Median , 0, 0);

[ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(Subjects_IDs, all_labels', Xtrain_dualsen_lsm_adxl , 0, 0);
[ ConfMtx_MULTICLASS_svmk_hfsig, ConfMtx_MULTICLASS_svmk_entity_hfsig, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(subID, labels, Xtrain_Median_dualsen , 0, 0);



%%%%%%%%%%%%%%%%%%%%%%
[results,missclassfied]= majorityVoting_Classifier(Xtrain_sys,all_labels,Subjects_IDs,0);
[results,missclassfied]= majorityVoting_Classifier(Xtrain_SenFused_SysDias,all_labels,Subjects_IDs,0);
[results,missclassfied]= majorityVoting_Classifier(Xtrain_dualsen_lsm_adxl,all_labels,Subjects_IDs,0);


[results,missclassfied]= majorityVoting_Classifier(Xtrain_Median_dualsen,labels',subID,0); %% THIS IS THE BEST ONE


% [results,missclassfied]= majorityVotingNew(Xtrain_Median,labels',subID,0);
% [results,missclassfied]= majorityVotingNew(Xtrain_Median_Dias,labels',subID,0);
% [results,missclassfied]= majorityVotingNew(Xtrain_Median_SysDias,labels',subID,0);


ensemble=0;
statisticaltest=0;
if statisticaltest
    %     grp=labels;
    %         obs=Xtrain;
    %
    %         holdoutCVP = cvpartition(grp,'holdout',1)
    %         dataTrain = obs(holdoutCVP.training,:);
    %         grpTrain = grp(holdoutCVP.training);
    %
    %         grp(grp==0)=2;
    %         dataTrainG1 = dataTrain(grp2idx(grpTrain)==1,:);
    %         dataTrainG2 = dataTrain(grp2idx(grpTrain)==2,:);
    %         [h,p,ci,stat] = ttest2(dataTrainG1,dataTrainG2,'Vartype','unequal')
    [pWil,hWil,bestpvalsWil] = Wilcoxontest(Features_pre(:,1:end-2),Features_post(:,1:end-2))
    figure
    ecdf(pWil);
    xlabel('P value');
    ylabel('CDF value');
    [ap, bp]=(find(pWil<0.05));
    
    
    %  [bestfeat,featureIdxSortbyP] = sort(pWil(bestpvalsWil),2);
    %         [~,featureIdxSortbyP] = sort(p,2); % sort the features
    % %         [~,featureIdxSortbyP] = sort(p(p<0.01),2);
    %                 bestfeats=featureIdxSortbyP(1:10)
end

%%%%%%%%%%%%%%%%%%%%% Box Plots
boxplots=1;
if boxplots
    Gyro_Engergy_Complexity_MI=Features_pre(:,1);
    Gyro_Engergy_Complexity_PostMI=Features_post(:,1);
    Gyro_Engergy_Mobility_MI=Features_pre(:,18);
    Gyro_Engergy_Mobility_postMI=Features_post(:,18);
    col=@(x)reshape(x,numel(x),1);
    boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
        cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
        'Labels',{'STEMI','CONTROL'},'Whisker',1);
    figure
    boxplot2({Gyro_Engergy_Complexity_MI,Gyro_Engergy_Complexity_PostMI});
    % title('Rotational Complexity Ratio ')
    title(' Total Mechanical Complexity Dispersion  ')
    % title('Rotational Energy Complexity Index')
    
    figure
    % boxplot([Gyro_Engergy_Mobility_MI,Gyro_Engergy_Mobility_postMI],'Labels',{'MI','post MI'},'Whisker',1)
    boxplot2({Gyro_Engergy_Mobility_MI,Gyro_Engergy_Mobility_postMI});
    % title(' Mechanical Energy Mobility Deviation  ')
    % title('Rotational Complexity Dispersion ')
    title('Mechanical Energy Mobility Deviation ')
end

%%%%%%%%%%%%%%%%%%%%%%%%

% col=@(x)reshape(x,numel(x),1);
% boxplot2=@(C,varargin)boxplot(cell2mat(cellfun(col,col(C),'uni',0)),...
%     cell2mat(arrayfun(@(I)I*ones(numel(C{I}),1),col(1:numel(C)),'uni',0)),varargin{:},...
%     'Labels',{'STEMI','Post PCI','Post PCI (follow-up)','Control'},'Whisker',1);
%  boxplot2({TDMD_stemi,TDMD_poststemi,TDMD_postpost,TDMD_control});
% title(' Total Mechanical Complexity Dispersion  ')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CLASSIFICATION
% [results,missclassfied]= majorityVotingMIstudy(Xtrain_sys,labels,subID,0);
%%
if ensemble
    T = 500;
    treeStump = templateTree('MaxNumSplits',1);
    adaStump = fitcensemble(Xtrain_sys,labels,'Method','AdaBoostM1','NumLearningCycles',T, ...
        'Learners',treeStump);
    totalStump = fitcensemble(Xtrain_sys,labels,'Method','TotalBoost','NumLearningCycles',T, ...
        'Learners',treeStump);
    lpStump = fitcensemble(Xtrain_sys,labels,'Method','LPBoost','NumLearningCycles',T, ...
        'Learners',treeStump);
    
    
    figure;
    plot(resubLoss(adaStump,'Mode','Cumulative'));
    hold on
    plot(resubLoss(totalStump,'Mode','Cumulative'),'r');
    plot(resubLoss(lpStump,'Mode','Cumulative'),'g');
    hold off
    xlabel('Number of stumps');
    ylabel('Training error');
    legend('AdaBoost','TotalBoost','LPBoost','Location','NE');
    
    
    cvlp = crossval(lpStump,'leaveout', 'on');
    cvtotal = crossval(totalStump,'leaveout','on');
    cvada = crossval(adaStump,'leaveout','on');
    
    figure;
    plot(kfoldLoss(cvada,'Mode','Cumulative'));
    hold on
    plot(kfoldLoss(cvtotal,'Mode','Cumulative'),'r');
    plot(kfoldLoss(cvlp,'Mode','Cumulative'),'g');
    hold off
    xlabel('Ensemble size');
    ylabel('Cross-validated error');
    legend('AdaBoost','TotalBoost','LPBoost','Location','NE');
    
    %%%%%%%%%%%%%%%%%%%%%%%
    X=Xtrain_sys;
    Y=labels;
    rng(1); % For reproducibility
    MdlDeep = fitctree(X,Y,'CrossVal','on','MergeLeaves','off',...
        'MinParentSize',1);
    MdlStump = fitctree(X,Y,'MaxNumSplits',1,'CrossVal','on');
    
    n = size(X,1);
    m = floor(log(n - 1)/log(3));
    learnRate = [0.1 0.25 0.5 1];
    numLR = numel(learnRate);
    maxNumSplits = 3.^(0:m);
    numMNS = numel(maxNumSplits);
    numTrees = 150;
    Mdl = cell(numMNS,numLR);
    
    for k = 1:numLR;
        for j = 1:numMNS;
            t = templateTree('MaxNumSplits',maxNumSplits(j));
            Mdl{j,k} = fitcensemble(X,Y,'NumLearningCycles',numTrees,...
                'Learners',t,'leaveout','on','LearnRate',learnRate(k));
        end;
    end;
    
    kflAll = @(x)kfoldLoss(x,'Mode','cumulative');
    errorCell = cellfun(kflAll,Mdl,'Uniform',false);
    error = reshape(cell2mat(errorCell),[numTrees numel(maxNumSplits) numel(learnRate)]);
    errorDeep = kfoldLoss(MdlDeep);
    errorStump = kfoldLoss(MdlStump);
    mnsPlot = [1 round(numel(maxNumSplits)/2) numel(maxNumSplits)];
    figure;
    
    for k = 1:3;
        subplot(2,2,k);
        plot(squeeze(error(:,mnsPlot(k),:)),'LineWidth',2);
        axis tight;
        hold on;
        h = gca;
        plot(h.XLim,[errorDeep errorDeep],'-.b','LineWidth',2);
        plot(h.XLim,[errorStump errorStump],'-.r','LineWidth',2);
        plot(h.XLim,min(min(error(:,mnsPlot(k),:))).*[1 1],'--k');
        h.YLim = [0 0.8];
        xlabel 'Number of trees';
        ylabel 'Cross-validated misclass. rate';
        title(sprintf('MaxNumSplits = %0.3g', maxNumSplits(mnsPlot(k))));
        hold off;
    end
    
    hL = legend([cellstr(num2str(learnRate','Learning Rate = %0.2f'));...
        'Deep Tree';'Stump';'Min. misclass. rate']);
    hL.Position(1) = 0.6;
    
end

% [p,h,stats]=ranksum(all_mccr_preMI,all_mccr_postMI,'method','approximate')
% [p,h,stats]=ranksum(all_mcct_preMI,all_mcct_postMI)
% [p,h,stats]=ranksum(all_mccrSD_preMI,all_mccrSD_postMI)
% [p,h,stats]=ranksum(all_mcctSD_preMI,all_mcctSD_postMI)

%%%%%
% X_TEST_DL=Features_SYS_DIAS_SenFused(1180:2224,:);
% X_TRAIN_DL=Features_SYS_DIAS_SenFused;
% X_TRAIN_DL(1180:2224,:)=[];
% X_TRAIN_DL_norm=zscore(X_TRAIN_DL);
% X_TEST_DL_norm=(X_TEST_DL-mean(X_TRAIN_DL))./std(X_TRAIN_DL);
% Y_TEST_DL=Features_SYS_DIAS(1180:2224,end-1);
% Y_TRAIN_DL=Features_SYS_DIAS(:,end-1);
% Y_TRAIN_DL(1180:2224,:)=[];
