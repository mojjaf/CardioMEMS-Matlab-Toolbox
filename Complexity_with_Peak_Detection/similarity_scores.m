function [similarity_Scores,pvals,delays] = similarity_scores(all_ensembles_pre,all_ensembles_post,ccflag)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

[r,c]=size(all_ensembles_post)

 
for i=1:c
    if ccflag
    [accX_syn_pre,accX_syn_post,D1(i,:)] = alignsignals(all_ensembles_pre{1,i}{1,1}',all_ensembles_post{1,i}{1,1}');
    [accY_syn_pre,accY_syn_post,D2(i,:)] = alignsignals(all_ensembles_pre{1,i}{2,1}',all_ensembles_post{1,i}{2,1}');
    [accZ_syn_pre,accZ_syn_post,D3(i,:)] = alignsignals(all_ensembles_pre{1,i}{3,1}',all_ensembles_post{1,i}{3,1}');
    [GyrX_syn_pre,GyrX_syn_post,D4(i,:)] = alignsignals(all_ensembles_pre{1,i}{4,1}',all_ensembles_post{1,i}{4,1}');
    [GyrY_syn_pre,GyrY_syn_post,D5(i,:)] = alignsignals(all_ensembles_pre{1,i}{5,1}',all_ensembles_post{1,i}{5,1}');
    [GyrZ_syn_pre,GyrZ_syn_post,D6(i,:)] = alignsignals(all_ensembles_pre{1,i}{6,1}',all_ensembles_post{1,i}{6,1}');
%     
     
       jitter=1000;
    [Sim_Score_AccX(i,:),Pval_accX(i,:)]=corr(accX_syn_pre(abs(D1(i,:))+1:jitter),accX_syn_post(abs(D1(i,:))+1:jitter),'Type','Pearson');
    [Sim_Score_AccY(i,:),Pval_accY(i,:)]=corr(accY_syn_pre(abs(D2(i,:))+1:jitter),accY_syn_post(abs(D2(i,:))+1:jitter),'Type','Pearson');
    [Sim_Score_AccZ(i,:),Pval_accZ(i,:)]=corr(accZ_syn_pre(abs(D3(i,:))+1:jitter),accZ_syn_post(abs(D3(i,:))+1:jitter),'Type','Pearson');
    [Sim_Score_gyroX(i,:),Pval_gyroX(i,:)]=corr(GyrX_syn_pre(abs(D4(i,:))+1:jitter),GyrX_syn_post(abs(D4(i,:))+1:jitter),'Type','Pearson');
    [Sim_Score_gyroY(i,:),Pval_gyroY(i,:)]=corr(GyrY_syn_pre(abs(D5(i,:))+1:jitter),GyrY_syn_post(abs(D5(i,:))+1:jitter),'Type','Pearson');
    [Sim_Score_gyroZ(i,:),Pval_gyroZ(i,:)]=corr(GyrZ_syn_pre(abs(D6(i,:))+1:jitter),GyrZ_syn_post(abs(D6(i,:))+1:jitter),'Type','Pearson');
               delays=struct('Delay_accX',D1,'Delay_accY',D2,'Delay_accZ',D3,'Delay_gyroX',D4,'Delay_gyroY',D5,'Delay_gyroZ',D6);

    else
        jitter_Stop=1000;
       jitter_Start=1;
    [Sim_Score_AccX(i,:),Pval_accX(i,:)]=corr(all_ensembles_pre{1,i}{1,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{1,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    [Sim_Score_AccY(i,:),Pval_accY(i,:)]=corr(all_ensembles_pre{1,i}{2,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{2,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    [Sim_Score_AccZ(i,:),Pval_accZ(i,:)]=corr(all_ensembles_pre{1,i}{3,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{3,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    [Sim_Score_gyroX(i,:),Pval_gyroX(i,:)]=corr(all_ensembles_pre{1,i}{4,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{4,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    [Sim_Score_gyroY(i,:),Pval_gyroY(i,:)]=corr(all_ensembles_pre{1,i}{5,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{5,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    [Sim_Score_gyroZ(i,:),Pval_gyroZ(i,:)]=corr(all_ensembles_pre{1,i}{6,1}(jitter_Start:jitter_Stop)',all_ensembles_post{1,i}{6,1}(jitter_Start:jitter_Stop)','Type','Pearson');
    delays=0;
    end
    similarity_Scores=struct('SimScor_accX',Sim_Score_AccX,'SimScor_accY',Sim_Score_AccY,'SimScor_accZ',Sim_Score_AccZ,'SimScor_gyroX',Sim_Score_gyroX,'SimScor_gyroY',Sim_Score_gyroY,'SimScor_gyroZ',Sim_Score_gyroZ);
      pvals=struct('Pval_accX',Pval_accX,'Pval_accY',Pval_accY,'Pval_accZ',Pval_accZ,'Pval_gyroX',Pval_gyroX,'Pval_gyroY',Pval_gyroY,'Pval_gyroZ',Pval_gyroZ);


end

    
end

