
function [ ConfMtx_MULTICLASS_svmk, ConfMtx_MULTICLASS_svmk_entity, accuracy, accuracyV, sensitivityV, specificityV] = classify_matrix_generator(LABELIT, Y_TULOS, X , Xpca, drawPCAplot)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%  'Starting classification'
 
%*************************************************************************
% Preparations for classification
%*************************************************************************

% --Not needed if classes start from 0 and are in increasing order--
% Make Y_TULOS_UUSI start from 0 in increasing order of Y_TULOS
% Y_TULOS_UUSI=0;
% Y_TULOS_UUSI(1:length(Y_TULOS))=-1;
% running_sum=0;
% for inc_loop=1:6, % starting from the smallest condition
%     if any(Y_TULOS==inc_loop),
%     Y_TULOS_UUSI(Y_TULOS==inc_loop)=running_sum; % running_sum: 0,1,2,3 ...
%     running_sum=running_sum+1;
%     end;
% end;

Y_TULOS_UUSI=0;
Y_TULOS_UUSI=Y_TULOS;

% Make LABELIT as increasing order with respect to Y_TULOS_UUSI 

LABELIT_UUSI=0;
LABELIT_UUSI(1:length(LABELIT))=-1;
inc_loop2_prev=0;

for inc_loop2=unique(Y_TULOS_UUSI), 
    % At first iteration (Y_TULOS_UUSI=0) do not change anything
    LABELIT_UUSI(Y_TULOS_UUSI==inc_loop2)=LABELIT(Y_TULOS_UUSI==inc_loop2)+sum(inc_loop2_prev); 
    % At succeeding iterations, add always the maximum of previous iterations
    % to LABELIT_UUSI to make LABELIT_UUSI unique and increasing
    inc_loop2_prev=cat(2,inc_loop2_prev, [max(LABELIT(Y_TULOS_UUSI==inc_loop2))]); % use all of the previous max:es
end

X_UUSI=0;
X_UUSI=X;% all attributes
LABELIT_UUSI=LABELIT_UUSI';%  all subjects
Y_TULOS_UUSI=Y_TULOS_UUSI';%  all classes assigned to subjects

dt = datestr(now,'mmmm dd yyyy HH-MM-SS-FFF');
%dlmwrite(strcat(dt,'-x.txt'),X);

% Do various experiments here
% X(isnan(X)) = 0; 
% X_UUSI=fastica(X');
% X_UUSI=X_UUSI';

%*************************************************************************
%-------------------------------------------------------------------------
% CLASSIFY ACCORDING TO [X] AND [Y_TULOS_UUSI] [LABELIT_UUSI]
%-------------------------------------------------------------------------
% Leave-one-out person specific cross-validation 
%-------------------------------------------------------------------------
%*************************************************************************

Y_TULOS_UUSI'
LABELIT_UUSI'
%pause;
score_kw=0;
w_tulos_all=0;

for wloop=1:max(LABELIT_UUSI), %  all subjects
    
    'Cross validation iteration (wloop)'
    wloop

    WX_TEST=0;
    WX_TEST=X_UUSI(LABELIT_UUSI==wloop,:); % Choose test set
    WY_TULOS_TEST=0;
    WY_TULOS_TEST=Y_TULOS_UUSI(LABELIT_UUSI==wloop,1); % The correct classes for test set 

    WX=0;    
    WX=X_UUSI(LABELIT_UUSI~=wloop,:); % Choose training set (i.e. all the others that we are not testing)
    WY_TULOS=0;
    WY_TULOS=Y_TULOS_UUSI(LABELIT_UUSI~=wloop,1); % The correct classes for training set 

    
        
    % Multiclass SVM classifier
    
    svm_w=0;

    if max(Y_TULOS_UUSI)<2, % two classes
       % Without cost function

%       svm_w=fitcsvm(WX,WY_TULOS,'Standardize',true,'KernelFunction','gaussian','KernelScale','auto');    

      % With random forest:       
      svm_w=fitensemble(WX,WY_TULOS,'Bag',350,'tree','type','classification');
         
      % With cost function (e.g. try to get as few false classified ami->healthy as possible
%       svm_w=fitcsvm(WX,WY_TULOS,'Standardize',true,'KernelFunction','RBF','KernelScale','auto','Cost',[0 6; 1 0]);    
    else % more than two classes

%----------- ORIG: -------------------------------------------------------
% __DEFAULT__ WHICH IS LINEAR KERNEL:
        t = templateSVM('Standardize',1);        
        svm_w = fitcecoc(WX,WY_TULOS,'Learners',t);
%--------------------------------------------------------------------------
% OR:
%        t = templateSVM('Standardize',1,'KernelFunction','gaussian','KernelScale','auto');        
%        svm_w = fitcecoc(WX,WY_TULOS,'Learners',t);

        % templateSVM('Standardize',1,'KernelFunction','gaussian'); 
        %(also with KERNEL) -> <<'gaussian'>>, 'linear', and 'polynomial')

        % -==RANDOM FOREST==- multiclass:        
        % svm_w=fitensemble(WX,WY_TULOS,'Bag',200,'tree','type','classification'); 
        
        %--------------------------------------------------------------------------
        % Testing polynomial kernel, the default results are with linear kernel for multiclass case and for rbf 
        %(gaussian) kernel in this code in two-class case

%        t = templateSVM('Standardize',1,'KernelFunction','polynomial','KernelScale','auto');        
%        svm_w = fitcecoc(WX,WY_TULOS,'Learners',t);
%--------------------------------------------------------------------------


                    
    end;
    
    score_svm_kw=0;
    w_tulos=0;
    [w_tulos score_svm_kw]=svm_w.predict(WX_TEST);

    size(score_svm_kw)
    score_svm_kw
        
    wlabels=0;
    wlabels=find(LABELIT_UUSI==wloop); % test set %  subjectID equals to  wloop
            
        for wloop2=1:length(wlabels),                    
          score_kw(wlabels(wloop2),1:2)=score_svm_kw(wloop2,1:2);                
          w_tulos_all(wlabels(wloop2))=w_tulos(wloop2);
        end;

end; % wloop 

ConfMtx_MULTICLASS_svmk=0;
ConfMtx_MULTICLASS_svmk = confusionmat(Y_TULOS_UUSI, w_tulos_all);

%*************************************************************************
%-------------------------------------------------------------------------
%                      Perform majority voting 
%-------------------------------------------------------------------------
%*************************************************************************

w_tulos_all_multiclass_entity=0; % "majority voted" w_tulos_all_multiclass
Y_TULOS_UUSI_multiclass_entity=0; % correct classes for "majority voted"

valid_LABEL_vec_indexes=0;

sum11=0;

for aajj=1:max(LABELIT_UUSI),  %  all subjects   
     if (sum(LABELIT_UUSI==aajj))>0,
     sum11=sum11+1;
     valid_LABEL_vec_indexes(sum11)=aajj; % There must be at least one index
     end; 
end;

index_xsum=0;

for joku_wloop=valid_LABEL_vec_indexes, % For all with at least one label index
    
    index_xsum=index_xsum+1;
    
     www_labels=0;
     www_labels=find(LABELIT_UUSI==joku_wloop); % Check this label index       

     Y_TULOS_UUSI_multiclass_entity(index_xsum)=Y_TULOS_UUSI(www_labels(1)); % Assign correct "majority voted" class                  
     
     summa_tot=0;
     summa_tot=length(www_labels);
     summa_odd=0;

     summa_odd=summa_tot; % because of old version of code
     % "summa_odd" is the sum of all labels with current index joku_wloop
     
     %***********************
     %***********************
     %***********************
     
     summa_nollat=0;
     summa_ykkoset=0;
     summa_kakkoset=0;
     summa_kolmoset=0;
     summa_neloset=0;
     summa_vitoset=0;
          
     % 0-1
     summa_nollat=sum(w_tulos_all(www_labels(1:summa_odd))==(0));  
     summa_ykkoset=sum(w_tulos_all(www_labels(1:summa_odd))==(1)); 
     
     if max(Y_TULOS_UUSI)>1, % 2
        summa_kakkoset=sum(w_tulos_all(www_labels(1:summa_odd))==(2));  
     if max(Y_TULOS_UUSI)>2, % 3
        summa_kolmoset=sum(w_tulos_all(www_labels(1:summa_odd))==(3));      
     if max(Y_TULOS_UUSI)>3, % 4
         summa_neloset=sum(w_tulos_all(www_labels(1:summa_odd))==(4));      
     if max(Y_TULOS_UUSI)>4, % 5
         summa_vitoset=sum(w_tulos_all(www_labels(1:summa_odd))==(5));      
     end;
     end;
     end;
     end;
          
     indd_vec=0;
          
     if max(Y_TULOS_UUSI)<2,
     indd_vec=[summa_nollat summa_ykkoset];     
     end;

     if max(Y_TULOS_UUSI)==2,
     indd_vec=[summa_nollat summa_ykkoset summa_kakkoset];            
     end;

     if max(Y_TULOS_UUSI)==3,
     indd_vec=[summa_nollat summa_ykkoset summa_kakkoset summa_kolmoset];                     
     end;

     if max(Y_TULOS_UUSI)==4,
     indd_vec=[summa_nollat summa_ykkoset summa_kakkoset summa_kolmoset summa_neloset];                          
     end;

     if max(Y_TULOS_UUSI)==5,
     indd_vec=[summa_nollat summa_ykkoset summa_kakkoset summa_kolmoset summa_neloset summa_vitoset];                                   
     end;
         
     [indd_a indd_b]=max(indd_vec); % indd_b is the correct class using a majority vote!
     
     w_tulos_all_multiclass_entity(index_xsum)=indd_b-1;
     
end;

accuracy = 0.0;
accuracyV = 0.0;
sensitivityV = 0.0;
specificityV = 0.0;
ConfMtx_MULTICLASS_svmk_entity=0;
ConfMtx_MULTICLASS_svmk_entity = confusionmat(Y_TULOS_UUSI_multiclass_entity, w_tulos_all_multiclass_entity);
[rows,cols] = size(ConfMtx_MULTICLASS_svmk_entity);

%output variables are mentioned also
'Accuracy without majority voting'
accuracy = sum(diag(ConfMtx_MULTICLASS_svmk))/sum(sum(ConfMtx_MULTICLASS_svmk))

'Accuracy with majority voting'
accuracyV = sum(diag(ConfMtx_MULTICLASS_svmk_entity))/sum(sum(ConfMtx_MULTICLASS_svmk_entity))
if rows == 2    
    'Specificity with majority voting'
    specificityV = ConfMtx_MULTICLASS_svmk_entity(1,1)/(ConfMtx_MULTICLASS_svmk_entity(1,1)+ConfMtx_MULTICLASS_svmk_entity(1,2))
    'Sensitivity with majority voting'
    sensitivityV = ConfMtx_MULTICLASS_svmk_entity(2,2)/(ConfMtx_MULTICLASS_svmk_entity(2,1)+ConfMtx_MULTICLASS_svmk_entity(2,2))

    'Confusion matrix without majority voting'
    ConfMtx_MULTICLASS_svmk

    'Confusion matrix with majority voting'
    ConfMtx_MULTICLASS_svmk_entity
elseif rows > 2    
    %A = [30 20 10; 50 60 10;20 20 80]
    % A = [25 5 2;3 32 4;1 0 15]
    rowSum = zeros(1,rows);
    precision = zeros(1,rows);
    colSum = zeros(1,cols);
    sensitivity = zeros(1,cols); % recall
    specificity = zeros(1,rows);
    accuracyY = zeros(1,rows);
    rowSum3 = zeros(1,rows);

    diagonals = diag(ConfMtx_MULTICLASS_svmk_entity);% true positives
    for i=1:rows
        rowSum(i) = sum(ConfMtx_MULTICLASS_svmk_entity(i,:));
    end
    for i=1:cols
        colSum(i) = sum(ConfMtx_MULTICLASS_svmk_entity(:,i));
    end
    for i=1:rows
        sensitivity(i) = diagonals(i)/rowSum(i);
    end
    for i=1:cols
        precision(i) = diagonals(i)/colSum(i);
    end
    rowSum2 = rowSum;
    for i=1:rows
        rowSum2(i) = rowSum2(i)-diagonals(i);%false negative
    end    
    colSum2 = colSum;    
    for i=1:cols
        colSum2(i) = colSum2(i)-diagonals(i);%false positive
    end    
    for i=1:rows
        rowSum3(i) = sum(diagonals)-diagonals(i);%true negative
    end    
    for i=1:rows
         accuracyY(i) = (diagonals(i)+rowSum3(i))/(diagonals(i)+rowSum2(i)+rowSum3(i)+colSum2(i));
    end
    for i=1:rows
         specificity(i) = rowSum3(i)/(rowSum3(i)+colSum2(i));
    end    
    accuracyY
    sensitivity
    specificity
    %sum(diagonals) = true positives + true negatives
    
%     % Another way to calculate specificiity
%     true_neg = zeros(1,rows);
%     avoidRow=0;avoidCol=0;
%     for k=1:rows
%         avoidRow=k;avoidCol=k;
%         for i=1:rows
%             if i ~= avoidRow
%                 for j=1:cols 
%                     if j ~= avoidCol
%                         true_neg(k) = true_neg(k) + A(i,j);
%                     end
%                 end
%             end
%         end
%     end
%     for i=1:rows
%          specificity(i) = true_neg(i)/(true_neg(i)+colSum2(i));
%     end
    
end; 


