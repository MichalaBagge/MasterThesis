%% Features and other ML methods
% Restart session
close all; clc;
clearvars -except subjects

% Define directories
dir_data = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Data\Eye movements con-glo Julie\MatlabData';
dir_annotations = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Data\Eye movements con-glo Julie\EM manual scorings 2019';
dir_feature_funcs = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Code\MATLAB\Feature_functions';
dir_scripts = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Code\MATLAB\Workspace_201224';

% Load subject IDs
cd(dir_data)
F = dir(dir_data); F = F(3:end);
F = char(F(:).name);

% Load subject feature table
cd(dir_scripts)
dir_subject = fullfile(dir_scripts,"AllSubjects_REM_SEM_Feature_Table.csv");
EM_Feature_Table = readtable(dir_subject);
attributeNames = EM_Feature_Table.Properties.VariableNames;
attributeNames = attributeNames(:,3:end-1);
attributeNames = strrep(attributeNames, '_', '\_');

% Classification
% Extract y as response variable
y = EM_Feature_Table.y;
X = table2array(EM_Feature_Table(:,3:end-1));
X_IDs = char(EM_Feature_Table.Subejct_ID);
[~,M] = size(X);
N = length(F);

% Initialize variables
rng(0);
K = 5;
CV = cvpartition(N, 'Kfold', K);

%accuracy = zeros(K, 1);
FeatureImportance_ave = zeros(K, M);

% Initiate figures
%figure;

% For each crossvalidation fold
for k = 1:K
    fprintf('Crossvalidation fold %d/%d\n', k, K);
    
    % Extract the training and test set
    train_set = F(CV.training(k), :); 
    test_set = F(CV.test(k), :);

    idx_train = ismember(string(X_IDs), string(train_set));
    idx_test = ismember(string(X_IDs), string(test_set));
    
    X_train = X(idx_train, :); 
    y_train = y(idx_train);
    X_test = X(idx_test, :); 
    y_test = y(idx_test);

    % Train a Random Forest model on the training set
    model = fitcensemble(X_train, y_train, ...
        'Method', 'Bag', ...
        'NumLearningCycles', 100);

    % Test the model on the test set
    [y_est, scores] =  predict(model, X_test);
    

    % Compute feature importance for the current model
    FeatureImportanceCV(k, :) = predictorImportance(model);

    % Assign binary classes
    y_test_nonem = zeros(size(y_test));
    y_test_rem = zeros(size(y_test));
    y_test_sem = zeros(size(y_test));
    y_test_nonem(y_test==0) = 1;
    y_test_rem(y_test==1) = 1;
    y_test_sem(y_test==2) = 1;
   
    % Customize threshold to be more realistic
    custom_thresholds = 0:0.1:1;

%     [FPR_rem, TPR_rem, T_rem] = perfcurve(y_test_rem, scores(:,2), 1, 'xcrit', 'reca', 'ycrit', 'prec'); % 1 is the positive class
%     [FPR_nonem, TPR_nonem, T_nonem] = perfcurve(y_test_nonem, scores(:,1), 1, 'xcrit', 'reca', 'ycrit', 'prec'); % 0 is the positive class
%     [FPR_sem, TPR_sem, T_sem] = perfcurve(y_test_sem, scores(:,3), 1, 'xcrit', 'reca', 'ycrit', 'prec'); % 2 is the positive class

    [FPR_rem, TPR_rem,T_rem, ~,~] = perfcurve(y_test_rem, scores(:,2), 1, 'XVals', custom_thresholds); % 1 is the positive class
    [FPR_nonem, TPR_nonem,T_nonem, ~,~] = perfcurve(y_test_nonem, scores(:,1), 1, 'XVals', custom_thresholds); % 0 is the positive class
    [FPR_sem, TPR_sem, T_sem, ~,~] = perfcurve(y_test_sem, scores(:,3), 1, 'XVals', custom_thresholds); % 2 is the positive class

    var_save = {'scores','y_est','y_test','FPR_rem','TPR_rem','FPR_sem','TPR_sem','FPR_nonem','TPR_nonem','T_nonem','T_rem','T_sem'}; 
    varsave = {'scores','y_est','y_test','FPR_rem','TPR_rem','FPR_sem','TPR_sem','FPR_nonem','TPR_nonem','T_nonem','T_rem','T_sem'};
    for ii = 1:length(var_save)
        eval(sprintf('results_kfolds.%s.k%d = %s ',char(var_save(ii)),k,char(varsave(ii))));
    end

    % Clear temporary variables
    clearvars -except dir_scripts dir_annotations dir_data dir_subject ...
           F k ii results_kfolds y X X_IDs M N K CV ...
           FeatureImportanceCV attributeNames scores


end
% clear command window
clc;
% Average feature importance across all folds
FeatureImportance_ave = mean(FeatureImportanceCV, 1);

%%
% Plot Feature Importance
figure;
bar(FeatureImportance_ave);
title('Average Feature Importance Across Folds');
%xlabel('Feature Index');
ylabel('Importance');
xticks(1:length(attributeNames));
xticklabels(attributeNames);
xtickangle(45);
grid on;


%%

[fig_roc_nonem, opt_th_nonem, best_non_fpr, best_non_tpr] = generat_RF_EM_Part2_plots(results_kfolds.FPR_nonem, results_kfolds.TPR_nonem,results_kfolds.T_nonem,K);
%[fig_roc_rem, best_rem_fpr, best_rem_tpr]  = generat_RF_EM_Part2_PR_plots(results_kfolds.FPR_rem, results_kfolds.TPR_rem,results_kfolds.T_rem,results_kfolds.y_test,1,K);
%[fig_roc_sem, best_sem_fpr, best_sem_tpr]= generat_RF_EM_Part2_PR_plots(results_kfolds.FPR_sem, results_kfolds.TPR_sem,results_kfolds.T_sem,results_kfolds.y_test,2,K);

%%

fig_confmatrix = figure;

for k = 1:K
    subplot(1,5,k)
    eval(sprintf('confusionchart(results_kfolds.y_test.k%d, results_kfolds.y_est.k%d);',k,k));
    eval(sprintf('confusion_matrix = confusionmat(results_kfolds.y_test.k%d, results_kfolds.y_est.k%d);',k,k));
    % Class 0
    TP0 = confusion_matrix(1, 1);
    FP0 = confusion_matrix(2, 1) + confusion_matrix(3, 1);
    FN0 = confusion_matrix(1, 2) + confusion_matrix(1, 3);
    TN0 = confusion_matrix(2, 2) + confusion_matrix(2, 3) + confusion_matrix(3, 2) + confusion_matrix(3, 3);
    % Class 1
    TP1 = confusion_matrix(2, 2);
    FP1 = confusion_matrix(1, 2) + confusion_matrix(3, 2);
    FN1 = confusion_matrix(2, 1) + confusion_matrix(2, 3);
    TN1 = confusion_matrix(1, 1) + confusion_matrix(1, 3) + confusion_matrix(3, 1) + confusion_matrix(3, 3);
    % Class 2
    TP2 = confusion_matrix(3, 3);
    FP2 = confusion_matrix(1, 3) + confusion_matrix(2, 3);
    FN2 = confusion_matrix(3, 1) + confusion_matrix(3, 2);
    TN2 = confusion_matrix(1, 1) + confusion_matrix(1, 2) + confusion_matrix(2, 1) + confusion_matrix(2, 2);

    eval(sprintf('EM_metrics.k%d = EM_performance(TP2, TN2, FP2, FN2);',k));  
    clear TP TN FP FN
end

%

EM_metrics_all = [[EM_metrics.k1(:,1); EM_metrics.k2(:,1); EM_metrics.k3(:,1); EM_metrics.k4(:,1); EM_metrics.k5(:,1)], ...
     [EM_metrics.k1(:,2); EM_metrics.k2(:,2); EM_metrics.k3(:,2); EM_metrics.k4(:,2); EM_metrics.k5(:,2)], ...   
     [EM_metrics.k1(:,3); EM_metrics.k2(:,3); EM_metrics.k3(:,3); EM_metrics.k4(:,3); EM_metrics.k5(:,3)], ...   
     [EM_metrics.k1(:,4); EM_metrics.k2(:,4); EM_metrics.k3(:,4); EM_metrics.k4(:,4); EM_metrics.k5(:,4)], ...   
     [EM_metrics.k1(:,5); EM_metrics.k2(:,5); EM_metrics.k3(:,5); EM_metrics.k4(:,5); EM_metrics.k5(:,5)], ...   
     [EM_metrics.k1(:,6); EM_metrics.k2(:,6); EM_metrics.k3(:,6); EM_metrics.k4(:,6); EM_metrics.k5(:,6)]];   

EM_metrics_mean_all = [mean(EM_metrics_all); std(EM_metrics_all)]


%%
%%% Select featur
n_feat = 4; % 1 2 4 7

% Define example data
groups = {'Non', 'REM', 'SEM'};

X_oneMostIF = X(:,n_feat);
data_size = {X_oneMostIF(y == 0),...
             X_oneMostIF(y == 1),...
             X_oneMostIF(y ==2)};

g = [repmat(groups(1),size(data_size{1,1}));...
    repmat(groups(2),size(data_size{1,2}));...
    repmat(groups(3),size(data_size{1,3}))];

data = [X_oneMostIF(y == 0);...
             X_oneMostIF(y == 1);...
             X_oneMostIF(y ==2)];


figure;
boxplot(data,g)
ylabel(sprintf('%s [Watt]',char(attributeNames(n_feat))));

p=[ranksum(X_oneMostIF(y == 0),X_oneMostIF(y ==1)),...
ranksum(X_oneMostIF(y == 0),X_oneMostIF(y == 2)),...
ranksum(X_oneMostIF(y == 1),X_oneMostIF(y ==2))]


