%% Features and other ML methods
% Restart session
close all; clc;
clearvars -except subjects

% Define directories
dir_scripts = '';

% Load subject feature table
cd(dir_scripts)
dir_subject = fullfile(dir_scripts,"FeatureTable_RBDPD_Part3___.csv");
dir_labels = fullfile(dir_scripts,"keys_anonymized.csv");
y_temp = readtable(dir_labels);
idx_y = ismember(y_temp.SubjectID,EM_Feature_Table.SubjectID);
%y = sum([y_temp.PDwoRBD(idx_y), y_temp.PDwRBD(idx_y)],2);
C = sum([y_temp.Control(idx_y),y_temp.PLM(idx_y)],2);
y = y_temp.RBD(idx_y); % y
PDw = y_temp.PDwRBD(idx_y);
PDwo = y_temp.PDwoRBD(idx_y);
% %
% EM_Feature_Table(RBD ==1,:) = [];
% y(RBD ==1) = [];
%%% CvsRBD
eli = sum([y_temp.PDwoRBD(idx_y), y_temp.PDwRBD(idx_y)],2);
 EM_Feature_Table(eli ==1,:) = [];
 y(eli ==1,:) = [];
%%% 
% EM_Feature_Table(PDwo ==1,:) = [];
% y = y_temp.PDwRBD(idx_y);
% y(PDwo ==1,:) = [];
% C(PDwo ==1,:) = [];
% PDw(PDwo ==1,:) = [];
% EM_Feature_Table(RBD ==1,:) = [];
% y(RBD ==1,:) = [];
% C(RBD ==1,:) = [];
% PDwo(RBD ==1,:) = [];
% %%%
% EM_Feature_Table(PDw ==1,:) = [];
% y = y_temp.PDwoRBD(idx_y);
% y(PDw ==1,:) = [];
% C(PDw ==1,:) = [];
% PDwo(PDw ==1,:) = [];
% EM_Feature_Table(RBD ==1,:) = [];
% y(RBD ==1,:) = [];
% C(RBD ==1,:) = [];
% PDwo(RBD ==1,:) = [];
%%%
attributeNames = EM_Feature_Table.Properties.VariableNames;
attributeNames = attributeNames(:,2:end);
attributeNames = strrep(attributeNames, '_', '\_');
% Clear temporary variables
%clear idx_y y_temp

% Classification
% Extract y as response variable
X = table2array(EM_Feature_Table(:,2:end));
X_IDs = char(EM_Feature_Table.SubjectID);
[N,M] = size(X);


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
    train_set = X_IDs(CV.training(k), :); 
    test_set = X_IDs(CV.test(k), :);

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
    [FPR, TPR, T,~,~] = perfcurve(y_test, scores(:,2), true); % 1 is the positive class

    var_save = {'scores','y_est','y_test','FPR','TPR','T'}; 
    varsave = {'scores','y_est','y_test','FPR','TPR','T'}; 
    for ii = 1:length(var_save)
        eval(sprintf('results_kfolds.%s.k%d = %s ',char(var_save(ii)),k,char(varsave(ii))));
    end

    % Clear temporary variables
    clearvars -except dir_scripts dir_annotations dir_data dir_subject ...
           F k ii results_kfolds y X X_IDs M N K CV ...
           FeatureImportanceCV attributeNames scores ...
           C RBD PDwo PDw EM_Feature_Table

end

%% Average feature importance across all folds
FeatureImportance_ave = mean(FeatureImportanceCV, 1);

[FeatureImportance_sorted, idx_FI]= sort(FeatureImportance_ave, 'descend');

% Plot Feature Importance
figure;
bar(FeatureImportance_sorted(1:10));
title('Average Feature Importance Across Folds');
ax = gca; 
ax.FontSize = 14;
xlabel('Feature Index');
ylabel('Importance');
xticks(1:length(idx_FI));
xticklabels(attributeNames(idx_FI(1:10)));
xtickangle(45);
grid on;


%%

[fig_roc, opt_th_nonem, best_fpr, best_tpr] = generat_RF_EM_Part1_plots(results_kfolds.FPR, results_kfolds.TPR,results_kfolds.T,K);

%%

fig_confmatrix = figure;

for k = 1:K
    subplot(1,5,k)
    eval(sprintf('confusionchart(results_kfolds.y_test.k%d, results_kfolds.y_est.k%d);',k,k));
    eval(sprintf('confusion_matrix = confusionmat(results_kfolds.y_test.k%d, results_kfolds.y_est.k%d);',k,k));
    TP = confusion_matrix(2, 2);
    FP = confusion_matrix(1, 2);
    FN = confusion_matrix(2, 1);
    TN = confusion_matrix(1, 1);
    eval(sprintf('EM_metrics.k%d = EM_performance(TP, TN, FP, FN);',k));  
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

correlationmatrix = corr(X);

% Create the colormap
num_colors = 256; % Number of colors in the colormap
cmap = [linspace(0, 1, num_colors)', linspace(0, 1, num_colors)', ones(num_colors, 1)]; % Blue to White
cmap = [cmap; [ones(num_colors, 1), linspace(1, 0, num_colors)', linspace(1, 0, num_colors)']]; % White to Red

figure;
% Create Heatmap
h = heatmap(attributeNames,attributeNames,correlationmatrix);
colormap(cmap)
caxis([-1 1]);


%%
%%% Select featur
n_feat = 8; % 1,2, 

% Define example data
groups = {'C', 'RBD', 'PDwo', 'PDw'};

X_oneMostIF = X(:,idx_FI(n_feat));
data_size = {X_oneMostIF(C == 1),...
             X_oneMostIF(RBD == 1),... 
             X_oneMostIF(PDwo ==1),...
             X_oneMostIF(PDw ==1)};

g = [repmat(groups(1),size(data_size{1,1}));...
    repmat(groups(2),size(data_size{1,2}));...
    repmat(groups(3),size(data_size{1,3}));...
    repmat(groups(4),size(data_size{1,4}))];

data = [X_oneMostIF(C == 1);... 
        X_oneMostIF(RBD == 1);...
        X_oneMostIF(PDwo ==1);...
        X_oneMostIF(PDw ==1)];



figure;
boxplot(data,g)
ylabel(attributeNames(idx_FI(n_feat)));
ylabel(sprintf('%s [Hz]',char(attributeNames(idx_FI(n_feat)))));
ax = gca; 
ax.FontSize = 12;


p=[ranksum(X_oneMostIF(C == 1),X_oneMostIF(RBD ==1)),...
ranksum(X_oneMostIF(C == 1),X_oneMostIF(PDwo ==1)),...
ranksum(X_oneMostIF(C == 1),X_oneMostIF(PDw ==1)),...
ranksum(X_oneMostIF(RBD == 1),X_oneMostIF(PDwo ==1)),...
ranksum(X_oneMostIF(RBD == 1),X_oneMostIF(PDw ==1)),...
ranksum(X_oneMostIF(PDwo == 1),X_oneMostIF(PDw ==1))]









