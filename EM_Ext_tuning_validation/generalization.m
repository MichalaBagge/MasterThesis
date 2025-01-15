% generalization
clear; clc;
close all

% Define directories
dir_epoch_factor = '...Workspace_Final_Validation\epoch_factor_new';
dir_corr_factor = '...Workspace_Final_Validation\corr_factor';

opt_factor = 1;
names = {'epoch_factor','corr_factor'};
% Set directory
eval(sprintf('cd(dir_%s)',char(names(opt_factor))))


for k = 1:12

    if opt_factor == 1
        eval(sprintf('load("epoch_factor_config_fold%d.mat")',k))

        eval(sprintf('EM_metrics_REMs_testfolds.k%d = EM_metrics_REMs_test;',k))
        eval(sprintf('EM_metrics_SEMs_testfolds.k%d = EM_metrics_SEMs_test;',k))

    elseif opt_factor == 2
        eval(sprintf('load("corr_factor_config_fold%d.mat")',k)) 
    elseif opt_factor == 3
        eval(sprintf('load("Pth_factor_config_fold%d.mat")',k)) 
    end




    eval(sprintf('EM_metrics_EM_testfolds.k%d = EM_metrics_EM_test;',k))
    eval(sprintf('EM_metrics_EMs_testfolds.k%d = EM_metrics_EMs_test;',k))
    eval(sprintf('EMs_metrics_factor_train.k%d = EMs_metrics_factor;',k))
    eval(sprintf('idx_th_train.k%d = idx_th;',k))



    clear EM_metrics_EM_test EM_metrics_EMs_test ...
          EMs_metrics_factor idx_th ...
          EM_metrics_REMs_test EM_metrics_SEMs_test
        

end

% Clear temporary signals for next subject
clearvars -except ...
 EM_metrics_EM_testfolds EM_metrics_EMs_testfolds ...
 EM_metrics_REMs_testfolds EM_metrics_SEMs_testfolds ...
 EMs_metrics_factor_train idx_th_train ...
 epoch_factor fs corr_factor opt_factor names

%%

EMs_metrics_factor = [mean([EMs_metrics_factor_train.k1(:,1), EMs_metrics_factor_train.k2(:,1), ...
                          EMs_metrics_factor_train.k3(:,1), EMs_metrics_factor_train.k4(:,1), ...
                          EMs_metrics_factor_train.k5(:,1), EMs_metrics_factor_train.k6(:,1), ...
                          EMs_metrics_factor_train.k7(:,1), EMs_metrics_factor_train.k8(:,1), ...
                          EMs_metrics_factor_train.k9(:,1), EMs_metrics_factor_train.k10(:,1), ...
                          EMs_metrics_factor_train.k11(:,1), EMs_metrics_factor_train.k12(:,1)],2), ...
                          mean([EMs_metrics_factor_train.k1(:,2), EMs_metrics_factor_train.k2(:,2), ...
                          EMs_metrics_factor_train.k3(:,2), EMs_metrics_factor_train.k4(:,2), ...
                          EMs_metrics_factor_train.k5(:,2), EMs_metrics_factor_train.k6(:,2), ...
                          EMs_metrics_factor_train.k7(:,2), EMs_metrics_factor_train.k8(:,2), ...
                          EMs_metrics_factor_train.k9(:,2), EMs_metrics_factor_train.k10(:,2), ...
                          EMs_metrics_factor_train.k11(:,2), EMs_metrics_factor_train.k12(:,2)],2), ...
                          mean([EMs_metrics_factor_train.k1(:,3), EMs_metrics_factor_train.k2(:,3), ...
                          EMs_metrics_factor_train.k3(:,3), EMs_metrics_factor_train.k4(:,3), ...
                          EMs_metrics_factor_train.k5(:,3), EMs_metrics_factor_train.k6(:,3), ...
                          EMs_metrics_factor_train.k7(:,3), EMs_metrics_factor_train.k8(:,3), ...
                          EMs_metrics_factor_train.k9(:,3), EMs_metrics_factor_train.k10(:,3), ...
                          EMs_metrics_factor_train.k11(:,3), EMs_metrics_factor_train.k12(:,3)],2), ...
                          mean([EMs_metrics_factor_train.k1(:,4), EMs_metrics_factor_train.k2(:,4), ...
                          EMs_metrics_factor_train.k3(:,4), EMs_metrics_factor_train.k4(:,4), ...
                          EMs_metrics_factor_train.k5(:,4), EMs_metrics_factor_train.k6(:,4), ...
                          EMs_metrics_factor_train.k7(:,4), EMs_metrics_factor_train.k8(:,4), ...
                          EMs_metrics_factor_train.k9(:,4), EMs_metrics_factor_train.k10(:,4), ...
                          EMs_metrics_factor_train.k11(:,4), EMs_metrics_factor_train.k12(:,4)],2), ...
                          mean([EMs_metrics_factor_train.k1(:,5), EMs_metrics_factor_train.k2(:,5), ...
                          EMs_metrics_factor_train.k3(:,5), EMs_metrics_factor_train.k4(:,5), ...
                          EMs_metrics_factor_train.k5(:,5), EMs_metrics_factor_train.k6(:,5), ...
                          EMs_metrics_factor_train.k7(:,5), EMs_metrics_factor_train.k8(:,5), ...
                          EMs_metrics_factor_train.k9(:,5), EMs_metrics_factor_train.k10(:,5), ...
                          EMs_metrics_factor_train.k11(:,5), EMs_metrics_factor_train.k12(:,5)],2), ...
                          mean([EMs_metrics_factor_train.k1(:,6), EMs_metrics_factor_train.k2(:,6), ...
                          EMs_metrics_factor_train.k3(:,6), EMs_metrics_factor_train.k4(:,6), ...
                          EMs_metrics_factor_train.k5(:,6), EMs_metrics_factor_train.k6(:,6), ...
                          EMs_metrics_factor_train.k7(:,6), EMs_metrics_factor_train.k8(:,6), ...
                          EMs_metrics_factor_train.k9(:,6), EMs_metrics_factor_train.k10(:,6), ...
                          EMs_metrics_factor_train.k11(:,6), EMs_metrics_factor_train.k12(:,6)],2)];


EMs_metrics_factor_std = [std([EMs_metrics_factor_train.k1(:,1), EMs_metrics_factor_train.k2(:,1), ...
                          EMs_metrics_factor_train.k3(:,1), EMs_metrics_factor_train.k4(:,1), ...
                          EMs_metrics_factor_train.k5(:,1), EMs_metrics_factor_train.k6(:,1), ...
                          EMs_metrics_factor_train.k7(:,1), EMs_metrics_factor_train.k8(:,1), ...
                          EMs_metrics_factor_train.k9(:,1), EMs_metrics_factor_train.k10(:,1), ...
                          EMs_metrics_factor_train.k11(:,1), EMs_metrics_factor_train.k12(:,1)],0,2), ...
                          std([EMs_metrics_factor_train.k1(:,2), EMs_metrics_factor_train.k2(:,2), ...
                          EMs_metrics_factor_train.k3(:,2), EMs_metrics_factor_train.k4(:,2), ...
                          EMs_metrics_factor_train.k5(:,2), EMs_metrics_factor_train.k6(:,2), ...
                          EMs_metrics_factor_train.k7(:,2), EMs_metrics_factor_train.k8(:,2), ...
                          EMs_metrics_factor_train.k9(:,2), EMs_metrics_factor_train.k10(:,2), ...
                          EMs_metrics_factor_train.k11(:,2), EMs_metrics_factor_train.k12(:,2)],0,2), ...
                          std([EMs_metrics_factor_train.k1(:,3), EMs_metrics_factor_train.k2(:,3), ...
                          EMs_metrics_factor_train.k3(:,3), EMs_metrics_factor_train.k4(:,3), ...
                          EMs_metrics_factor_train.k5(:,3), EMs_metrics_factor_train.k6(:,3), ...
                          EMs_metrics_factor_train.k7(:,3), EMs_metrics_factor_train.k8(:,3), ...
                          EMs_metrics_factor_train.k9(:,3), EMs_metrics_factor_train.k10(:,3), ...
                          EMs_metrics_factor_train.k11(:,3), EMs_metrics_factor_train.k12(:,3)],0,2), ...
                          std([EMs_metrics_factor_train.k1(:,4), EMs_metrics_factor_train.k2(:,4), ...
                          EMs_metrics_factor_train.k3(:,4), EMs_metrics_factor_train.k4(:,4), ...
                          EMs_metrics_factor_train.k5(:,4), EMs_metrics_factor_train.k6(:,4), ...
                          EMs_metrics_factor_train.k7(:,4), EMs_metrics_factor_train.k8(:,4), ...
                          EMs_metrics_factor_train.k9(:,4), EMs_metrics_factor_train.k10(:,4), ...
                          EMs_metrics_factor_train.k11(:,4), EMs_metrics_factor_train.k12(:,4)],0,2), ...
                          std([EMs_metrics_factor_train.k1(:,5), EMs_metrics_factor_train.k2(:,5), ...
                          EMs_metrics_factor_train.k3(:,5), EMs_metrics_factor_train.k4(:,5), ...
                          EMs_metrics_factor_train.k5(:,5), EMs_metrics_factor_train.k6(:,5), ...
                          EMs_metrics_factor_train.k7(:,5), EMs_metrics_factor_train.k8(:,5), ...
                          EMs_metrics_factor_train.k9(:,5), EMs_metrics_factor_train.k10(:,5), ...
                          EMs_metrics_factor_train.k11(:,5), EMs_metrics_factor_train.k12(:,5)],0,2), ...
                          std([EMs_metrics_factor_train.k1(:,6), EMs_metrics_factor_train.k2(:,6), ...
                          EMs_metrics_factor_train.k3(:,6), EMs_metrics_factor_train.k4(:,6), ...
                          EMs_metrics_factor_train.k5(:,6), EMs_metrics_factor_train.k6(:,6), ...
                          EMs_metrics_factor_train.k7(:,6), EMs_metrics_factor_train.k8(:,6), ...
                          EMs_metrics_factor_train.k9(:,6), EMs_metrics_factor_train.k10(:,6), ...
                          EMs_metrics_factor_train.k11(:,6), EMs_metrics_factor_train.k12(:,6)],0,2)];


%%
%%% %%% Plot %%% %%%
figure;
subplot(1,5,1)
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,1),'LineWidth', 1.2); 
hold on 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,1) + EMs_metrics_factor_std(:,1) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,1) - EMs_metrics_factor_std(:,1) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
hold off
ax = gca; 
ax.FontSize = 12;
ylabel('Accuracy Score','fontsize',14);
%xlabel('epoch\_factor','fontsize',14);
%xlabel('corr\_factor','fontsize',14);
%xlabel('Pth','fontsize',14);

%%%
subplot(1,5,2)
%fill(EMs_metrics_factor(:,4) + EMs_metrics_factor_std(:,4), EMs_metrics_factor(:,4) - EMs_metrics_factor_std(:,4), 'Color',[0.8, 0.8894, 0.9482])
%hold on
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,4),'LineWidth', 1.2); 
%hold off
 hold on 
 plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,4) + EMs_metrics_factor_std(:,4) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
 plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,4) - EMs_metrics_factor_std(:,4) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
 hold off
ax = gca; 
ax.FontSize = 12;
ylabel('Sensitivity')
%xlabel('epoch\_factor','fontsize',14);
%xlabel('corr\_factor','fontsize',14);
%xlabel('Pth','fontsize',14);

%%%
subplot(1,5,3)
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,5),'LineWidth', 1.2); 
hold on 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,5) + EMs_metrics_factor_std(:,5) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,5) - EMs_metrics_factor_std(:,5) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
hold off
ax = gca; 
ax.FontSize = 12;
ylabel('Specificity Score')
xlabel('epoch\_factor','fontsize',14);
%xlabel('corr\_factor','fontsize',14);
%xlabel('Pth','fontsize',14);

%%%
subplot(1,5,4)
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,2),'LineWidth', 1.2); 
hold on 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,2) + EMs_metrics_factor_std(:,2) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,2) - EMs_metrics_factor_std(:,2) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
hold off
ax = gca; 
ax.FontSize = 12;
ylabel('F1 Score')
%xlabel('epoch\_factor','fontsize',14);
%xlabel('corr\_factor','fontsize',14);
%xlabel('Pth','fontsize',14);

%%%
subplot(1,5,5)
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,6),'LineWidth', 1.2); 
hold on 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,6) + EMs_metrics_factor_std(:,6) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
plot(eval(char(names(opt_factor))),EMs_metrics_factor(:,6) - EMs_metrics_factor_std(:,6) ,'LineWidth', 1.2, 'Color',[0.8, 0.8894, 0.9482]); 
hold off
ax = gca; 
ax.FontSize = 12;
ylabel('Balanced Accuracy')
%xlabel('epoch\_factor','fontsize',14);
%xlabel('corr\_factor','fontsize',14);
%xlabel('Pth','fontsize',14);




%%
all_EMs_test_metrics = [EM_metrics_EMs_testfolds.k1; EM_metrics_EMs_testfolds.k10; ...
 EM_metrics_EMs_testfolds.k11; EM_metrics_EMs_testfolds.k12;...
 EM_metrics_EMs_testfolds.k2; EM_metrics_EMs_testfolds.k3; ...
 EM_metrics_EMs_testfolds.k4; EM_metrics_EMs_testfolds.k5; ...
 EM_metrics_EMs_testfolds.k6; EM_metrics_EMs_testfolds.k7; ...
 EM_metrics_EMs_testfolds.k8; EM_metrics_EMs_testfolds.k9];

all_REMs_test_metrics = [EM_metrics_REMs_testfolds.k1; EM_metrics_REMs_testfolds.k10; ...
 EM_metrics_REMs_testfolds.k11; EM_metrics_REMs_testfolds.k12;...
 EM_metrics_REMs_testfolds.k2; EM_metrics_REMs_testfolds.k3; ...
 EM_metrics_REMs_testfolds.k4; EM_metrics_REMs_testfolds.k5; ...
 EM_metrics_REMs_testfolds.k6; EM_metrics_REMs_testfolds.k7; ...
 EM_metrics_REMs_testfolds.k8; EM_metrics_REMs_testfolds.k9];

all_SEMs_test_metrics = [EM_metrics_SEMs_testfolds.k1; EM_metrics_SEMs_testfolds.k10; ...
 EM_metrics_SEMs_testfolds.k11; EM_metrics_EMs_testfolds.k12;...
 EM_metrics_EMs_testfolds.k2; EM_metrics_EMs_testfolds.k3; ...
 EM_metrics_EMs_testfolds.k4; EM_metrics_EMs_testfolds.k5; ...
 EM_metrics_EMs_testfolds.k6; EM_metrics_EMs_testfolds.k7; ...
 EM_metrics_EMs_testfolds.k8; EM_metrics_EMs_testfolds.k9];

all_EM_test_metrics= [EM_metrics_EM_testfolds.k1; EM_metrics_EM_testfolds.k10; ...
 EM_metrics_EM_testfolds.k11; EM_metrics_EM_testfolds.k12;...
 EM_metrics_EM_testfolds.k2; EM_metrics_EM_testfolds.k3; ...
 EM_metrics_EM_testfolds.k4; EM_metrics_EM_testfolds.k5; ...
 EM_metrics_EM_testfolds.k6; EM_metrics_EM_testfolds.k7; ...
 EM_metrics_EM_testfolds.k8; EM_metrics_EM_testfolds.k9];

% Mean

EM_test_final = mean(all_EM_test_metrics)
EMs_test_final = [mean(all_EMs_test_metrics); std(all_EMs_test_metrics)]

REMs_test_final = [mean(all_REMs_test_metrics); std(all_REMs_test_metrics)]
SEMs_test_final = [mean(all_SEMs_test_metrics); std(all_SEMs_test_metrics)]

%[min(all_EMs_test_metrics); max(all_EMs_test_metrics)]
%diff([min(all_EMs_test_metrics); max(all_EMs_test_metrics)])

