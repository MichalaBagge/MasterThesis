% Main
% Restart session
close all; clear; clc;

% Define directories
dir_scripts = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Code\MATLAB\EM_tuning_scripts_RH';


%% Load subjects 

% Load data
% file_path = fullfile(dir_scripts,'Fraction_Optimization_Parameters_C.csv');
%file_path = fullfile(dir_scripts,'Fraction_Optimization_Parameters_A0.csv');
file_path = fullfile(dir_scripts,'Fraction_Optimization_Parameters_SoA.csv');
% file_path = fullfile(dir_scripts,'Fraction_Optimization_Parameters.csv');
data_table = readtable(file_path);
attributeNames = data_table.Properties.VariableNames;
attributeNames = attributeNames(2:end);
X =data_table.X;


% Fraction for evaluation
for n = 1:length(attributeNames)
    eval(sprintf('idx_inf = data_table.%s < Inf;',char(attributeNames(n))));   
    eval(sprintf('evaluation_fraction_train(n) = mean(data_table.%s(idx));',char(attributeNames(n)))); 
    eval(sprintf('evaluation_fraction_std(n) = std(data_table.%s(idx));',char(attributeNames(n)))); 
end



%% Plot
% Define range for adaptive threshold factors
Pth_factor = [90, 92, 94, 96, 98];
T_angle_factor = [0.50*pi, 0.55*pi, 0.60*pi, 0.65*pi, 0.70*pi, 0.75*pi, 0.80*pi, 0.85*pi, 0.90*pi, 0.95*pi];
dur_hole_factor = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
dur_EM_factor = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50];
epoch_factor = 2:2:30;  
corr_factor = -0.1:-0.2:-1; 
%%%%%%%%%%%%%%%%%%%%%%%

% %%% Plot for A0 Model %%%
% Display result 
figure;
plot(Pth_factor, evaluation_fraction_train(:,1:5),'.-','LineWidth', 1.2);  
hold on
plot(Pth_factor, evaluation_fraction_train(:,1:5)-evaluation_fraction_std(:,1:5),'-', 'Color',[0.5, 0.7235, 0.8705]);   
plot(Pth_factor, evaluation_fraction_train(:,1:5)+evaluation_fraction_std(:,1:5),'-', 'Color',[0.5, 0.7235, 0.8705]);  
hold off
%legend({'Training','Test'},'Location','NorthWest','fontsize',14);
ax = gca; % Get current axis
ax.FontSize = 14;
xlabel('Pth','fontsize',14);
ylabel('Fraction','fontsize',14);
xlim([90,98])

figure;
plot(T_angle_factor, evaluation_fraction_train(:,6:15),'.-','LineWidth', 1.2);   
hold on
plot(T_angle_factor,evaluation_fraction_train(:,6:15)-evaluation_fraction_std(:,6:15),'-', 'Color',[0.5, 0.7235, 0.8705]);   
plot(T_angle_factor, evaluation_fraction_train(:,6:15)+evaluation_fraction_std(:,6:15),'-', 'Color',[0.5, 0.7235, 0.8705]);  
hold off
%legend({'Training','Test'},'Location','NorthWest','fontsize',14);
ax = gca; % Get current axis
ax.FontSize = 14;
xlabel('T\_angle [rad]','fontsize',14);
xticklabels({'0.50\pi', '0.55\pi', '0.60\pi', '0.65\pi', '0.70\pi','0.75\pi', '0.80\pi', '0.85\pi', '0.90\pi', '0.95\pi'}); 
xticklabels({'0.50\pi', '0.65\pi', '0.80\pi', '0.95\pi'}); 
ylabel('Fraction','fontsize',14);

%%% Plots for C Model %%%
% figure;
% plot(epoch_factor, mean(evaluation_fraction_train(:,16:30)),'.-');   
% hold on
% plot(epoch_factor, mean(evaluation_fraction_test(:,16:30)),'.-');   
% hold off
% legend({'Training','Test'},'Location','NorthWest','fontsize',14);
% xlabel('epoch','fontsize',14);
% ylabel('Fraction','fontsize',14);
% 
% figure;
% plot(corr_factor, mean(evaluation_fraction_train(:,31:end)),'.-'); 
% hold on
% plot(corr_factor, mean(evaluation_fraction_test(:,31:end)),'.-'); 
% hold off
% legend({'Training','Test'},'Location','NorthWest','fontsize',14);
% xlabel('correlation coefficient','fontsize',14);
% ylabel('Fraction','fontsize',14);


% Plot for dur parametre

figure;
%plot(dur_hole_factor(3:end), evaluation_fraction_train(:,18:25),'.-');   
plot(dur_hole_factor, evaluation_fraction_train(:,16:25),'.-','LineWidth', 1.2);
hold on
plot(dur_hole_factor,evaluation_fraction_train(:,16:25)-evaluation_fraction_std(:,16:25),'-', 'Color',[0.5, 0.7235, 0.8705]);   
plot(dur_hole_factor, evaluation_fraction_train(:,16:25)+evaluation_fraction_std(:,16:25),'-', 'Color',[0.5, 0.7235, 0.8705]);  
hold off
%hold on
%plot(dur_hole_factor, mean(evaluation_fraction_test(:,1:10)),'.-');   
%hold off
%legend({'Training','Test'},'Location','NorthWest','fontsize',14);
ax = gca; % Get current axis
ax.FontSize = 14;
xlabel('dur\_hole [s]','fontsize',14);
ylabel('Fraction','fontsize',14);

figure;
plot(dur_EM_factor, evaluation_fraction_train(:,26:35),'.-','LineWidth', 1.2);   
hold on
plot(dur_EM_factor,evaluation_fraction_train(:,26:35)-evaluation_fraction_std(:, 26:35),'-', 'Color',[0.5, 0.7235, 0.8705]);   
plot(dur_EM_factor, evaluation_fraction_train(:,26:35)+evaluation_fraction_std(:,26:35),'-', 'Color',[0.5, 0.7235, 0.8705]);  
hold off
%legend({'Training','Test'},'Location','NorthWest','fontsize',14);
%set(gcf, 'Color', [0.5, 0.7235, 0.8705]); % Light gray figure background
ax = gca; % Get current axis
ax.FontSize = 14;
%set(ax, 'Color', [0.9500, 0.9724, 0.9870]); % Light grayish blue background for the graph area
xlabel('dur\_EM [s]','fontsize',14);
ylabel('Fraction','fontsize',14);



%% 

Pth_fracall = evaluation_fraction_train(:,1:5);

T_angle_fracall = evaluation_fraction_train(:,6:15);

dur_hole_fracall = evaluation_fraction_train(:,16:25);

dur_EM_fracall = evaluation_fraction_train(:,26:35);

%%

summar = [min(Pth_fracall), max(Pth_fracall); ...
 min(T_angle_fracall), max(T_angle_fracall);...
 min(dur_hole_fracall), max(dur_hole_fracall); ...
 min(dur_EM_fracall), max(dur_EM_fracall)];



