%% Main
% Restart session
close all; clc;
clearvars -except subjects

% Define directories
dir_scripts = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Code\MATLAB\EM_tuning_scripts_RH';
dir_data_RBDPDcontrols = 'C:\Users\micha\OneDrive - Danmarks Tekniske Universitet\Semester\Data\RBD PD';
cd(dir_data_RBDPDcontrols);

%% Load subjects 

% Control subjects
X = {'RBD-PD-098', 'RBD-PD-099-b', 'RBD-PD-100', 'RBD-PD-101', 'RBD-PD-102', ...
      'RBD-PD-103', 'RBD-PD-104', 'RBD-PD-105', 'RBD-PD-106', 'RBD-PD-107', ...
     'RBD-PD-109-b', 'RBD-PD-110', 'RBD-PD-111', 'RBD-PD-112', 'RBD-PD-113', ...
     'RBD-PD-114', 'RBD-PD-115', 'RBD-PD-116', 'RBD-PD-117', 'RBD-PD-118', ...
     'RBD-PD-119', 'RBD-PD-120', 'RBD-PD-121', 'RBD-PD-122', 'RBD-PD-123', ...
     'RBD-PD-124', 'RBD-PD-125', 'RBD-PD-126', 'RBD-PD-127', 'RBD-PD-128', ...
     'RBD-PD-129', 'RBD-PD-130', 'RBD-PD-131', 'RBD-PD-132', 'RBD-PD-133', ...
     'RBD-PD-134', 'RBD-PD-135', 'RBD-PD-166', 'RBD-PD-168', 'RBD-PD-169', ... 
     'RBD-PD-173', 'RBD-PD-174', 'RBD-PD-175', 'RBD-PD-176', 'RBD-PD-177', ...
     'RBD-PD-178', 'RBD-PD-179', 'RBD-PD-183', 'RBD-PD-184', 'RBD-PD-186', ...
     'RBD-PD-189', 'RBD-PD-191', 'RBD-PD-192', 'RBD-PD-193', 'RBD-PD-194', ...
     'RBD-PD-195', 'RBD-PD-196', 'RBD-PD-197', 'RBD-PD-198', 'RBD-PD-199', ...
     'RBD-PD-201', 'RBD-PD-202', 'RBD-PD-203', 'RBD-PD-204', 'RBD-PD-207', ...
     'RBD-PD-208', 'RBD-PD-209', 'RBD-PD-210', 'RBD-PD-211', 'RBD-PD-212', ...
     'RBD-PD-213', 'RBD-PD-214', 'RBD-PD-215'};
X = X';

% Define range for adaptive threshold factors
% Pth_factor = [90, 92, 94, 96, 98];
% T_angle_factor = [0.50*pi, 0.55*pi, 0.60*pi, 0.65*pi, 0.70*pi, 0.75*pi, 0.80*pi, 0.85*pi, 0.90*pi, 0.95*pi];
dur_hole_factor = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0];
dur_EM_factor = [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50];

% Preallocate for speed
% Pth = zeros(length(X),length(Pth_factor));
% T_angle = zeros(length(X),length(T_angle_factor));
dur_hole = zeros(length(X),length(dur_hole_factor));
dur_EM = zeros(length(X),length(dur_EM_factor));
%fraction_OptParameters = table(X,Pth, T_angle, dur_hole, dur_EM);
fraction_OptParameters = table(X, dur_hole, dur_EM);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     LIGHT ANNOTATIONS      %%%%%%%%%%%%%%%%
% Set light off/on annotations
cd(dir_data_RBDPDcontrols);
light_annotations = cell2table(readcell('Light epochs.txt'));
light_annotations = light_annotations(375:end,1:3);
light_annotations(15,2:3) = table({250},{2000}); 
light_annotations(26,2:3) = table({950}, {1866}); 
light_annotations(148,2:3) = table({0}, {2000}); 
logical_idx_names = ismember(table2cell(light_annotations(:,1)), X);
light_annotations = light_annotations(logical_idx_names,:);
for n = 1:length(table2cell(light_annotations))
    logical_idx_doubles(n) = isnumeric(light_annotations.Var2{n})&& isnumeric(light_annotations.Var3{n}) &&...
                             ~isnan(light_annotations.Var2{n})&& ~isnan(light_annotations.Var3{n}); 
end
light_annotations = light_annotations(logical_idx_doubles,:);
clearvars logical_idx_doubles logical_idx_names n   



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%         TUNING LOOP        %%%%%%%%%%%%%%%%

% Parameter tuning 
% %%%%% Pth %%%%% 
% for jj = 1:length(Pth_factor)
%     % Set directory
%     cd(dir_scripts)
%     % Print CV fold of 
%     fprintf('Parameter optimization (Pth) %d/%d\n', jj, length(Pth_factor));
%     % Define parameter for optimization
%     factor_mode = 'Pth';
%     % Calculate fraction
%     fraction = eval_fraction(X, factor_mode, Pth_factor(jj),light_annotations);
%     % Determine optimal parameter of fraction 
%     fraction_OptParameters.Pth(:,jj) = fraction.value;
%     % Clear temporary variabels
%     clear fraction
% end 
% % Clear temporary variables for speed
% clear jj factor_mode
% %%%%% T_angle %%%%% 
% for jj = 1:length(T_angle_factor)
%     % Set directory
%     cd(dir_scripts)
%     % Print CV fold of 
%     fprintf('Parameter optimization (T_angle) %d/%d\n', jj, length(T_angle_factor));
%     % Define parameter for optimization
%     factor_mode = 'T_angle';
%     % Calculate fraction
%     fraction = eval_fraction(X, factor_mode, T_angle_factor(jj),light_annotations);
%     % Determine optimal parameter of fraction 
%     fraction_OptParameters.T_angle(:,jj) = fraction.value;
%     % Clear temporary variabels
%     clear fraction
% end 
% Clear temporary variables for speed
% clear jj factor_mode
% Parameter tuning 

%%%%% dur_hole %%%%% 
% Set directory
cd(dir_scripts)
for jj = 1:length(dur_hole_factor)
    % Set directory
    cd(dir_scripts)
    % Print CV fold of 
    fprintf('Parameter optimization (dur_hole) %d/%d\n', jj, length(dur_hole_factor));
    % Define parameter for optimization
    factor_mode = 'dur_hole';
    % Calculate fraction
    fraction = eval_fraction_A(X, factor_mode, dur_hole_factor(jj),light_annotations);
    % Determine optimal parameter of fraction 
    fraction_OptParameters.dur_hole(:,jj) = fraction.value;
    % Clear temporary variabels
    clear fraction
end 
% Clear temporary variables for speed
clear jj factor_mode
%%%%% dur_EM %%%%% 
for jj = 1:length(dur_EM_factor)
    % Set directory
    cd(dir_scripts)
    % Print CV fold of 
    fprintf('Parameter optimization (dur_EM) %d/%d\n', jj, length(dur_EM_factor));
    % Define parameter for optimization
    factor_mode = 'dur_EM';
    % Calculate fraction
    fraction = eval_fraction_A(X, factor_mode, dur_EM_factor(jj),light_annotations);
    % Determine optimal parameter of fraction 
    fraction_OptParameters.dur_EM(:,jj) = fraction.value;
    % Clear temporary variabels
    clear fraction
end 
% Clear temporary variables for speed
clear jj factor_mode



%% Save parameter tuning
% Set directory
cd(dir_scripts)
% Save feature table to csv-file 
writetable(fraction_OptParameters, 'Fraction_Optimization_Parameters.csv');  



