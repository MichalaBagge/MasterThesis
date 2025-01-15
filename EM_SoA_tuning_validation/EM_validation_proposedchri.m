%% Main
% Restart session
close all; clc;
clearvars -except subjects

% Define directories
dir_data = '';
dir_annotations = '';
dir_scripts = '';
cd(dir_scripts)
 
% Load subjects 
F=dir(dir_data);
F=F(3:end); % FOR USE
% F=F(13:end); % FOR TEST

% Preallocate for speed
EM_metrics_EMs_vec = zeros(length(F),6);                                    % OBS!
EM_metrics_EM_vec = zeros(length(F),6); 
EM_metrics_EM_vec_SoA = zeros(length(F),6); 


for ii =1:length(F)
    tic

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% LOAD POLYSOMNOGRAM SIGNALS %%%%%%%%%%%%%%%%
   
    % Load polysomnogram data
    [eog_l, eog_r, hypnogram_ext] = load_polysomnogram(dir_data,char(F(ii).name));
    % Load light off / on targets 
    cd(dir_scripts)
    [~,~,~,start_idx,stop_idx] = load_annotations(char(F(ii).name));
    eog_l = eog_l(start_idx:stop_idx);
    eog_r = eog_r(start_idx:stop_idx);
    hypnogram_ext = hypnogram_ext(start_idx:stop_idx);

    % Create subject time axis and epochs
    fs = 256;
    t = (0:length(eog_r)-1) / fs;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
    % Set directory
    cd(dir_scripts)
    [~, EM_vec, EMs_vec, REMs_vec, SEMs_vec, EMs_dur, ~, ~] ...
        = EM_detector_Ext(eog_r, eog_l, fs);

    EM_vec_SoA = EM_detector_SoA(eog_r, eog_l, fs);



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%     LOAD EM ANNOTATIONS    %%%%%%%%%%%%%%%%
    [EM_annotations, REM_annotations, SEM_annotations,start_idx,stop_idx] = load_annotations(char(F(ii).name));

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%   PERFORMANCE MEASURE      %%%%%%%%%%%%%%%%

    % Set directory
    cd(dir_scripts)
    % Performance measures for EM detection
    EM_metrics_EM_vec_SoA(ii, :) = EM_performance_val(EM_vec_SoA, EM_annotations);   % OBS!
    EM_metrics_EM_vec(ii, :) = EM_performance_val(EM_vec, EM_annotations);   % OBS!
    EM_metrics_EMs_vec(ii, :) = EM_performance_val(EMs_vec, EM_annotations); % OBS!

    % Clear temporary signals for next subject
    clearvars -except dir_scripts dir_data dir_annotations ...
    F ii ...
    EM_metrics_EM_vec EM_metrics_EMs_vec EM_metrics_EM_vec_SoA ...
    EM_metrics_EMs_vec EM_metrics_REMs_vec EM_metrics_SEMs_vec

    toc

end    

%%

A0 = [mean(EM_metrics_EM_vec_SoA(:, 1)), ...
 mean(EM_metrics_EM_vec_SoA(:, 2)), ...
 mean(EM_metrics_EM_vec_SoA(:, 3)), ...
 mean(EM_metrics_EM_vec_SoA(:, 4)), ...
 mean(EM_metrics_EM_vec_SoA(:, 5)), ...
 mean(EM_metrics_EM_vec_SoA(:, 6)); ... %
 std(EM_metrics_EM_vec_SoA(:, 1)), ...
 std(EM_metrics_EM_vec_SoA(:, 2)), ...
 std(EM_metrics_EM_vec_SoA(:, 3)), ...
 std(EM_metrics_EM_vec_SoA(:, 4)), ...
 std(EM_metrics_EM_vec_SoA(:, 5)), ...
 std(EM_metrics_EM_vec_SoA(:, 6))];
    
A = [mean(EM_metrics_EM_vec(:, 1)), ...
 mean(EM_metrics_EM_vec(:, 2)), ...
 mean(EM_metrics_EM_vec(:, 3)), ...
 mean(EM_metrics_EM_vec(:, 4)), ...
 mean(EM_metrics_EM_vec(:, 5)), ...
 mean(EM_metrics_EM_vec(:, 6)); ... %
 std(EM_metrics_EM_vec(:, 1)), ...
 std(EM_metrics_EM_vec(:, 2)), ...
 std(EM_metrics_EM_vec(:, 3)), ...
 std(EM_metrics_EM_vec(:, 4)), ...
 std(EM_metrics_EM_vec(:, 5)), ...
 std(EM_metrics_EM_vec(:, 6))];


C = [mean(EM_metrics_EMs_vec(:, 1)), ...
 mean(EM_metrics_EMs_vec(:, 2)), ...
 mean(EM_metrics_EMs_vec(:, 3)), ...
 mean(EM_metrics_EMs_vec(:, 4)), ...
 mean(EM_metrics_EMs_vec(:, 5)), ...
 mean(EM_metrics_EMs_vec(:, 6)); ... %
 std(EM_metrics_EMs_vec(:, 1)), ...
 std(EM_metrics_EMs_vec(:, 2)), ...
 std(EM_metrics_EMs_vec(:, 3)), ...
 std(EM_metrics_EMs_vec(:, 4)), ...
 std(EM_metrics_EMs_vec(:, 5)), ...
 std(EM_metrics_EMs_vec(:, 6))];

