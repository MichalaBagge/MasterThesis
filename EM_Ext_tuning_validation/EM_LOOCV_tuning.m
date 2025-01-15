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

% Classification
X = char(F(:).name);
N = length(F);
CV = cvpartition(N, 'Leaveout');

% For each crossvalidation fold
for k = 1:N
    fprintf('Crossvalidation fold %d/%d\n', k, N);
    
    % Extract the training and test set
    train_set = cellstr(X(CV.training(k), :)); 
    test_set = cellstr(X(CV.test(k), :));

% Define
%epoch_factor = 0.05:0.05:2;
epoch_factor = 0.5:0.5:3;

% Preallocate for speed
EM_metrics_EMs_vec = zeros(length(train_set),6);                                    % OBS!
EM_metrics_EM_vec = zeros(length(train_set),6); 
EM_metrics_REMs_vec = zeros(length(train_set),6);  
EM_metrics_SEMs_vec = zeros(length(train_set),6);  
EM_metrics_factor = zeros(length(epoch_factor),6);
EMs_metrics_factor = zeros(length(epoch_factor),6);
REMs_metrics_factor = zeros(length(epoch_factor),6);
SEMs_metrics_factor = zeros(length(epoch_factor),6);


for n = 1:length(epoch_factor)        

        tic
    
    for ii =1:length(train_set)
        factor_mode = 'epoch_factor';
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% LOAD POLYSOMNOGRAM SIGNALS %%%%%%%%%%%%%%%%
       
        % Load polysomnogram data
        [eog_l, eog_r, hypnogram_ext] = load_polysomnogram(dir_data,char(train_set(ii)));
        % Load light off / on targets 
        cd(dir_scripts)
        [~,~,~,start_idx,stop_idx] = load_annotations(char(train_set(ii)));
        eog_l = eog_l(start_idx:stop_idx);
        eog_r = eog_r(start_idx:stop_idx);
        hypnogram_ext = hypnogram_ext(start_idx:stop_idx);
    
        % Create subject time axis and epochs
        fs = 256;
        t = (0:length(eog_r)-1) / fs;
        epochs = 1:fs*30:length(eog_r);
        epochs_mini = 1:fs*epoch_factor(n)*30:length(eog_r);                       % OBS!
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
        % Set directory
        cd(dir_scripts)

        % C2 Use EM detector
%         [EM_epoch_vec, EM_vec, EMs, REMs, SEMs,  EMs_vec, eog_r_recon, ...           % OBS!
%         eog_l_recon, A_diff,th] = EM_detector_val_C2(eog_r, eog_l, ...          % OBS!
%         fs,epochs_mini);                                                        % OBS!  
%     

        [EM_epoch_vec, EM_vec, EMs_vec, REMs_vec, SEMs_vec, ...        
        eog_r_recon, eog_l_recon, eog_r_mov, eog_l_mov]  ...
            = EM_detector_Ext(eog_r, eog_l, fs, factor_mode, epoch_factor(n));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%     LOAD EM ANNOTATIONS    %%%%%%%%%%%%%%%%
        [EM_annotations,REM_annotations, SEM_annotations,start_idx,stop_idx] = load_annotations(char(train_set(ii)));
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%   PERFORMANCE MEASURE      %%%%%%%%%%%%%%%%
    
        % Set directory
        cd(dir_scripts)
        % Performance measures for EM detection
        EM_metrics_EM_vec(ii, :) = EM_performance_val(EM_vec, EM_annotations);   % OBS!
        EM_metrics_EMs_vec(ii, :) = EM_performance_val(EMs_vec, EM_annotations); % OBS!
        EM_metrics_REMs_vec(ii, :) = EM_performance_val(REMs_vec, REM_annotations); % OBS!
        EM_metrics_SEMs_vec(ii, :) = EM_performance_val(SEMs_vec, SEM_annotations); % OBS!

    
        % Clear temporary signals for next subject
        clearvars -except dir_scripts dir_data dir_annotations ... 
        dir_figures dir_feature_funcs F ...
        t hypnogram_ext EM_cand EM_vec EMs_vec EM_annotations ...
        eog_r_recon eog_l_recon A_diff th plot_factor ...
        EM_metrics_EM_vec EM_metrics_EMs_vec ...
        EM_metrics_factor EMs_metrics_factor ...
        epoch_factor n train_set test_set  ...
        factor_mode k N CV X ...
        EM_metrics_REMs_vec EM_metrics_SEMs_vec
                                                          % OBS!
        
    
    
    end    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%     TEST SET METRICES      %%%%%%%%%%%%%%%%
                      
    EM_metrics_factor(n,:) = [mean(EM_metrics_EM_vec(:,1)), mean(EM_metrics_EM_vec(:,2)) ...
                      mean(EM_metrics_EM_vec(:,3)), mean(EM_metrics_EM_vec(:,4)), ...
                      mean(EM_metrics_EM_vec(:,5)), mean(EM_metrics_EM_vec(:,6))];
    
    EMs_metrics_factor(n,:) = [mean(EM_metrics_EMs_vec(:,1)), mean(EM_metrics_EMs_vec(:,2)) ...
                      mean(EM_metrics_EMs_vec(:,3)), mean(EM_metrics_EMs_vec(:,4)), ...
                      mean(EM_metrics_EMs_vec(:,5)), mean(EM_metrics_EMs_vec(:,6))];

    REMs_metrics_factor(n,:) = [mean(EM_metrics_REMs_vec(:,1)), mean(EM_metrics_REMs_vec(:,2)) ...
                      mean(EM_metrics_REMs_vec(:,3)), mean(EM_metrics_REMs_vec(:,4)), ...
                      mean(EM_metrics_REMs_vec(:,5)), mean(EM_metrics_REMs_vec(:,6))];

    SEMs_metrics_factor(n,:) = [mean(EM_metrics_SEMs_vec(:,1)), mean(EM_metrics_SEMs_vec(:,2)) ...
                      mean(EM_metrics_SEMs_vec(:,3)), mean(EM_metrics_SEMs_vec(:,4)), ...
                      mean(EM_metrics_SEMs_vec(:,5)), mean(EM_metrics_SEMs_vec(:,6))];
        toc   
end

        % Clear temporary signals for next subject
        clearvars -except dir_scripts dir_data dir_annotations ... 
        dir_figures dir_feature_funcs F ...
        EM_metrics_EM_vec EM_metrics_EMs_vec ...
        EM_metrics_factor EMs_metrics_factor ...
        epoch_factor n train_set test_set   ...
        factor_mode k N CV X ...
        EM_metrics_REMs_vec EM_metrics_SEMs_vec ...
        REMs_metrics_factor SEMs_metrics_factor
        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%    OPTIMAL FACTOR K-FOLD   %%%%%%%%%%%%%%%%

% Compute parameters
idx_th = [find(EMs_metrics_factor(:,1) == max(EMs_metrics_factor(:,1))); ...
          find(EMs_metrics_factor(:,2) == max(EMs_metrics_factor(:,2))); ...
          find(EMs_metrics_factor(:,3) == max(EMs_metrics_factor(:,3))); ...
          find(EMs_metrics_factor(:,4) == max(EMs_metrics_factor(:,4))); ...
          find(EMs_metrics_factor(:,5) == max(EMs_metrics_factor(:,5))); ...
          find(EMs_metrics_factor(:,6) == max(EMs_metrics_factor(:,6)))];

% Compute parameters
opt_epoch = epoch_factor(idx_th(6));


% Load polysomnogram data
[eog_l, eog_r, hypnogram_ext] = load_polysomnogram(dir_data,char(test_set));
% Load light off / on targets 
cd(dir_scripts)
[~,~,~,start_idx,stop_idx] = load_annotations(char(test_set));
eog_l = eog_l(start_idx:stop_idx);
eog_r = eog_r(start_idx:stop_idx);
hypnogram_ext = hypnogram_ext(start_idx:stop_idx);
% Create subject time axis and epochs
fs = 256;
t = (0:length(eog_r)-1) / fs;
epochs = 1:fs*30:length(eog_r);
epochs_mini = 1:fs*epoch_factor(n)*30:length(eog_r);                       % OBS!


% EM detector
cd(dir_scripts)
[EM_epoch_vec, EM_vec, EMs_vec, REMs_vec, SEMs_vec, ...        
eog_r_recon, eog_l_recon, eog_r_mov, eog_l_mov] ...
   = EM_detector_Ext(eog_r, eog_l, fs, factor_mode, opt_epoch);

% load EM annotations
[EM_annotations, REM_annotations, SEM_annotations] = load_annotations(char(test_set));

% Performance measures for EM detection
cd(dir_scripts)
EM_metrics_EM_test = EM_performance_val(EM_vec, EM_annotations);   % OBS!
EM_metrics_EMs_test = EM_performance_val(EMs_vec, EM_annotations); % OBS!
EM_metrics_REMs_test = EM_performance_val(REMs_vec, REM_annotations);   % OBS!
EM_metrics_SEMs_test = EM_performance_val(SEMs_vec, SEM_annotations); % OBS!



% Visualize per fold
% factor_visualize
% close all
eval(sprintf('save("epoch_factor_config_fold%d.mat")',k))


end


% Clear temporary signals for next subject
clearvars -except dir_annotations dir_data ...
dir_scripts k N CV X 



