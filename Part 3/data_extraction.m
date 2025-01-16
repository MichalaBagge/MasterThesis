% Data extraction script for controls, RBD and PD patients
clear; clc;

% Define and set directory
cd(dir_data_RBDPDcontrols);
% Get overview of subject folders
F = dir(dir_data_RBDPDcontrols); F = F(3:end);
F = {F.name};
F_idx = isfolder(fullfile(dir_data_RBDPDcontrols,F));
F = F(F_idx);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%     LIGHT ANNOTATIONS      %%%%%%%%%%%%%%%%
% Set light off/on annotations
light_annotations = cell2table(readcell('Light epochs.txt'));
light_annotations = light_annotations(375:end,1:3);
light_annotations(15,2:3) = table({250},{2000});
light_annotations(26,2:3) = table({950}, {1866});
light_annotations(148,2:3) = table({0}, {2000}); 
logical_idx_names = ismember(table2cell(light_annotations(:,1)), F);
light_annotations = light_annotations(logical_idx_names,:);
for n = 1:length(table2cell(light_annotations))
    logical_idx_doubles(n) = isnumeric(light_annotations.Var2{n})&& isnumeric(light_annotations.Var3{n}) &&...
                             ~isnan(light_annotations.Var2{n})&& ~isnan(light_annotations.Var3{n}); 
end
light_annotations = light_annotations(logical_idx_doubles,:);
clearvars logical_idx_doubles logical_idx_names F_idx n   


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        PREALLOCATE         %%%%%%%%%%%%%%%%
% Initialize subject IDs in feature table
subject_ID = F.';
Feature_Table = table(subject_ID, 'VariableNames', {'SubjectID'});

%%%%%%%%%%%%%%%%    Preallocate features    %%%%%%%%%%%%%%%%
% Preallocate EM covarge feature in percents  
Feature_Table.EMcovarge_overall = zeros(size(F.')); % zeros(length(F.'),1);
Feature_Table.EMcovarage_N3 = zeros(size(F.'));
Feature_Table.EMcovarage_N2 = zeros(size(F.'));
Feature_Table.EMcovarage_N1 = zeros(size(F.'));
Feature_Table.EMcovarage_REM = zeros(size(F.'));
Feature_Table.EMcovarage_Wake = zeros(size(F.'));
% Number of all EMs total
Feature_Table.NoEMs = zeros(size(F.'));
Feature_Table.NoREMs = zeros(size(F.'));
Feature_Table.NoSEMs  = zeros(size(F.'));
% Number of all EMs per sleep stages
Feature_Table.NoEMs_N3 = zeros(size(F.'));
Feature_Table.NoEMs_N2 = zeros(size(F.'));
Feature_Table.NoEMs_N1 = zeros(size(F.'));
Feature_Table.NoEMs_REM = zeros(size(F.'));
Feature_Table.NoEMs_Wake = zeros(size(F.'));
Feature_Table.NoREMs_N3 = zeros(size(F.'));
Feature_Table.NoREMs_N2 = zeros(size(F.'));
Feature_Table.NoREMs_N1 = zeros(size(F.'));
Feature_Table.NoREMs_REM = zeros(size(F.'));
Feature_Table.NoREMs_Wake = zeros(size(F.'));
Feature_Table.NoSEMs_N3 = zeros(size(F.'));
Feature_Table.NoSEMs_N2 = zeros(size(F.'));
Feature_Table.NoSEMs_N1 = zeros(size(F.'));
Feature_Table.NoSEMs_REM = zeros(size(F.'));
Feature_Table.NoSEMs_Wake = zeros(size(F.'));
% Average duration of all EMs
Feature_Table.EMAverageDuration = zeros(size(F.'));
Feature_Table.REMAverageDuration = zeros(size(F.'));
Feature_Table.SEMAverageDuration = zeros(size(F.'));
% Remaining Time series features
Feature_Table.EM_corr_coef = zeros(size(F.'));
Feature_Table.EM_deriv_Entropy = zeros(size(F.'));
Feature_Table.No_EM_Extrema = zeros(size(F.'));
Feature_Table.EM_Slope = zeros(size(F.'));
Feature_Table.REM_corr_coef = zeros(size(F.'));
Feature_Table.REM_deriv_Entropy = zeros(size(F.'));
Feature_Table.No_REM_Extrema = zeros(size(F.'));
Feature_Table.REM_Slope = zeros(size(F.'));
Feature_Table.SEM_corr_coef = zeros(size(F.'));
Feature_Table.SEM_deriv_Entropy = zeros(size(F.'));
Feature_Table.No_SEM_Extrema = zeros(size(F.'));
Feature_Table.REM_Slope = zeros(size(F.'));
% % % WT Energy of all EMs
Feature_Table.EM_WT_TotalAverageEnergy = zeros(size(F.'));
Feature_Table.REM_WT_TotalAverageEnergy = zeros(size(F.'));
Feature_Table.SEM_WT_TotalAverageEnergy = zeros(size(F.'));
% Central frequency and power EM features
Feature_Table.EM_CF = zeros(size(F.'));
Feature_Table.EM_Power = zeros(size(F.'));
Feature_Table.REM_CF = zeros(size(F.'));
Feature_Table.REM_Power = zeros(size(F.'));
Feature_Table.SEM_CF = zeros(size(F.'));
Feature_Table.SEM_Power = zeros(size(F.'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%        FEATURE LOOP        %%%%%%%%%%%%%%%%
% Specify subject
for ii = 1:length(F)
    tic

    % Set directory to specific subject folder
    dir_subject_RBDPDcontrols = fullfile(dir_data_RBDPDcontrols,char(F(ii)));
    cd(dir_subject_RBDPDcontrols)
    
    if exist('mat-Matteo')
        % Check Matteo and Lykke folder for data extracts
        dir_cleaned_M =fullfile(dir_subject_RBDPDcontrols,'mat-Matteo');
        cd(dir_cleaned_M); M = dir(dir_cleaned_M); M = M(3:end);  
    elseif exist('mat-Lykke')
        dir_cleaned_L =fullfile(dir_subject_RBDPDcontrols,'mat-Lykke'); 
        cd(dir_cleaned_L); L = dir(dir_cleaned_L); L = L(3:end);
        M = [];
    else
        Feature_Table(ii,:) =[];
        continue
    end

    if ~isempty(M) 
        
        % Set directory
        cd(dir_cleaned_M)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% LOAD POLYSOMNOGRAM SIGNALS %%%%%%%%%%%%%%%%
        % Load data
        subjectID = sprintf('%s.mat',char(F(ii)));
        load(subjectID);

        % Identify EOG channels
        for i = 1:length(channelNames)
            channelNames(i) = strtrim(cellstr(channelNames{i,1}));
        end

        r = find(strcmp(cellstr(channelNames(:,1)),'EOGH'));
        l = find(strcmp(cellstr(channelNames(:,1)),'EOGV'));

        if isempty(r) && isempty(l)
            r = find(strcmp(cellstr(channelNames(:,1)),'EOGH-A1'));
            l = find(strcmp(cellstr(channelNames(:,1)),'EOGV-A2')); 
        end

        % Extract data
        eog_r_extract = PSG_data{r,1};
        eog_l_extract = PSG_data{l,1};
        
        % Convert to convension used in main script
        eog_r = eog_r_extract(:);
        eog_l = eog_l_extract(:);

        % Extract sampling frequency
        fs_initial = samplingRate(r);

        % Investigate hypnogram for NaNs and correct
        if sum(isnan(hyp)) >= 1
            hypnan_idx = find(isnan(hyp));
            % Correct to same epoch as preceding
            hyp(hypnan_idx)=hyp(hypnan_idx(1)-1);    
        end

        % Create extended hypnogram
        hypnogram_ext = reshape(repmat(hyp, fs_initial*30, 1), [], 1);
        hypnogram_ext = round(hypnogram_ext);

        % Adjust inclusion of subjects based on hypnograms 
        TST_temp = (length(eog_r)/fs_initial)/3600;
        if ((sum(hypnogram_ext == 1))/fs_initial)/3600 > TST_temp*0.6
            continue
        end

        % Resample to fs = 256 Hz for EM detection procedure 
        [p,q] = rat(256/fs_initial);
        if fs_initial ~= 256
        eog_r = resample(eog_r,p,q);
        eog_l = resample(eog_l,p,q);
        hypnogram_ext =  resample(hypnogram_ext,p,q);
        end
        hypnogram_ext = round(hypnogram_ext);


        % Create subject time axis and epochs
        fs = 256;
        t = (0:length(eog_r)-1) / fs;
        epochs = 1:fs*30:length(eog_r);
        epochs_mini = 1:fs*5:length(eog_r);
        TST = (length(eog_r)/fs)/3600;



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
        % Set directory
        cd(dir_scripts)
        
        % Use EM detector
        [EM_epoch_vec, EM_vec, EMs_vec, EMs, REMs, SEMs, ...         
         eog_r_recon, eog_l_recon, eog_r_mov, eog_l_mov] ...
        = EM_detector_Ext(eog_r, eog_l, fs);
         
    
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%           FEATURES         %%%%%%%%%%%%%%%%
    
        % Calculate EM coverage features to table
        Feature_Table.EMcovarge_overall(ii) = (sum(EM_vec) / length(EM_vec)) * 100;
        
        % Calculation EM coverage for each stage to table
         covarage_stages = zeros(1,5);
        for n = -3:1
            stage = EM_vec(hypnogram_ext == n);    
            if ~isempty(stage)
                covarage_stages(n+4) = (sum(stage) / length(stage))  * 100;
            else
                covarage_stages(n+4) = NaN; 
            end
        end
        clear stage
        Feature_Table(ii,3:7) = array2table(covarage_stages);
    
    
        % Calculate number of EMs
        Feature_Table.NoEMs(ii) = size(EMs,1);
        Feature_Table.NoREMs(ii) = size(REMs,1);
        Feature_Table.NoSEMs(ii) = size(SEMs,1);

        % Set directories to feature function folder
        cd(dir_feature_funcs)
        
        % Calculate for each sleep stage to the table 
        NoEMsSleepStages = NumberofEMspersSleepStage(EMs,hypnogram_ext);
        Feature_Table(ii,9:13) = array2table(diag(NoEMsSleepStages)');
        NoREMsSleepStages = NumberofEMspersSleepStage(REMs,hypnogram_ext);
        Feature_Table(ii,14:18) = array2table(diag(NoREMsSleepStages)');
        NoSEMsSleepStages = NumberofEMspersSleepStage(SEMs,hypnogram_ext);
        Feature_Table(ii,19:23) = array2table(diag(NoSEMsSleepStages)');
    
        % Calculate average duration per EM (hours per TST)
        Feature_Table.EMAverageDuration(ii) = mean((diff(EMs,1,2)/fs)/3600,1)/TST;
        Feature_Table.REMAverageDuration(ii) = mean((diff(REMs,1,2)/fs)/3600,1)/TST;
        Feature_Table.SEMAverageDuration(ii) = mean((diff(SEMs,1,2)/fs)/3600,1)/TST;
        
        % Remaining Time series features
       [Feature_Table.EM_corr_coef(ii),Feature_Table.EM_deriv_Entropy(ii), ...
        Feature_Table.No_EM_Extrema(ii),Feature_Table.EM_Slope(ii)] = EMTimeSeriesAndSlopes(EMs, eog_r_recon, eog_l_recon);
       [Feature_Table.REM_corr_coef(ii),Feature_Table.REM_deriv_Entropy(ii), ...
        Feature_Table.No_REM_Extrema(ii),Feature_Table.REM_Slope(ii)] = EMTimeSeriesAndSlopes(REMs, eog_r_recon, eog_l_recon);
       [Feature_Table.SEM_corr_coef(ii),Feature_Table.SEM_deriv_Entropy(ii), ...
        Feature_Table.No_SEM_Extrema(ii),Feature_Table.SEM_Slope(ii)] = EMTimeSeriesAndSlopes(SEMs, eog_r_recon, eog_l_recon);
        
        % Calculate EM WT energy features
        [Feature_Table.EM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, EMs);   
        [Feature_Table.REM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, REMs);
        [Feature_Table.SEM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, SEMs);
    
        % Central frequency and power EM features
        % Average CF and powers across all EMs
        [EM_CF_temp, EM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, EMs, fs);
        Feature_Table.EM_CF(ii) = mean(EM_CF_temp);
        Feature_Table.EM_Power(ii) = mean(EM_Power_temp);
    
        [REM_CF_temp, REM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, REMs, fs);
        Feature_Table.REM_CF(ii) = mean(REM_CF_temp);
        Feature_Table.REM_Power(ii) = mean(REM_Power_temp);

        [SEM_CF_temp, SEM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, SEMs, fs);
        Feature_Table.SEM_CF(ii) = mean(SEM_CF_temp);
        Feature_Table.SEM_Power(ii) = mean(SEM_Power_temp);
        
    
        % Set directories to the main script folder
        cd(dir_scripts)   

        % Clear temporary signals for next subject
        clearvars -except dir_scripts dir_data_RBDPDcontrols dir_subject ...
                dir_feature_funcs dir_figures F Feature_Table light_annotations


    elseif  ~isempty(L) 
     


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%% LOAD POLYSOMNOGRAM SIGNALS %%%%%%%%%%%%%%%%
    % Load eog (left,right) and hypnogram
    load('eogl-a2.mat'); eog_l = data;
    load('eogr-a1.mat'); eog_r = data;
    load('hypnogram.mat'); 
    % Verify readme-file to get sampling frequency
    if exist('#ANONYMOUSreadme#.mat')
        load('#ANONYMOUSreadme#.mat'); fs_initial = readme.fs(1);
    elseif exist('#readme#.mat')
        load('#readme#.mat'); fs_initial = readme.fs(1);
    else
        disp(F(ii))
    end
    clear data readme


    % Investigate hypnogram for NaNs and correct
    if sum(isnan(hypnogram)) >= 1
        hypnan_idx = find(isnan(hypnogram));
        % Correct to same epoch as preceding
        hypnogram(hypnan_idx)=hypnogram(hypnan_idx(1)-1);    
    end


    % Create extended hypnogram
    hypnogram_ext = reshape(transpose(repmat(hypnogram, 1, fs_initial*30)), [], 1);
    hypnogram_ext = round(hypnogram_ext);

    % Adjust all signal length from hypnogram
    light_idx = find(string(light_annotations.Var1) == string(F(ii)));
    start_idx =  cell2mat(light_annotations.Var2(light_idx))*fs_initial*30;
    stop_idx =  cell2mat(light_annotations.Var3(light_idx))*fs_initial*30;
    
    % Create subject fields for eog and hypnogram extended in subject
    % struct
    eog_l = eog_l(start_idx:stop_idx);
    eog_r = eog_r(start_idx:stop_idx);
    hypnogram_ext = hypnogram_ext(start_idx:stop_idx);


    % Resample to fs = 256 Hz for EM detection procedure 
    [p,q] = rat(256/fs_initial);
    if fs_initial ~= 256
    eog_r = resample(eog_r,p,q);
    eog_l = resample(eog_l,p,q);
    hypnogram_ext =  resample(hypnogram_ext,p,q);
    end
    hypnogram_ext = round(hypnogram_ext);


    % Create subject time axis and epochs
    fs = 256;
    t = (0:length(eog_r)-1) / fs;
    epochs = 1:fs*30:length(eog_r);
    epochs_mini = 1:fs*5:length(eog_r);
    TST = (length(eog_r)/fs)/3600;


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
    % Set directory
    cd(dir_scripts)
    
    % Use EM detector
      [EM_epoch_vec, EM_vec, EMs_vec, EMs, REMs, SEMs, ...         
      eog_r_recon, eog_l_recon, eog_r_mov, eog_l_mov] ...
       = EM_detector_Ext(eog_r, eog_l, fs);

           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%           FEATURES         %%%%%%%%%%%%%%%%
    
        % Calculate EM coverage features to table
        Feature_Table.EMcovarge_overall(ii) = (sum(EM_vec) / length(EM_vec)) * 100;
        
        % Calculation EM coverage for each stage to table
         covarage_stages = zeros(1,5);
        for n = -3:1
            stage = EM_vec(hypnogram_ext == n);    
            if ~isempty(stage)
                covarage_stages(n+4) = (sum(stage) / length(stage))  * 100;
            else
                covarage_stages(n+4) = NaN; 
            end
        end
        clear stage
        Feature_Table(ii,3:7) = array2table(covarage_stages);
    
    
        % Calculate number of EMs
        Feature_Table.NoEMs(ii) = size(EMs,1);
        Feature_Table.NoREMs(ii) = size(REMs,1);
        Feature_Table.NoSEMs(ii) = size(SEMs,1);

        % Set directories to feature function folder
        cd(dir_feature_funcs)
        
        % Calculate for each sleep stage to the table 
        NoEMsSleepStages = NumberofEMspersSleepStage(EMs,hypnogram_ext);
        Feature_Table(ii,9:13) = array2table(diag(NoEMsSleepStages)');
        NoREMsSleepStages = NumberofEMspersSleepStage(REMs,hypnogram_ext);
        Feature_Table(ii,14:18) = array2table(diag(NoREMsSleepStages)');
        NoSEMsSleepStages = NumberofEMspersSleepStage(SEMs,hypnogram_ext);
        Feature_Table(ii,19:23) = array2table(diag(NoSEMsSleepStages)');
    
        % Calculate average duration per EM (hours per TST)
        Feature_Table.EMAverageDuration(ii) = mean((diff(EMs,1,2)/fs)/3600,1)/TST;
        Feature_Table.REMAverageDuration(ii) = mean((diff(REMs,1,2)/fs)/3600,1)/TST;
        Feature_Table.SEMAverageDuration(ii) = mean((diff(REMs,1,2)/fs)/3600,1)/TST;
        
        % Remaining Time series features
       [Feature_Table.EM_corr_coef(ii),Feature_Table.EM_deriv_Entropy(ii), ...
        Feature_Table.No_EM_Extrema(ii),Feature_Table.EM_Slope(ii)] = EMTimeSeriesAndSlopes(EMs, eog_r_recon, eog_l_recon);
       [Feature_Table.REM_corr_coef(ii),Feature_Table.REM_deriv_Entropy(ii), ...
        Feature_Table.No_REM_Extrema(ii),Feature_Table.REM_Slope(ii)] = EMTimeSeriesAndSlopes(REMs, eog_r_recon, eog_l_recon);
       [Feature_Table.SEM_corr_coef(ii),Feature_Table.SEM_deriv_Entropy(ii), ...
        Feature_Table.No_SEM_Extrema(ii),Feature_Table.SEM_Slope(ii)] = EMTimeSeriesAndSlopes(SEMs, eog_r_recon, eog_l_recon);
        

        % Calculate EM WT energy features
        [Feature_Table.EM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, EMs);   
        [Feature_Table.REM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, REMs);
        [Feature_Table.SEM_WT_TotalAverageEnergy(ii), ~,~] = ...
         EM_WT_Energy_FeatureFunction(eog_r_recon, SEMs);
    
        % Central frequency and power EM features
        % Average CF and powers across all EMs
        [EM_CF_temp, EM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, EMs, fs);
        Feature_Table.EM_CF(ii) = mean(EM_CF_temp);
        Feature_Table.EM_Power(ii) = mean(EM_Power_temp);
    
        [REM_CF_temp, REM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, REMs, fs);
        Feature_Table.REM_CF(ii) = mean(REM_CF_temp);
        Feature_Table.REM_Power(ii) = mean(REM_Power_temp);

        [SEM_CF_temp, SEM_Power_temp] ...
        = EM_Spectral_FeatureFunction(eog_r_recon, SEMs, fs);
        Feature_Table.SEM_CF(ii) = mean(SEM_CF_temp);
        Feature_Table.SEM_Power(ii) = mean(SEM_Power_temp);
        

    % Set directories to the main script folder
    cd(dir_scripts)

    % Clear temporary signals for next subject
    clearvars -except dir_scripts dir_figures dir_data_RBDPDcontrols dir_subject ...
               dir_feature_funcs F Feature_Table light_annotations


    else
        continue
    end
    toc
end


%% Save feature table
% Set directory
cd(dir_scripts)
% Save feature table to csv-file 
writetable(Feature_Table, 'FeatureTable_RBDPD_Part3.csv');  




