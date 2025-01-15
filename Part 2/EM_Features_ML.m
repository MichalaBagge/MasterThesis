%% Features and other ML methods
% Restart session
close all; clc;
clearvars -except subjects

% Define directories
%...
cd(dir_scripts)
 
% Load subjects 
F=dir(dir_data);
F=F(3:end); % FOR USE
% F=F(13:end); % FOR TEST

% Initiate table
 EM_subject_Feature_Table = table();

% Subject loop
for ii =1:length(F)
% Load subject     
dir_subject = fullfile(dir_data,F(ii).name);
cd(dir_subject)

%%%  Load PSG data and light off/on annotations
load('eogl-m2.mat'); load('eogr-m2.mat'); load('vec_hypnogram.mat');  
eog_l = eoglm2; eog_r = eogrm2;
clear eogrm2 eoglm2

% Create extended hypnogram
hypnogram_ext = reshape(transpose(repmat(hypnogram, 1, 256*30)), [], 1);

%%% Adjust all signal length from hypnogram and eog
% Load light off / on targets 
cd(dir_annotations)
% EM stop annotation
stop = sprintf('%s_Fabio_StopEMscoring_targetvector.mat',char(F(ii).name));
load(stop);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
stop_temp = find(targetVector_R == 1);
stop_idx = stop_temp(end);
clear targetVector_R stop
% EM start annotation
start = sprintf('%s_Fabio_StartEMscoring_targetvector.mat',char(F(ii).name));
load(start);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
start_temp = find(targetVector_R == 1);
start_idx = start_temp(1);
clear targetVector_R start

% Create subject fields for eog and hypnogram extended in subject
% struct
eog_l = eog_l(start_idx:stop_idx);
eog_r = eog_r(start_idx:stop_idx);
hypnogram_ext = hypnogram_ext(start_idx:stop_idx);

%%% Load EM annotationer
% Set directory
cd(dir_annotations)

% EM REM annotation
rem = sprintf('%s_Fabio_REM_targetvector.mat', char(F(ii).name));
load(rem);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
REM_annotations= targetVector_R(start_idx:stop_idx);
clear targetVector_R rem
% EM SEM annotation
sem = sprintf('%s_Fabio_SEM_targetvector.mat', char(F(ii).name));
load(sem);
targetVector_R = double(targetVector_R.x); clear targetVector_R.x;
SEM_annotations = targetVector_R(start_idx:stop_idx);
clear targetVector_R sem
% EM annotation vector
EM_annotations = SEM_annotations + REM_annotations;
EM_annotations(EM_annotations==2) = 1;

%%% Create subject time axis and epochs
fs = 256;
t = (0:length(eog_r)-1) / fs;
epochs = 1:fs*30:length(eog_r);
epochs_mini = 1:fs*3*60:length(eog_r);

cd(dir_scripts)
% EM function
[EM_epoch_vec, EM_vec, EMs, eog_r_recon_chri, eog_l_recon_chri, eog_r_recon, eog_l_recon] ...
    = EM_detector_AllOptions(eog_r, eog_l, fs); 


% Prealloacte Features 
EM_ID = transpose(1:length(EMs));
EM_Duration = zeros(length(EMs),1);
EM_corr_coef = zeros(length(EMs),1);
EM_CF = zeros(length(EMs),1);
EM_Power = zeros(length(EMs),1);
EM_Slope = zeros(length(EMs),1);
EM_deriv_Entropy = zeros(length(EMs),1);
No_EM_Extrema = zeros(length(EMs),1);
y = zeros(length(EMs),1);
Subejct_ID = repmat(char(F(ii).name),length(EMs),1);

% Preallocate EM Feature Table
EM_Feature_Table = table(Subejct_ID, EM_ID, EM_Duration, EM_corr_coef, EM_CF, EM_Power, ...
    EM_deriv_Entropy, No_EM_Extrema, EM_Slope, y);

% Make count
count = 0;

     %%%% Extract Features Per Eye Movement (EM) %%%%
    for n =1:length(EMs)
    
        %%%% Time Series Features %%%%
        % Estimate EM duration
        EM_Feature_Table.EM_Duration(n) = (EMs(n,2) - EMs(n,1))/fs;
        % Estimate correlation coefficient
        EM_Feature_Table.EM_corr_coef(n) = corr(eog_r_recon(EMs(n,1):EMs(n,2)), eog_l_recon(EMs(n,1):EMs(n,2)));
        % Estimate the power
        EM_Feature_Table.EM_Power(n) = mean(eog_r_recon(EMs(n,1):EMs(n,2)).^2);
        %%%%  Estimate First Derivative Curvature Features %%%% 
        t_segment = t(EMs(n,1):15:EMs(n,2)); % Time for segment
        eog_segment = eog_r_recon(EMs(n,1):15:EMs(n,2))'; % Signal segment
        slope_segment = diff(eog_segment) ./ diff(t_segment); % First derivative
        slope_time = t_segment(1:end-1); 

        % First derivative - Entropy feature
        [N,edges,bin] = histcounts(abs(slope_segment)); 
        P = N / sum(N); % Normalize to probabilities
        EM_Feature_Table.EM_deriv_Entropy(n) = -sum(P(P>0) .* log2(P(P>0)));

        % First derivative - Detect curvature shifts
        critical_idxs = find(diff(sign(slope_segment)) ~= 0); % Locate sign changes
        critical_t = slope_time(critical_idxs+1);%%%%%%%%%%
        critical_segs = eog_segment(critical_idxs+1);
        if isempty(diff(critical_idxs))
            critical_seg = eog_segment(critical_idxs+1);
        else
        count = count +1;
            critical_idx = critical_idxs((diff(critical_idxs) > 5));
            critical_idx = [critical_idx+1, critical_idxs(end)+1];
            critical_seg = eog_segment(critical_idx);
        end
        EM_Feature_Table.No_EM_Extrema(n) = length(critical_seg);

        % First derivative - Initiate slope feature
        
        if isempty(critical_idxs)
            EM_Feature_Table.EM_Slope(n) = NaN;
        elseif isempty(diff(critical_idxs))            
            x1 = t_segment(1);
            y1 = eog_segment(1);
            x2 = critical_t;
            y2 = critical_segs;
            EM_Feature_Table.EM_Slope(n) = (y2 - y1) / (x2 - x1);
        else
            critical_time =slope_time(critical_idx); %
            x1 = t_segment(1);
            y1 = eog_segment(1);
            x2 = critical_time(1);
            y2 = critical_seg(1);
            EM_Feature_Table.EM_Slope(n) = (y2 - y1) / (x2 - x1);
        end


        
        %%%% Spectral Features %%%%
        % Perform Welch Method
        [S_welch, f_welch] = pwelch(eog_r_recon(EMs(n,1):EMs(n,2)),[], [], [], fs);
        % Calculate central frequency
        EM_Feature_Table.EM_CF(n) = sum(S_welch .* f_welch) / sum(S_welch);
       

        %%%%%%%% Decide label - REM/SEM type %%%%%%%%
        if sum(REM_annotations(EMs(n,1):EMs(n,2))) == 0 && sum(SEM_annotations(EMs(n,1):EMs(n,2))) == 0
            EM_Feature_Table.y(n) = 0;
        else
            %if sum(REM_annotations(EMs(n,1):EMs(n,2))) > length(REM_annotations(EMs(n,1):EMs(n,2)))*0.5
            if sum(REM_annotations(EMs(n,1):EMs(n,2))) >= 1 || ...
               sum(SEM_annotations(EMs(n,1):EMs(n,2))) >= 1     
                if sum(REM_annotations(EMs(n,1):EMs(n,2))) >= ...
                   sum(SEM_annotations(EMs(n,1):EMs(n,2))) 
    
                    EM_Feature_Table.y(n) = 1;
    
                elseif sum(REM_annotations(EMs(n,1):EMs(n,2))) < ...
                   sum(SEM_annotations(EMs(n,1):EMs(n,2))) 
                    
                    EM_Feature_Table.y(n) = 2;
                    
                elseif  sum(REM_annotations(EMs(n,1):EMs(n,2))) >= 1
         
                    EM_Feature_Table.y(n) = 1;
    
                else
                    
                    EM_Feature_Table.y(n) = 2;
    
                end
            end
        end


        
    end

        % Save length of EM_ID to get length for subject_ID
        EM_subject_Feature_Table = [EM_subject_Feature_Table; EM_Feature_Table];



        % Clear temporary variables
        clearvars -except dir_scripts dir_annotations dir_data dir_subject ...
                   F ii EM_subject_Feature_Table

end 



%%%% Save feature table %%%%
% Set directory
cd(dir_scripts)
% Save feature table to csv-file 
writetable(EM_subject_Feature_Table, "AllSubjects_REM_SEM_Feature_Table.csv");








