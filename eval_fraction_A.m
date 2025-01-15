function fraction = eval_fraction_A(X, factor_mode, factor,light_annotations)

% Define directories
dir_scripts 
dir_data_RBDPDcontrols
cd(dir_scripts)

% Preallocate for speed
value = zeros(length(X),1);
fraction = table(X,value);
    

% Looping over subjects
for ii =1:length(X)
    dir_subject = char(fullfile(dir_data_RBDPDcontrols,X(ii)));
    cd(dir_subject)
    
        % Set directory to specific subject folder
        dir_subject_RBDPDcontrols = fullfile(dir_data_RBDPDcontrols,char(X(ii)));
        cd(dir_subject_RBDPDcontrols)
        
        if exist('mat-Matteo') || exist('mat-Lykke')
            % Check Matteo and Lykke folder for data extracts
            dir_cleaned_M =fullfile(dir_subject_RBDPDcontrols,'mat-Matteo');
            dir_cleaned_L =fullfile(dir_subject_RBDPDcontrols,'mat-Lykke');
            cd(dir_cleaned_M); M = dir(dir_cleaned_M); M = M(3:end);   
            cd(dir_cleaned_L); L = dir(dir_cleaned_L); L = L(3:end);
        else
            continue
        end
    
        tic
        if ~isempty(M) 
        
        % Set directory
        cd(dir_cleaned_M) 


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%% LOAD POLYSOMNOGRAM SIGNALS %%%%%%%%%%%%%%%%
        % Load data
        subjectID = sprintf('%s.mat',char(X(ii)));
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

        % Clear temporary variables
        clearvars -except eog_r eog_l fs factor_mode factor ...
            dir_scripts dir_subject dir_data_RBDPDcontrols dir_subject_RBDPDcontrols ...
            light_annotations X hypnogram_ext ii fraction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
    % Set directory
    cd(dir_scripts)

    % EM function
     EM_vec = EM_detector_A(eog_r, eog_l, fs, factor_mode, factor);                                                      


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%   PERFORMANCE MEASURE      %%%%%%%%%%%%%%%%

    % Set directory
    cd(dir_scripts)

    % Validation/tuning measures for EM detection
    [em_W, em_R, em_N1, em_N2, em_N3] = PerSleepStage(hypnogram_ext, EM_vec);
    EM_SleepStage = {'em_W', 'em_R', 'em_N1', 'em_N2', 'em_N3'};

    for n = 1:5
        eval(sprintf('percent_%s = 100 * sum(eval(char(EM_SleepStage(n)))) / numel(eval(char(EM_SleepStage(n))));',char(EM_SleepStage(n))));
    end 
    
    % Investigate em percents for NaNs and correct
    if sum(isnan([percent_em_W, percent_em_R, percent_em_N1, percent_em_N3])) >= 1
        temp = [percent_em_W, percent_em_R, percent_em_N1, percent_em_N3];
        em_idx = find(isnan(temp));
        temp(em_idx) = 10^(-3);

        fraction.value(ii) = (temp(1) + temp(2) + temp(3)) / temp(4);


    else
    
        fraction.value(ii) = (percent_em_W + percent_em_R + percent_em_N1) / percent_em_N3;

    end


    
    elseif ~isempty(L) 
        
    % Set directory
    cd(dir_cleaned_L) 

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
        disp(X(ii))
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

    % Adjust all signal length from hypnogram
    light_idx = find(string(light_annotations.Var1) == string(X(ii)));
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

    % Clear temporary variables
    clearvars -except eog_r eog_l fs factor_mode factor ...
        dir_scripts dir_subject dir_data_RBDPDcontrols dir_subject_RBDPDcontrols ...
        light_annotations X hypnogram_ext ii fraction

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%         EM DETECTOR        %%%%%%%%%%%%%%%%
    % Set directory
    cd(dir_scripts)

    % EM function
    EM_vec = EM_detector_A(eog_r, eog_l, fs, factor_mode, factor);                                                      


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%   PERFORMANCE MEASURE      %%%%%%%%%%%%%%%%

    % Set directory
    cd(dir_scripts)

    % Validation/tuning measures for EM detection
    [em_W, em_R, em_N1, em_N2, em_N3] = PerSleepStage(hypnogram_ext, EM_vec);
    EM_SleepStage = {'em_W', 'em_R', 'em_N1', 'em_N2', 'em_N3'};

    for n = 1:5
        eval(sprintf('percent_%s = 100 * sum(eval(char(EM_SleepStage(n)))) / numel(eval(char(EM_SleepStage(n))));',char(EM_SleepStage(n))));
    end 


    % Investigate em percents for NaNs and correct
    if sum(isnan([percent_em_W, percent_em_R, percent_em_N1, percent_em_N3])) >= 1
        temp = [percent_em_W, percent_em_R, percent_em_N1, percent_em_N3];
        em_idx = find(isnan(temp));
        temp(em_idx) = 10^(-3);

        fraction.value(ii) = (temp(1) + temp(2) + temp(3)) / temp(4);
    else
        fraction.value(ii) = (percent_em_W + percent_em_R + percent_em_N1) / percent_em_N3;
    end
        
        
    end

    toc
end
