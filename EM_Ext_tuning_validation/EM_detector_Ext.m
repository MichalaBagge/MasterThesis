function [EM_epoch_vec, EM_vec, EMs_vec, REMs_vec, SEMs_vec, EMs_dur, ...         
    eog_r_recon, eog_l_recon, eog_r_mov, eog_l_mov] ...
    = EM_detector_Ext(eog_r, eog_l, fs)
% This function is eye movement (EM) detection algorithm


%%%%%%% Cleaning procedure %%%%%%%
    % Highpass Butterworth IIR filter
    fc = 0.1;
    [b,a] = butter(4,fc/(fs/2),'high');
    % figure;
    % freqz(b,a)
    eog_r_HP=filtfilt(b,a,eog_r);
    eog_l_HP=filtfilt(b,a,eog_l);
    
    % Extraxt signal parameters
    N = length(eog_r);
    
    % Calculate to nearst multiple of 2^(14) for the WT filter bank
    q = (floor(N/(2^(14))) + 1)*2^(14);
    
    % Create zeros for padding
    zeropad = zeros(q-N,1);
    
    % Zero-padded signal
    eog_r_pad = [eog_r_HP; zeropad];
    eog_l_pad = [eog_l_HP; zeropad];
    
    % EM detector: DTCWT
    J = 14;
    df = dtfilters('dtf3');
    wt_r = dddtree('cplxdt',eog_r_pad,J,df{1},df{2});
    wt_l = dddtree('cplxdt',eog_l_pad,J,df{1},df{2});
    %
    J_eli =  [1, 2, 3, 4, 5, 11, 12, 13, 14];
    for j = 1:length(J_eli)
       wt_r.cfs{J_eli(j)}(:,:,1) = 0;
       wt_r.cfs{J_eli(j)}(:,:,2) = 0;
       wt_l.cfs{J_eli(j)}(:,:,1) = 0;
       wt_l.cfs{J_eli(j)}(:,:,2) = 0;
    end
    
    T_angle = 0.90*pi;                                                       % 
    J = 1:14; J(J_eli)=[];
    for j = 1:length(J)
        R = atan2(wt_r.cfs{J(j)}(:,:,2),wt_r.cfs{J(j)}(:,:,1));
        L = atan2(wt_l.cfs{J(j)}(:,:,2),wt_l.cfs{J(j)}(:,:,1));
        D_angle = abs(R - L);
        D_angle =  min(D_angle, 2*pi -D_angle);
        angle_eli = (D_angle < T_angle);
        wt_r.cfs{J(j)}(angle_eli,:,:) = 0;
        wt_l.cfs{J(j)}(angle_eli,:,:) = 0;
    end
    
    eog_r_recon = idddtree(wt_r); 
    eog_l_recon = idddtree(wt_l); 
    eog_r_recon = eog_r_recon(1:end-length(zeropad));
    eog_l_recon = eog_l_recon(1:end-length(zeropad));
    eog_r_recon = eog_r_recon - mean(eog_r_recon);
    eog_l_recon = eog_l_recon - mean(eog_l_recon);

    %%%%%%% Thresholding Procedure %%%%%%%
    
    % Define absolute difference
    A_diff = abs(eog_r_recon - eog_l_recon);
    
    % Create EM cand
    EM_cand = zeros(size(A_diff));
    T_amp =  600; % 600 muV into microVolt                                  % Tuning
    EM_cand(A_diff > T_amp) = 0;
    % Define the fitted threshold
    Pth = 96;                                                               % Tuning
    th = prctile(A_diff, Pth);     
    EM_cand(A_diff > th) = 1;
    
    % Define parameters
    dur_hole = 1 * fs;                                                      % Tuning
    dur_EM = 2.5 * fs;      

    
    % Initiliaze
    EM_cand1 = EM_cand;
    
    % Finding transitions in the series using column vectors for diff operation
    transitions = diff([0; EM_cand1; 0]);
    starts = find(transitions == 1);
    ends = find(transitions == -1) - 1;
    
    % Close short gaps by finding gaps
    for i = 1:length(starts)-1
        if starts(i+1) - ends(i) <= dur_hole + 1
            EM_cand1((ends(i)+1):(starts(i+1)-1)) = 1;
        end
    end
    
    % Initiliaze
    EM_vec = EM_cand1;

    % Recalculate transitions after closing short gaps
    transitions = diff([0; EM_vec; 0]);
    starts = find(transitions == 1);
    ends = find(transitions == -1) - 1;
    
    % Eliminate short bursts following long gaps by setting them to 0
    for i = 1:(length(starts)-1)
        if ends(i) > starts(i)
            gap =  ends(i) - starts(i);
            
            % Check if the gap is larger than the elimination threshold and the burst is short
            if gap <  dur_EM                         
                EM_vec(starts(i):ends(i)) = 0;  
            end
            
        end
    end

    %%%%%% Extension %%%%%%
    % Initialize
    EM_epoch_vec = EM_vec;

    % Define parameters
    epoch_factor = 4;                                                      % Tuning
    epochs_mini = 1:fs*epoch_factor*30:length(eog_r); 
    
    %%%%%%% Define epochs for to search for individual EM detection %%%%%%%
    for jj = 1:length(epochs_mini)-1
        if sum(EM_epoch_vec(epochs_mini(jj):epochs_mini(jj+1)) == 1) > 10
            EM_epoch_vec(epochs_mini(jj):epochs_mini(jj+1)) = 1;
        else 
            EM_epoch_vec(epochs_mini(jj):epochs_mini(jj+1)) = 0;
        end
    end 

    % Use moving average signal to maintain signal curvature
    window = round(fs/4);
    eog_r_mov = movmean(eog_r, window);  
    eog_l_mov = movmean(eog_l, window);  

    % Identify individual EMs by determination of intersections
    % Find intersections
    intersect_idx = find(diff(sign(eog_r_mov - eog_l_mov))) + 1; 
    intersect_idx = intersect_idx([true; diff(intersect_idx) >= 50]);
    
    % Analyze phase behavior between valid crossovers
    EMs_temp = [];
    
    for k = 1:(length(intersect_idx)-1)
        segment_start = intersect_idx(k);
        segment_end = intersect_idx(k+1);
        segment_eog_r_mov = eog_r_mov(segment_start:segment_end);
        segment_eog_l_mov = eog_l_mov(segment_start:segment_end);
       

        % Calculate correlation to determine
        correlation = corr(segment_eog_r_mov(:), segment_eog_l_mov(:));
        
        % Threshold for correlation factor
        corr_factor = -0.3;  % Threshold for out-of-phase behavior          % Tuning     
%%%     eval(sprintf('%s = %d;',factor_mode, factor));
%%%    
    
        if correlation < corr_factor
            EMs_temp = [EMs_temp; segment_start segment_end];
        end
    end

    % Identify REM and SEM by AASM definitions 
    % Calculate EM duration
    % Determine duration
    EMs_dur = diff(EMs,1,2)/fs;
    EMs_dur = EMs_dur(EMs_dur<=3);
    % Determine EM categories
    REMs = EMs(EMs_dur<=0.5,:); %*2  
    SEMs = EMs(EMs_dur>0.5, :); %*2
     
    % Make REMs, SEMs and EMs of logical characteristics
    EMs_vec = zeros(size(eog_r_mov));
    for n = 1:size(EMs, 1)
        if all(EM_epoch_vec(EMs(n, 1):EMs(n, 2)) == 1) 
            EMs_vec(EMs(n, 1):EMs(n, 2)) = 1;  
        end
    end
    
    REMs_vec = zeros(size(eog_r_mov));
    for n = 1:length(REMs)
        REMs_vec(REMs(n,1):REMs(n,2)) = 1;
    end
    SEMs_vec = zeros(size(eog_r_mov));
    for n = 1:length(SEMs)
        SEMs_vec(SEMs(n,1):SEMs(n,2)) = 1;
    end


end


