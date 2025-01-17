function EM_vec = EM_detector_A(eog_r, eog_l, fs, factor_mode, factor)
% This function is a reconstruction of the eye movement (EM) detection
% algorithm proposed by Christensen et al. (2017).
% Model output: EM vector array for EM coverage percentage
%
% Written by Michala Bagge

%% Define apative parameter for tuning 

% Define the adaptive parameter
eval(sprintf('%s = %d;',factor_mode, factor));

% Initiate adaptive loop
if exist(factor_mode, 'var') % 
%% Cleaning procedure
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
    
    T_angle = 0.90*pi;
    eval(sprintf('%s = %d;',factor_mode, factor));                      
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

    %Thresholding Procedure
    % Define absolute difference
    A_diff = abs(eog_r_recon - eog_l_recon);
    
    % Create EM cand
    EM_cand = zeros(size(A_diff));
    T_amp =  600; % 600 muV into microVolt                                       
    EM_cand(A_diff > T_amp) = 0;
    % Define the fitted threshold
    Pth = 92;
    eval(sprintf('%s = %d;',factor_mode, factor));                               
    th = prctile(A_diff, Pth);       
    EM_cand(A_diff > th) = 1;

    
    % Thresholding
    % Define parameters
    dur_hole = 2 * fs;
    eval(sprintf('%s = %d*fs;',factor_mode, factor));                              
    dur_EM = 2.5 * fs;                                                           
    eval(sprintf('%s = %d*fs;',factor_mode, factor));                            

    
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
    
end
end