function [fig_roc, opt_th, best_fpr, best_tpr] = generat_RF_EM_Part2_plots(FPR, TPR,T,K)

    ksize = [size(FPR.k1,1),size(FPR.k2,1),size(FPR.k3,1)...
          size(FPR.k4,1), size(FPR.k5,1)];
    idx_ksize = find(ksize == min(ksize));
    kFPRTPR_size = ksize(idx_ksize);
    
    for k = 1:K
        eval(sprintf('[FPR_int_%d, ~, idx_fpr] = unique(FPR.k%d);',k,k));
        eval(sprintf('tpr_temp = accumarray(idx_fpr, TPR.k%d, [], @mean);',k));
        eval(sprintf('TPR_int_%d = interp1(FPR_int_%d,tpr_temp,linspace(0,1,1000));',k,k));
        clear tpr_temp idx_fpr]

        eval(sprintf('T_%d = interp1(linspace(0,1,length(T.k%d)),T.k%d,linspace(0,1,1000));',k,k,k));

       
    end
    
    
    FPR_ave = linspace(0,1,1000);
    TPR_ave = mean([TPR_int_1; TPR_int_2; TPR_int_3; TPR_int_4; TPR_int_5],1);
    T_ave = mean([T_1; T_2; T_3; T_4; T_5]);
    %
    % Compute Euclidean distance to (0, 1)
    distances = sqrt((FPR_ave).^2 + (1 - TPR_ave).^2);
    [~, idx] = min(distances); 
    best_fpr = FPR_ave(idx); 
    best_tpr = TPR_ave(idx); 
    opt_th = T_ave(idx);
    AUC = trapz(FPR_ave, TPR_ave);
    
    fig_roc = figure;
    hold on
    plot(linspace(0,1,1000), TPR_int_1, 'LineWidth', 1.2, 'Color', [0.878, 0.884, 0.921]) 
    plot(linspace(0,1,1000), TPR_int_2, 'LineWidth', 1.2, 'Color', [0.878, 0.884, 0.921]) 
    plot(linspace(0,1,1000), TPR_int_3, 'LineWidth', 1.2, 'Color', [0.878, 0.884, 0.921]) 
    plot(linspace(0,1,1000), TPR_int_4, 'LineWidth', 1.2, 'Color', [0.878, 0.884, 0.921]) 
    plot(linspace(0,1,1000), TPR_int_5, 'LineWidth', 1.2, 'Color', [0.878, 0.884, 0.921]) 
    plot(FPR_ave, TPR_ave, 'LineWidth', 1.5, 'Color', [0.267 0.306 0.525])% 
    title(sprintf('Random Forest Model (AUC = %.2f)', AUC));
    plot([0 1], [0 1], 'k--', 'DisplayName', 'Random Guess');
    plot(best_fpr, best_tpr, 'ro')
    hold off
    xlabel('False Positive Rate');
    ylabel('True Positive Rate');
    
    

end
