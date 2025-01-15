function EM_metrics = EM_performance_val(EM_vec, EM_annotations)

% Determine performance measures for EM detection
TP = sum((EM_vec == 1) & (EM_annotations == 1));                       
TN = sum((EM_vec == 0) & (EM_annotations == 0));                       
FP = sum((EM_vec == 1) & (EM_annotations == 0));                       
FN = sum((EM_vec == 0) & (EM_annotations == 1));                          

Precision = TP / (TP + FP);
Sensitivity = TP / (TP + FN); % Sensitivity / Recall
Specificity = TN / (TN + FP); 
Accuracy = (TP + TN) / (TP + TN + FP + FN);
Accuracy_Balanced = (Sensitivity + Specificity) / 2;
F1 = 2 * (Precision * Sensitivity) / (Precision + Sensitivity);
EM_metrics = [Accuracy, F1, Precision, Sensitivity, Specificity, Accuracy_Balanced];


