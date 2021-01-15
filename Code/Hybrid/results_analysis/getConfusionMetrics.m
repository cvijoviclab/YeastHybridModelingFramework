function [confMat,metrics] = getConfusionMetrics(observations,predictions)
TP = sum(observations>0 & predictions>0);
FP = sum(observations==0 & predictions>0);
TN = sum(observations==0 & predictions==0);
FN = sum(observations>0 & predictions==0);
%arrange results as confusion matrix
obs_positives = [TP;FN];
obs_negatives = [FP;TN];
confMat = [obs_positives,obs_negatives];
%compute metrics for predictive performance assesment
sensitivity = TP/(TP+FN);
specificity = TN/(TN+FP);
precision   = TP/(TP+FP);
accuracy    = (TP+TN)/(TP+TN+FP+FN);
F1_score    = 2*TP/(2*TP+FP+FN);
MCC         = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
term1        = TP/(TP+FP);
term2        = TP/(TP+FN);
FMI          = sqrt(term1*term2);
informedness = sensitivity+specificity-1; 
NPV          = TN/(TN+FN);
markedness   = precision+NPV-1; 
metrics = [sensitivity;specificity;precision;accuracy;F1_score;MCC;FMI;informedness;markedness];
end