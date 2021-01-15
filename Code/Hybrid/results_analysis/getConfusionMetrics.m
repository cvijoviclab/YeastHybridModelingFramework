function [sensitivity,specificity,precision,accuracy,MCC] = getConfusionMetrics(observations,predictions)
TP = sum(observations>0 & predictions>0);
FP = sum(observations==0 & predictions>0);
TN = sum(observations==0 & predictions==0);
FN = sum(observations>0 & predictions==0);

sensitivity = TP/(TP+FN);
specificity = TN/(TN+FP);
precision   = TP/(TP+FP);
accuracy    = (TP+TN)/(TP+TN+FP+FN);
MCC         = ((TP*TN)-(FP*FN))/sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
end