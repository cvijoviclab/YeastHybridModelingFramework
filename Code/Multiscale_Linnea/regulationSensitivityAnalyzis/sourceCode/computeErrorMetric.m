function [err_metric,errDist,rCoeff] = computeErrorMetric(data,predictions,cutoff)
%computeErrorMetric
%
% data          (double) vector with a distribution of experimental values
%               to use as a reference.
% predictions   (double) vector with a distribution of predicted values
% cutoff        (double) cutoff value to use for substituting zeros
%
% err_metric    (double) mean value of the absolute log10 transformed
%               ratios distribution. This is a measure of the average
%               deviation (in orders of magnitude) of predictions from
%               experimental values.
% errDist       (double) vector with the log10 transformed ratios between
%               predicted and experimental values. Values of zero indicate
%               a 0% of error in the prediction values higher than 1
%               indicate overpredictions by more than one order of
%               magnitude. Values lower than -1 indicate underpredictions
%               by more than one order of magnitude.
% rCoeff        Pearson correlation coefficient between experimental and
%               predicted values,the closer to one, the better the
%               predictions are.
% 
% Usage: [err_metric,errDist,rCoeff] = computeErrorMetric(data,predictions,cutoff)
% 
% Last modified: Ivan Domenzain 2020/06/29
%

if nargin<3
	cutoff = 1E-15;
end
%In case that vector lengths differ, shorten length to the one of the
%shortest one.
data_L = length(data);
pred_L = length(predictions);
minLength   = min(data_L,pred_L);
data        = data(1:minLength);
predictions = predictions(1:minLength);
%Remove zeros or extremely low values
data(data<cutoff) = cutoff;
predictions(predictions<cutoff) = cutoff;
%Get ratio between predicted values and experimental data
relErr     = (predictions./data);
errDist    = log10(relErr);
err_metric = mean(abs(errDist));
%Get pearson correlation coefficient between the two distributions
rCoeff = corrcoef(abs((predictions)),abs((data)));
rCoeff = rCoeff(1,2);
end