function [respiration]=prepRespirationData(modelProteomics,pU )
%Respiration
respiration=table();
respiration.simulated = table2array(pU{1}(modelProteomics.indexRespiration~=0,5));
% Correct for satturation
respiration.simulated = respiration.simulated/0.38;