function [fermentation]=prepFermentationData(modelProteomics, enzUsageSimulated,rellativeAbundance, geneExprG1)
fermentation=table();
fermentation.simulated = enzUsageSimulated{1}.D_04(modelProteomics.indexFermentation~=0,:);
% Correct for satturation
fermentation.simulated = fermentation.simulated/0.38;
% make to percentages
modelProteomics.indexRellativeAbundance = zeros(size(modelProteomics,1),1);
for i=1:size(modelProteomics,1)
    j = (find(ismember(rellativeAbundance.string_external_id,['4932.' char(geneExprG1.enzGenes(i))])));
    if ~isempty(j)
        modelProteomics.indexRellativeAbundance(i) = j;
    end
end
modelProteomics.RellativeAbundance = zeros(size(modelProteomics,1),1);
for i=1:size(modelProteomics,1)
    if modelProteomics.indexRellativeAbundance(i)~=0
        modelProteomics.RellativeAbundance(i) = table2array(rellativeAbundance(modelProteomics.indexRellativeAbundance(i),3));
    end
end
%%rellative abundance is in ppm /10000-->%
modelProteomics.RellativeAbundance=modelProteomics.RellativeAbundance/10000;
fermentation.simulated=fermentation.simulated/sum(modelProteomics.RellativeAbundance);
