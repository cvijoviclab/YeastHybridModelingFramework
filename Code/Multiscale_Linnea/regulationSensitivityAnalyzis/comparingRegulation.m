% comparingRegulation

addpath ('sourceCode');

% Load all data
load('../../../models/reduced_ecYeast_fermentation.mat');
[fermentationDataset1,respirationDataset1,rellativeAbundance, enzUsageSimulated,pUD_1Simulated, geneExprG0, geneExprG1]=loadData();
% Create input data table
modelProteomics = geneExprG0(:,[1 3]);
modelProteomics.G0regulation=geneExprG0.rank;
modelProteomics.G1regulation=geneExprG1.rank;

modelProteomics.indexFermentation = zeros(size(modelProteomics,1),1);
modelProteomics.indexRespiration = zeros(size(modelProteomics,1),1);
regulatedGenes=geneExprG1.enzGenes(geneExprG0.rank~=0 | geneExprG1.rank~=0 );
regulatedGenes(103)={'unregulated'};
for i=1:size(modelProteomics,1)
    j = (find(ismember(fermentationDataset1.ProteinId, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        modelProteomics.indexFermentation(i) = j;
    end
    j = (find(ismember(respirationDataset1.GeneID, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        modelProteomics.indexRespiration(i) = j;
    end
end
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration~=0,:);
respiration.proteomics = respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexRespiration~=0,:));
% Convert molecule/pgDW to mmol/gDW
avogadrosConstant=6.02214076*10^23;
respiration.proteomics=(respiration.proteomics/avogadrosConstant)*10^(12)*10^3;
respiration.proteomics(isnan(respiration.proteomics))=0;
for i=1:length(pUD_1Simulated)
    respiration=[respiration prepRespirationData(modelProteomics,pUD_1Simulated(i))];
    respiration.Properties.VariableNames{'simulated'} = regulatedGenes{i} ;
end
regulationAnalysisRespiration=table();
regulationAnalysisRespiration.Genes=regulatedGenes;
regulationAnalysisRespiration.err_metric=zeros(length(regulatedGenes),1);
regulationAnalysisRespiration.errDist=zeros(length(regulatedGenes),length(respiration.proteomics));
regulationAnalysisRespiration.rCoeff=zeros(length(regulatedGenes),1);
for i=1:length(regulatedGenes)
    [regulationAnalysisRespiration.err_metric(i),regulationAnalysisRespiration.errDist(i,:),regulationAnalysisRespiration.rCoeff(i)] = computeErrorMetric(respiration.proteomics,respiration{:,i+2});
end

%Fermentation
fermentation=table();
fermentation.enzNames = modelProteomics.enzNames(modelProteomics.indexFermentation~=0,:);
fermentation.proteomics = fermentationDataset1.glucose(modelProteomics.indexFermentation(modelProteomics.indexFermentation~=0,:));
fermentation.proteomics(isnan(fermentation.proteomics))=0;
% Remake to percentages
fermentation.proteomics=fermentation.proteomics/sum(fermentationDataset1.glucose(~isnan(fermentationDataset1.glucose)));

for i=1:length(enzUsageSimulated)
    fermentation=[fermentation prepFermentationData(modelProteomics,enzUsageSimulated(i),rellativeAbundance, geneExprG1)];
    fermentation.Properties.VariableNames{'simulated'} = regulatedGenes{i} ;
end

regulationAnalysisFermentation=table();
regulationAnalysisFermentation.Genes=regulatedGenes;
regulationAnalysisFermentation.err_metric=zeros(length(regulatedGenes),1);
regulationAnalysisFermentation.errDist=zeros(length(regulatedGenes),length(fermentation.proteomics));
regulationAnalysisFermentation.rCoeff=zeros(length(regulatedGenes),1);


for i=1:length(regulatedGenes)
    [regulationAnalysisFermentation.err_metric(i),regulationAnalysisFermentation.errDist(i,:),regulationAnalysisFermentation.rCoeff(i)] = computeErrorMetric(fermentation.proteomics,fermentation{:,i+2});
end