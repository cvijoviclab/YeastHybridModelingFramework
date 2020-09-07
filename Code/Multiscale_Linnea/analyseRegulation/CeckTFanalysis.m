close all
clear all
%add path to functions
addpath ('../sourceCode');
addpath ('functions');

[rellativeAbundance,proteomics, D,Dregulated, enzUsages_reducedYeast,enzUsages_reducedYeastRegulated, geneExprG0, geneExprG1]=loadData();
[fermentationDataset1,respirationDataset1]=importNewProteomics;

unregulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529/enzUsage_VarAnalysis_D_1.txt');
dataset1005D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr2/enzUsage_VarAnalysis_D_1.txt');
dataset1005D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr2/enzUsage_VarAnalysis_D_4.txt');

dataset2005D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr3/enzUsage_VarAnalysis_D_1.txt');
dataset2010D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr4/enzUsage_VarAnalysis_D_1.txt');
dataset2050D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr5/enzUsage_VarAnalysis_D_1.txt');
dataset2005D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr3/enzUsage_VarAnalysis_D_4.txt');
dataset2010D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr4/enzUsage_VarAnalysis_D_4.txt');
dataset2050D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr5/enzUsage_VarAnalysis_D_4.txt');

dataset3005D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr6/enzUsage_VarAnalysis_D_1.txt');
dataset3010D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr7/enzUsage_VarAnalysis_D_1.txt');
dataset3005D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr6/enzUsage_VarAnalysis_D_4.txt');
dataset3010D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr7/enzUsage_VarAnalysis_D_4.txt');

datasetJoint005D1=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr8/enzUsage_VarAnalysis_D_1.txt');
datasetJoint005D4=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/results/checkTF20200529nr8/enzUsage_VarAnalysis_D_4.txt');



modelProteomics = geneExprG0(:,[1 3]);
modelProteomics.indexFermentation = zeros(size(modelProteomics,1),1);
modelProteomics.indexRespiration = zeros(size(modelProteomics,1),1);
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

%Respiration
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration~=0,:);
respiration.proteomics = respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexRespiration~=0,:));
% Convert molecule/pgDW to mmol/gDW
avogadrosConstant=6.02214076*10^23;
respiration.proteomics=(respiration.proteomics/avogadrosConstant)*10^(12)*10^3;
% Correct for satturation
respiration.unregulated=unregulated{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset1005D1=dataset1005D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset2005D1=dataset1005D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset3005D1=dataset1005D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.datasetJoint005D1=dataset1005D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.proteomics(isnan(respiration.proteomics)) =0;
respiration.unregulated(isnan(respiration.unregulated)) =0;
respiration.dataset1005D1(isnan(respiration.dataset1005D1)) =0;
respiration.dataset2005D1(isnan(respiration.dataset2005D1)) =0;
respiration.dataset3005D1(isnan(respiration.dataset3005D1)) =0;
respiration.datasetJoint005D1(isnan(respiration.datasetJoint005D1)) =0;

respiration.dataset2010D1=dataset2010D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset3010D1=dataset3010D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset2010D1(isnan(respiration.dataset2010D1)) =0;
respiration.dataset3010D1(isnan(respiration.dataset3010D1)) =0;

respiration.dataset2050D1=dataset2050D1{modelProteomics.indexRespiration~=0,5}/0.38;
respiration.dataset2050D1(isnan(respiration.dataset2050D1)) =0;
err_metric=[];
errDist=[];
rCoeff=[];
for i=1:8
    [err_metric(i),errDist(:,i),rCoeff(i)] = computeErrorMetric(respiration.proteomics,respiration{:,i+2});
end

boxplot(errDist)
%% Fermentation
fermentation=table();
fermentation.enzNames = modelProteomics.enzNames(modelProteomics.indexFermentation~=0,:);
fermentation.proteomics = fermentationDataset1.glucose(modelProteomics.indexFermentation(modelProteomics.indexFermentation~=0,:));
%fermentation.simulated = enzUsages_reducedYeast.D_04(modelProteomics.indexFermentation~=0,:);

% Remake to percentages
fermentation.proteomics=fermentation.proteomics/sum(fermentationDataset1.glucose(~isnan(fermentationDataset1.glucose)));

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
fermentation.unregulated=(unregulated{modelProteomics.indexFermentation~=0,5}/(0.38*sum(modelProteomics.RellativeAbundance)*sum(unregulated{:,5})));
fermentation.dataset1005D4=(dataset1005D4{modelProteomics.indexFermentation~=0,5}/(0.38*sum(modelProteomics.RellativeAbundance)*sum(dataset1005D4{:,5})));
fermentation.dataset2005D4=(dataset2005D4{modelProteomics.indexFermentation~=0,5}/(0.38*sum(modelProteomics.RellativeAbundance)*sum(dataset2005D4{:,5})));
fermentation.dataset3005D4=(dataset3005D4{modelProteomics.indexFermentation~=0,5}/(0.38*sum(modelProteomics.RellativeAbundance)*sum(dataset3005D4{:,5})));
fermentation.datasetJoint005D4=(datasetJoint005D4{modelProteomics.indexFermentation~=0,5}/(0.38*sum(modelProteomics.RellativeAbundance)*sum(datasetJoint005D4{:,5})));

err_metric=[];
errDist=[];
rCoeff=[];
for i=1:5
    [err_metric(i),errDist(:,i),rCoeff(i)] = computeErrorMetric(fermentation.proteomics,fermentation{:,i+2});
end
boxplot(errDist)
% At 
