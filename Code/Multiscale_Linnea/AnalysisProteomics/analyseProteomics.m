% Desctiption: This scripts analyzes the enzyme usage from the multiscale 
% simulations and the proteomics data to see understand how the model works
% and how it is regulated. 
close all
clear all
%add path to functions
addpath ('../sourceCode');
addpath ('functions');
%% Data handeling
% Load unregulated and regulated simulation data and relative abundance.
[rellativeAbundance, D,Dregulated, enzUsages_reducedYeast ,enzUsages_reducedYeastRegulated, geneExprG0, geneExprG1]=loadData();
% Import proteomics data
[fermentationDataset1,respirationDataset1, respirationDataset2]=importNewProteomics;

% Import proteomics data for 
%% New proteomics set
modelProteomics = geneExprG0(:,[1 3]);
modelProteomics.G0regulation=geneExprG0.rank;
modelProteomics.G1regulation=geneExprG1.rank;

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
    j = (find(ismember(respirationDataset2.Protein, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        modelProteomics.indexRespiration2(i) = j;
    end
end
%% Respiration
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration~=0,:);
respiration.proteomics = respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexRespiration~=0,:));
% Convert molecule/pgDW to mmol/gDW
avogadrosConstant=6.02214076*10^23;
respiration.proteomics=(respiration.proteomics/avogadrosConstant)*10^(12)*10^3;
respiration.simulated = table2array(D{21}(modelProteomics.indexRespiration~=0,5));
% Correct for satturation
respiration.simulated = respiration.simulated/0.38;
respiration.simulatedRange = table2array(D{21}(modelProteomics.indexRespiration~=0,2));
respiration.simulatedMinU = table2array(D{21}(modelProteomics.indexRespiration~=0,3));
respiration.simulatedMaxU = table2array(D{21}(modelProteomics.indexRespiration~=0,4));
respiration.rank = modelProteomics.G0regulation(modelProteomics.indexRespiration~=0,:);

figure()
colorsG0 = zeros(size(respiration,1),3);
colorsG0(respiration.rank>0 ,2)= 1;
colorsG0(respiration.rank<0,1)= 1;
scatter(log10(respiration.proteomics),log10(respiration.simulated),[],colorsG0,'filled')
hold on
plot(linspace(0,-9,50),linspace(0,-9,50),'LineWidth',3)
plot(linspace(-9,-1,50),linspace(-8,-0,50),'LineWidth',3)
plot(linspace(0,-8,50),linspace(0,-8,50)-1,'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)
respiration.proteomics(isnan(respiration.proteomics)) =0;

upregulated=zeros(sum(respiration.rank>0),1);
vector = 1:length(respiration.proteomics)
vector = vector(respiration.rank>0)
for i = 1:sum(respiration.rank>0)
    if respiration.proteomics(vector(i))<respiration.simulatedMaxU(vector(i)) %&& respiration.proteomics(vector(i))>respiration.simulatedMinU(vector(i))
        upregulated(i)=(respiration.proteomics(vector(i))-respiration.simulatedMinU(vector(i)))/respiration.simulatedRange(vector(i));
    else
        upregulated(i)=1
    end
end
downregulated=[zeros(sum(respiration.rank<0),1)];
vector = 1:length(respiration.proteomics)
vector = vector(respiration.rank<0)
for i=1:sum(respiration.rank<0)
    if respiration.proteomics(vector(i))>respiration.simulatedMinU(vector(i)) %&& respiration.proteomics(vector(i))<respiration.simulatedMaxU(vector(i))
        downregulated(i)=(respiration.simulatedMaxU(vector(i))-respiration.proteomics(vector(i)))/respiration.simulatedRange(vector(i));
    else
        downregulated(i)=1
    end
end
respirationUR=respiration;

%% Fermentation
fermentation=table();
fermentation.enzNames = modelProteomics.enzNames(modelProteomics.indexFermentation~=0,:);
fermentation.proteomics = fermentationDataset1.glucose(modelProteomics.indexFermentation(modelProteomics.indexFermentation~=0,:));
fermentation.simulated = enzUsages_reducedYeast.D_04(modelProteomics.indexFermentation~=0,:);
% Correct for satturation
fermentation.simulated = fermentation.simulated/0.38;
fermentation.rank = modelProteomics.G1regulation(modelProteomics.indexFermentation~=0,:);
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
fermentation.simulated=fermentation.simulated/sum(modelProteomics.RellativeAbundance);

figure()
colorsG1 = zeros(size(fermentation,1),3);
colorsG1(fermentation.rank>0 ,2)= 1;
colorsG1(fermentation.rank<0,1)= 1;
scatter(log10(fermentation.proteomics),log10(fermentation.simulated),[],colorsG1,'filled')
hold on
plot(linspace(0,-6,50),linspace(0,-6,50),'LineWidth',3)
plot(linspace(-1,-6,50),linspace(0,-5,50),'LineWidth',3)
plot(linspace(0,-5,50),linspace(-1,-6,50),'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)
fermentation.proteomics(isnan(fermentation.proteomics)) =0;

[err_metricRespiration,errDistRespiration,rCoeffRespiration] = computeErrorMetric(respiration.proteomics,respiration.simulated);
[err_metricFermentation,errDistFermentation,rCoeffFermentation] = computeErrorMetric(fermentation.proteomics,fermentation.simulated);


%% redo proteomics with simulated data with regulation

%Respiration
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration~=0,:);
respiration.proteomics = respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexRespiration~=0,:));
% Convert molecule/pgDW to mmol/gDW
avogadrosConstant=6.02214076*10^23;
respiration.proteomics=(respiration.proteomics/avogadrosConstant)*10^(12)*10^3;
respiration.simulated = table2array(Dregulated{21}(modelProteomics.indexRespiration~=0,5));
% Correct for satturation
respiration.simulated = respiration.simulated/0.38;
respiration.rank = modelProteomics.G0regulation(modelProteomics.indexRespiration~=0,:);


figure()
colorsG0 = zeros(size(respiration,1),3);
colorsG0(respiration.rank>0 ,2)= 1;
colorsG0(respiration.rank<0,1)= 1;
scatter(log10(respiration.proteomics),log10(respiration.simulated),[],colorsG0,'filled')
hold on
plot(linspace(0,-9,50),linspace(0,-9,50),'LineWidth',3)
plot(linspace(-9,-1,50),linspace(-8,-0,50),'LineWidth',3)
plot(linspace(0,-8,50),linspace(0,-8,50)-1,'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)
respiration.proteomics(isnan(respiration.proteomics)) =0;

%Fermentation
fermentation=table();
fermentation.enzNames = modelProteomics.enzNames(modelProteomics.indexFermentation~=0,:);
fermentation.proteomics = fermentationDataset1.glucose(modelProteomics.indexFermentation(modelProteomics.indexFermentation~=0,:));
fermentation.simulated = enzUsages_reducedYeastRegulated.D_04(modelProteomics.indexFermentation~=0,:);
% Correct for satturation
fermentation.simulated = fermentation.simulated/0.38;
fermentation.rank = modelProteomics.G1regulation(modelProteomics.indexFermentation~=0,:);
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
fermentation.simulated=fermentation.simulated/sum(modelProteomics.RellativeAbundance);

figure()
colorsG1 = zeros(size(fermentation,1),3);
colorsG1(fermentation.rank>0 ,2)= 1;
colorsG1(fermentation.rank<0,1)= 1;
scatter(log10(fermentation.proteomics),log10(fermentation.simulated),[],colorsG1,'filled')
hold on
plot(linspace(0,-8,50),linspace(0,-8,50),'LineWidth',3)
plot(linspace(-8,-1,50),linspace(-7,-0,50),'LineWidth',3)
plot(linspace(0,-7,50),linspace(0,-7,50)-1,'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)
fermentation.proteomics(isnan(fermentation.proteomics)) =0;

[err_metricRespirationR,errDistRespirationR,rCoeffRespirationR] = computeErrorMetric(respiration.proteomics,respiration.simulated);
[err_metricFermentationR,errDistFermentationR,rCoeffFermentationR] = computeErrorMetric(fermentation.proteomics,fermentation.simulated);
figure()
boxplot([errDistRespiration errDistRespirationR]);
figure()
boxplot([errDistFermentation errDistFermentationR]);

figure()
scatter([1:length(errDistRespirationR)],errDistRespirationR);
respiration(errDistRespirationR>5|errDistRespirationR<-5,:) 
vorstPredicted=respiration(errDistRespirationR>5|errDistRespirationR<-5,1) 
badlyPredicted=respiration(errDistRespirationR>1|errDistRespirationR<-1,1) 
respiration(errDistRespirationR>1|errDistRespirationR<-1,:) 
figure()
scatter([1:length(errDistFermentationR)],errDistFermentationR);
fermentation(errDistFermentationR>5|errDistFermentationR<-5,:)
vorstPredicted=[vorstPredicted; fermentation(errDistFermentationR>5|errDistFermentationR<-5,1)];


[h,p,ks2stat] = kstest2(respiration.simulated,respirationUR.simulated);
[h,p,ks2stat] = kstest2(errDistRespiration,errDistRespirationR);
[h,p,ks2stat] = kstest2(errDistFermentationR,errDistFermentation);

% respiration 2
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration2~=0,:);
respiration.proteomics2 = respirationDataset2.FinalQuant(modelProteomics.indexRespiration2(modelProteomics.indexRespiration2~=0,:));
% Remake to percentages
respiration.proteomics2=respiration.proteomics2/sum(respirationDataset2.FinalQuant(~isnan(respirationDataset2.FinalQuant)));
respiration.simulated = table2array(Dregulated{21}(modelProteomics.indexRespiration2~=0,5));
% Correct for satturation
respiration.simulated = respiration.simulated/0.38;
respiration.rank = modelProteomics.G0regulation(modelProteomics.indexRespiration2~=0,:);
respiration.simulated=respiration.simulated/(sum(table2array(Dregulated{21}(:,5)))*sum(modelProteomics.RellativeAbundance));
respiration.proteomics2(isnan(respiration.proteomics2))=0;
[err_metricRespirationR,errDistRespirationR,rCoeffRespirationR] = computeErrorMetric(respiration.proteomics2,respiration.simulated);
figure()
boxplot([errDistRespirationR]);
figure()
scatter([1:length(errDistRespirationR)],errDistRespirationR);
respiration(errDistRespirationR>5|errDistRespirationR<-5,:) 
respiration(errDistRespirationR>1|errDistRespirationR<-1,:) 
% The null hypothesis: that the data in vectors x1 and x2 are from the same continuous distribution, using the two-sample Kolmogorov-Smirnov test.
%     If h = 1, this indicates the rejection of the null hypothesis at the Alpha significance level.
%     If h = 0, this indicates a failure to reject the null hypothesis at the Alpha significance level.
[h,p,ks2stat] = kstest2(errDistRespiration,errDistRespirationR)
%% Deletion strains
[rapamycin,SNF1deletion]=loadDEletionProteomics();
TORDeletionSimulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/AnalysisProteomics/Simulations/multiscaleDTor1Tor220200601/enzUsage_VarAnalysis_D_4.txt');
Snf1DeletionSimulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/AnalysisProteomics/Simulations/multiscaleDSnf120200601/enzUsage_VarAnalysis_D_1.txt');

modelProteomics.indexTORDeletion= zeros(size(modelProteomics,1),1);
modelProteomics.indexSnf1Deletion = zeros(size(modelProteomics,1),1);
for i=1:size(modelProteomics,1)
    j = (find(ismember(rapamycin.ORF, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        modelProteomics.indexTORDeletion(i) = j;
    end
    j = (find(ismember(SNF1deletion.YGL205W, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        modelProteomics.indexSnf1Deletion(i) = j;
    end
end

TORDeletion=table();
TORDeletion.enzNames = modelProteomics.enzNames(modelProteomics.indexTORDeletion~=0,:);
TORDeletion.proteomics = rapamycin.RAP700(modelProteomics.indexTORDeletion(modelProteomics.indexTORDeletion~=0,:));
TORDeletion.proteomics(isnan(TORDeletion.proteomics))=0;
TORDeletion.simulated=TORDeletionSimulated.pU(modelProteomics.indexTORDeletion~=0,:);
TORDeletion.unregulated= table2array(D{21}(modelProteomics.indexTORDeletion~=0,5));
% Convert to percentages 
TORDeletion.proteomics=TORDeletion.proteomics/sum(rapamycin.RAP700(~isnan(rapamycin.RAP700)));
TORDeletion.simulated=TORDeletion.simulated/(sum(TORDeletionSimulated.pU)*sum(modelProteomics.RellativeAbundance));
TORDeletion.unregulated=TORDeletion.unregulated/(sum(table2array(D{21}(:,5)))*sum(modelProteomics.RellativeAbundance));


% Correct for satturation
TORDeletion.simulated = TORDeletion.simulated/0.38;
TORDeletion.unregulated = TORDeletion.unregulated/0.38;

figure()
scatter(log10(TORDeletion.proteomics),log10(TORDeletion.simulated))
hold on
plot(linspace(0,-8,50),linspace(0,-8,50),'LineWidth',3)
plot(linspace(-8,-1,50),linspace(-7,-0,50),'LineWidth',3)
plot(linspace(0,-7,50),linspace(0,-7,50)-1,'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)

figure()
scatter(log10(TORDeletion.proteomics),log10(TORDeletion.unregulated))
hold on
plot(linspace(0,-8,50),linspace(0,-8,50),'LineWidth',3)
plot(linspace(-8,-1,50),linspace(-7,-0,50),'LineWidth',3)
plot(linspace(0,-7,50),linspace(0,-7,50)-1,'LineWidth',3)
hold off
xlabel('Proteomics Data','FontSize',18)
ylabel('Simulated Data','FontSize',18)

[err_metricTORDeletion,errDistTORDeletion,rCoeffTORDeletion] = computeErrorMetric(TORDeletion.proteomics,TORDeletion.simulated);
[err_metricTORDeletionUR,errDistTORDeletionUR,rCoeffTORDeletionUR] = computeErrorMetric(TORDeletion.proteomics,TORDeletion.unregulated);
figure()
boxplot([errDistTORDeletion errDistTORDeletionUR]);
scatter([1:length(errDistTORDeletion)],errDistTORDeletion);
vorstPredicted=[vorstPredicted; TORDeletion(errDistTORDeletion>2,1)];

Snf1Deletion=table();
Snf1Deletion.enzNames = modelProteomics.enzNames(modelProteomics.indexSnf1Deletion~=0 & modelProteomics.indexRespiration~=0,:);
Snf1Deletion.proteomics = SNF1deletion.VarName7(modelProteomics.indexSnf1Deletion(modelProteomics.indexSnf1Deletion~=0 & modelProteomics.indexRespiration~=0,:));
Snf1Deletion.simulation=Snf1DeletionSimulated.pU(modelProteomics.indexSnf1Deletion~=0 & modelProteomics.indexRespiration~=0,:);
Snf1Deletion.wt=respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexSnf1Deletion~=0 & modelProteomics.indexRespiration~=0,:));
%Correct for unit
avogadrosConstant=6.02214076*10^23;
Snf1Deletion.wt=(Snf1Deletion.wt/avogadrosConstant)*10^(12)*10^3;
% Correct for satturation
Snf1Deletion.simulation = Snf1Deletion.simulation/0.38;
Snf1Deletion.simulation(Snf1Deletion.simulation==0)=10^-9;
Snf1Deletion.log2diff=log2(Snf1Deletion.wt)-log2(Snf1Deletion.simulation);

[err_metricSnf1Deletion,errDistSnf1Deletion,rCoeffSnf1Deletion] = computeErrorMetric(Snf1Deletion.proteomics,Snf1Deletion.log2diff);
scatter([1:length(errDistSnf1Deletion)],errDistSnf1Deletion);
Snf1Deletion(errDistSnf1Deletion>2,:)
Snf1Deletion.rank=zeros(size(Snf1Deletion,1),1);
for i=1:size(Snf1Deletion,1)
    j = (find(ismember(geneExprG0.enzNames, Snf1Deletion.enzNames(i))));
    Snf1Deletion.rank(i)=j;
end
Snf1Deletion.rank=geneExprG0.rank(Snf1Deletion.rank);
Snf1Deletion(errDistSnf1Deletion>2,:)
vorstPredicted=[vorstPredicted; Snf1Deletion(errDistSnf1Deletion>2,1)];

[C,ia,ic] = unique(vorstPredicted);
a_counts = array2table(accumarray(ic,1));
value_counts = [C, a_counts]
value_counts(value_counts.Var1>1,:)