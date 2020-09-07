% Desctiption: This scripts analyzes the enzyme usage from the multiscale 
% simulations and the proteomics data to see understand how the model works
% and how it should be regulated. 
close all
clear all
%add path to functions
addpath ('../sourceCode');
addpath ('functions');
%% Data handeling
% Load data
[rellativeAbundance,proteomics, D,Dregulated, enzUsages_reducedYeast,enzUsages_reducedYeastRegulated, geneExprG0, geneExprG1]=loadData()
% Preprocess data 
[proteomics, D, enzUsages_reducedYeast, geneExprG0, geneExprG1]=preProcess(proteomics, D, enzUsages_reducedYeast, geneExprG0, geneExprG1);

%% Visulize and analyze simulated data without regluation
dilutionrate=0:(0.4/80):0.4; %i=59, D=0.29
[geneExprG0]= addSimData(geneExprG0,D);
[geneExprG1]= addSimData(geneExprG1,D);

% Visualiza ranges
G0regulatedGenes = geneExprG0(table2array(geneExprG0(:,2)) ~= 0,[6 7 1:3 9:12]);
G1regulatedGenes = geneExprG1(table2array(geneExprG1(:,2)) ~= 0,[6 7 1:3 9:12]);
% plot the parameter space statistics first spaning all dilution rates and
% then only for G0 or G1. Right plot spcae is fermentation left is
% respiration.
plotfigure1(G0regulatedGenes,G1regulatedGenes)
% plot pU statistics for all enzymes and regulated enzymes spaning all
% dilutionrates and the specific G0 and G1.
plotfigure2(geneExprG0,G0regulatedGenes,G1regulatedGenes)

% When vizualizing pU for the regulated enzymes in
disp('Some means of pU are 0.')
disp(['Of those ' num2str(sum(mean(G0regulatedGenes.pU')'==0 & G0regulatedGenes.rank>0)) ' are upregulated and ' num2str(sum(mean(G0regulatedGenes.pU')'==0 & G0regulatedGenes.rank<0)) ' are downregulated'])
disp(G0regulatedGenes(mean(G0regulatedGenes.pU')'==0  & G0regulatedGenes.rank>0,[3]));
disp(G0regulatedGenes(mean(G0regulatedGenes.pU')'==0  & G0regulatedGenes.rank<0,[3]));

% look in what range of the parameter space the pU is
% picking 2 dilution rates D_0.1 = respiration(20) D_0.35 = fermentation(71). (randomly picked) figure(3)
plotfigure3(G0regulatedGenes,G1regulatedGenes)
%all dilution rates
plotfigure4(G0regulatedGenes,G1regulatedGenes)

%% New proteomics set
[fermentationDataset1,respirationDataset1]=importNewProteomics;
modelProteomics = geneExprG0(:,[6 1 7]);
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
end
%Respiration
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

% Remake to percentages
% respiration.proteomics=respiration.proteomics/sum(respiration.proteomics(~isnan(respiration.proteomics),:));
% respiration.simulated=respiration.simulated/sum(respiration.simulated);

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
% disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(respiration.proteomics-respiration.simulated))*100) '%'])
% disp(['Mean difference between data simulated and proteomics data is ' ...
% num2str((mean(respiration.proteomics(respiration.rank>0)-respiration.simulated(respiration.rank>0))+ mean(respiration.simulated(respiration.rank<0)-respiration.proteomics(respiration.rank<0)))/2/(sum(table2array(D{21}(:,5))))*100) '%'])

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
disp(['Median difference between bounds simulated and proteomics data is '...
    num2str(median([upregulated; downregulated])*100) '%'])

list4=table();
list4.Names= [respiration.enzNames(respiration.rank>0);respiration.enzNames(respiration.rank<0)];
list4.distance=[upregulated; downregulated];
list4.distance(list4.distance<0)=0;
list4.class=kmeans(list4.distance,3);
disp( mean(list4.distance)*100)
disp( num2str(median(list4.distance(list4.class==1))*100)) 
disp(num2str(median(list4.distance(list4.class==2))*100)) 
disp(num2str(median(list4.distance(list4.class==3))*100));

%disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(respiration.proteomics(respiration.rank~=0)-respiration.simulated(respiration.rank~=0)))*100) '%'])

%disp(['Mean foldchange between data simulated and proteomics data is ' num2str(mean(abs(respiration.simulated./respiration.proteomics))*100) '%'])
%disp(['Max difference between data simulated and proteomics data is ' num2str(max(abs(respiration.proteomics-respiration.simulated))*100) '%'])
list1=(respiration.proteomics);
list1(isinf(list1))=0;
list2=(respiration.simulated);
list2(isinf(list2))=0;
list3=(list2-list1)./list1;
list=table();
list.name=respiration.enzNames(list3>10 | list3 <-0.1);
list.logdifference=list3(list3>10| list3 < -0.1);

%Bhattacharyya distance
q = respiration.proteomics;
p = respiration.simulated;
BC=sum(sqrt(p.*q));
D_b=-log(BC);
disp(['The Bhattacharyya distance during respiration is ' num2str(D_b)])


%Fermentation
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
disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(fermentation.proteomics-fermentation.simulated))*100) '%'])
disp(['Mean difference between data simulated and proteomics data is ' ...
num2str((mean(fermentation.proteomics(fermentation.rank>0)-fermentation.simulated(fermentation.rank>0))+ mean(fermentation.simulated(fermentation.rank<0)-fermentation.proteomics(fermentation.rank<0)))/2*100) '%'])
%disp(['Mean foldchange between data simulated and proteomics data is ' num2str(median(abs(fermentation.simulated./fermentation.proteomics))*100) '%'])
%disp(['Max difference between data simulated and proteomics data is ' num2str(max(abs(fermentation.proteomics-fermentation.simulated))*100) '%'])

%Bhattacharyya distance
q = fermentation.proteomics;
p = fermentation.simulated;
BC=sum(sqrt(p.*q));
D_b=-log(BC);
disp(['The Bhattacharyya distance during respiration is ' num2str(D_b)])
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

% Remake to percentages
% respiration.proteomics=respiration.proteomics/sum(respiration.proteomics(~isnan(respiration.proteomics),:));
% respiration.simulated=respiration.simulated/sum(respiration.simulated);
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
disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(respiration.proteomics-respiration.simulated))*100) '%'])
disp(['Mean difference between data simulated and proteomics data is ' ...
num2str((mean(respiration.proteomics(respiration.rank>0)-respiration.simulated(respiration.rank>0))+ mean(respiration.simulated(respiration.rank<0)-respiration.proteomics(respiration.rank<0)))/2*100) '%'])

list1=respiration.proteomics;
list1(isinf(list1))=0;
list2=respiration.simulated;
list2(isinf(list2))=0;
list3=(list2-list1)./list1;
list2=table();
list2.name=respiration.enzNames(list3>10 | list3 <-0.1);
list2.logdifference=list3(list3>10 | list3 <-0.1);
%disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(respiration.proteomics(respiration.rank~=0)-respiration.simulated(respiration.rank~=0)))*100) '%'])


%disp(['Mean foldchange between data simulated and proteomics data is ' num2str(mean(abs(respiration.simulated./respiration.proteomics))*100) '%'])
%disp(['Max difference between data simulated and proteomics data is ' num2str(max(abs(respiration.proteomics-respiration.simulated))*100) '%'])

%Bhattacharyya distance
q = respiration.proteomics;
p = respiration.simulated;
BC=sum(sqrt(p.*q));
D_b=-log(BC);
disp(['The Bhattacharyya distance during respiration is ' num2str(D_b)])

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
disp(['Mean difference between data simulated and proteomics data is ' num2str(mean(abs(fermentation.proteomics-fermentation.simulated))*100) '%'])

%Bhattacharyya distance
q = fermentation.proteomics;
p = fermentation.simulated;
BC=sum(sqrt(p.*q));
D_b=-log(BC);
disp(['The Bhattacharyya distance during respiration is ' num2str(D_b)])


[err_metricRespirationR,errDistRespirationR,rCoeffRespirationR] = computeErrorMetric(respiration.proteomics,respiration.simulated);
[err_metricFermentationR,errDistFermentationR,rCoeffFermentationR] = computeErrorMetric(fermentation.proteomics,fermentation.simulated);


