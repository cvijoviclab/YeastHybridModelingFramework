% Figure1
% Desctiption: This scripts analyzes the enzyme usage from the multiscale 
% simulations and the proteomics data to see understand how the model works
% and how it is regulated. 

%add path to functions
addpath ('functions');
%% Data handeling
% Load unregulated and regulated simulation data and relative abundance.
[rellativeAbundance, D,Dregulated, enzUsages_reducedYeast ,enzUsages_reducedYeastRegulated, geneExprG0, geneExprG1]=loadData();
% Import proteomics data
[fermentationDataset1,respirationDataset1, respirationDataset2]=importProteomics;

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
end
%% Respiration

% UNREGULATED MODEL
respiration=table();
respiration.enzNames = modelProteomics.enzNames(modelProteomics.indexRespiration~=0,:);
respiration.proteomics = respirationDataset1.REF(modelProteomics.indexRespiration(modelProteomics.indexRespiration~=0,:));
% Convert molecule/pgDW to mmol/gDW
avogadrosConstant=6.02214076*10^23;
respiration.proteomics=(respiration.proteomics/avogadrosConstant)*10^(12)*10^3;
respiration.proteomics(isnan(respiration.proteomics)) =0;
respiration.simulated = table2array(D{21}(modelProteomics.indexRespiration~=0,5));
% Correct for satturation
respiration.simulated = respiration.simulated/0.48;
respiration.rank = modelProteomics.G0regulation(modelProteomics.indexRespiration~=0,:);


% REGULATED MODEL
respiration.simulatedR = table2array(Dregulated{21}(modelProteomics.indexRespiration~=0,5));
% Correct for saturation
respiration.simulatedR = respiration.simulatedR/0.48;

% STATISTICS
%Pearson correlation coefficients and significance
[err_metricRespiration,errDistRespiration,rCoeffRespiration] = computeErrorMetric(respiration.proteomics,respiration.simulated);
[err_metricRespirationR,errDistRespirationR,rCoeffRespirationR] = computeErrorMetric(respiration.proteomics,respiration.simulatedR);
disp(['The Person correlation coefficients for the unregulated model is ' num2str(rCoeffRespiration)]);
disp(['The Person correlation coefficients for the regulated model is ' num2str(rCoeffRespirationR)]);
UnregulatedCoeff=1:2000;
RegulatedCoeff=1:2000;
for i=1:2000
    [~,~,UnregulatedCoeff(i)] = computeErrorMetric(respiration.proteomics,randsample(respiration.simulated,size(respiration.simulated,1)));
    [~,~,RegulatedCoeff(i)] = computeErrorMetric(respiration.proteomics,randsample(respiration.simulatedR,size(respiration.simulatedR,1)));
end
disp(['P-value for PCC of the unregulated simulation: ' num2str(sum(UnregulatedCoeff>rCoeffRespiration)/2000)])
disp(['P-value for PCC of the regulated simulation: ' num2str(sum(RegulatedCoeff>rCoeffRespirationR)/2000)])

figure()
h=histogram(UnregulatedCoeff,'FaceColor','#0072BD');
hold on
xline(rCoeffRespiration,'Color','#A2142F','LineWidth',3);
text(rCoeffRespiration+0.01,max(h.Values)*1.05,'Observed PCC')
hold off
title('Permutation test for PCC - Unregulated','FontSize',16);
ylabel('Frequency','FontSize',12);
xlabel('PCC','FontSize',12);

figure()
h=histogram(RegulatedCoeff,'FaceColor','#0072BD');
hold on
xline(rCoeffRespirationR,'Color','#A2142F','LineWidth',3);
text(rCoeffRespirationR+0.01,max(h.Values)*0.85,'Observed PCC')
hold off
title('Permutation test for PCC - Regulated','FontSize',16);
ylabel('Frequency','FontSize',12);
xlabel('PCC','FontSize',12);

[h,p,ks2stat] = kstest2(respiration.proteomics,respiration.simulated);
disp(['The null hypothesis, equal distributions of the unregulated model and the experimental data, H: ' num2str(h)]);
disp([ 'P-value: ' num2str(p) ]);
disp(['KS statistic: ' num2str(ks2stat)]);

[h,p,ks2stat] = kstest2(respiration.proteomics,respiration.simulatedR);
disp(['The null hypothesis, equal distributions of the regulated model and the experimental data, H: ' num2str(h)]);
disp([ 'P-value: ' num2str(p) ]);
disp(['KS statistic: ' num2str(ks2stat)]);

% Visulizing summary statistics
figure()
boxplot([errDistRespiration errDistRespirationR],'label',{'Unregulated','Regulated'});
title('Respiration','FontSize',16);
ylabel('Log_{10}(Predicted/Experimental)','FontSize',12);
xlabel('Error distribution','FontSize',12);
text([1,2],[-1.5,-1],{['r=' num2str(round(err_metricRespiration,2))], ['r=' num2str(round(err_metricRespirationR,2))]}, 'HorizontalAlignment','center');

% Visulizing individual error statistics
figure()
scatter(1:length(errDistRespirationR),errDistRespirationR);
title('Respiration','FontSize',16);
ylabel('Log_{10}(Predicted/Experimental)','FontSize',12);
xticks(1:length(errDistRespirationR));
xticklabels(respiration.enzNames);
xtickangle(90);
disp('Following proteins are over- or underpredicted by the model with over 5 orders of magnitude:');
disp(respiration(errDistRespirationR>5|errDistRespirationR<-5,[1,2,4,5]));
disp('Following proteins are over- or underpredicted by the model with over 1 order of magnitude:');
disp(respiration(errDistRespirationR>1|errDistRespirationR<-1,[1,2,4,5]));

% GENERATING FIGURES FOR COMPARING PROTEINS
figure()
for i = 1:size(respiration,1)
b=bar(categorical(respiration{i,1}),respiration{i,[2,3,5]});
b(1).FaceColor = [153, 153, 153]/255;
b(2).FaceColor  = [230, 159, 0]/255;
b(3).FaceColor = [204, 121, 167]/255;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold');  
saveas(gcf,['figures/respiration/' respiration{i,1}{1,1} '.png']);
end

%% Fermentation
 
% UNREGULATED MODEL
fermentation=table();
fermentation.enzNames = modelProteomics.enzNames(modelProteomics.indexFermentation~=0,:);
fermentation.proteomics = fermentationDataset1.glucose(modelProteomics.indexFermentation(modelProteomics.indexFermentation~=0,:));
fermentation.simulated = enzUsages_reducedYeast.D_4(modelProteomics.indexFermentation~=0,:);
% Correct for satturation
fermentation.simulated = fermentation.simulated/0.48;
fermentation.rank = modelProteomics.G1regulation(modelProteomics.indexFermentation~=0,:);
% Remake to percentages
fermentation.proteomics=fermentation.proteomics/sum(fermentationDataset1.glucose(~isnan(fermentationDataset1.glucose)));
% Map index in model for rellative abundance
modelProteomics.indexRellativeAbundance = zeros(size(modelProteomics,1),1);
for i=1:size(modelProteomics,1)
    j = (find(ismember(rellativeAbundance.string_external_id,['4932.' char(geneExprG1.enzGenes(i))])));
    if ~isempty(j)
        modelProteomics.indexRellativeAbundance(i) = j;
    end
end
% Fetch relative abundance
modelProteomics.RellativeAbundance = zeros(size(modelProteomics,1),1);
for i=1:size(modelProteomics,1)
    if modelProteomics.indexRellativeAbundance(i)~=0
        modelProteomics.RellativeAbundance(i) = table2array(rellativeAbundance(modelProteomics.indexRellativeAbundance(i),3));
    end
end
%rellative abundance is in ppm /10000-->%
modelProteomics.RellativeAbundance=modelProteomics.RellativeAbundance/10000;
fermentation.simulated=fermentation.simulated/sum(modelProteomics.RellativeAbundance);

% REGULATED MODEL
fermentation.simulatedR = enzUsages_reducedYeastRegulated.D_4(modelProteomics.indexFermentation~=0,:);
fermentation.simulatedR = fermentation.simulatedR/0.48;
fermentation.simulatedR=fermentation.simulatedR/sum(modelProteomics.RellativeAbundance);

% STATISTICS
%Pearson correlation coefficients and significance
[err_metricFermentation,errDistFermentation,rCoeffFermentation] = computeErrorMetric(fermentation.proteomics,fermentation.simulated);
[err_metricFermentationR,errDistFermentationR,rCoeffFermentationR] = computeErrorMetric(fermentation.proteomics,fermentation.simulatedR);
disp(['The Person correlation coefficients for the unregulated model is ' num2str(rCoeffFermentation)]);
disp(['The Person correlation coefficients for the regulated model is ' num2str(rCoeffFermentationR)]);
UnregulatedCoeff=1:2000;
RegulatedCoeff=1:2000;
for i=1:2000
    [~,~,UnregulatedCoeff(i)] = computeErrorMetric(fermentation.proteomics,randsample(fermentation.simulated,size(fermentation.simulated,1)));
    [~,~,RegulatedCoeff(i)] = computeErrorMetric(fermentation.proteomics,randsample(fermentation.simulatedR,size(fermentation.simulatedR,1)));
end
disp(['P-value for PCC of the unregulated simulation: ' num2str(sum(UnregulatedCoeff>rCoeffFermentation)/2000)])
disp(['P-value for PCC of the regulated simulation: ' num2str(sum(RegulatedCoeff>rCoeffFermentationR)/2000)])

figure()
h=histogram(UnregulatedCoeff,'FaceColor','#0072BD');
hold on
xline(rCoeffFermentation,'Color','#A2142F','LineWidth',3);
text(rCoeffFermentation+0.01,max(h.Values)*1.05,'Observed PCC')
hold off
title('Permutation test for PCC - Unregulated','FontSize',16);
ylabel('Frequency','FontSize',12);
xlabel('PCC','FontSize',12);

figure()
h=histogram(RegulatedCoeff,'FaceColor','#0072BD');
hold on
xline(rCoeffFermentationR,'Color','#A2142F','LineWidth',3);
text(rCoeffFermentationR+0.01,max(h.Values)*0.85,'Observed PCC')
hold off
title('Permutation test for PCC - Regulated','FontSize',16);
ylabel('Frequency','FontSize',12);
xlabel('PCC','FontSize',12);


[h,p,ks2stat] = kstest2(fermentation.proteomics,fermentation.simulated);
disp(['The null hypothesis, equal distributions of the unregulated model and the experimental data, H: ' num2str(h)]);
disp([ 'P-value: ' num2str(p) ]);
disp(['KS statistic: ' num2str(ks2stat)]);

[h,p,ks2stat] = kstest2(fermentation.proteomics,fermentation.simulatedR);
disp(['The null hypothesis, equal distributions of the regulated model and the experimental data, H: ' num2str(h)]);
disp([ 'P-value: ' num2str(p) ]);
disp(['KS statistic: ' num2str(ks2stat)]);

% Visulizing summary statistics
figure()
boxplot([errDistFermentation, errDistFermentationR],'label',{'Unregulated','Regulated'});
title('Fermentation','FontSize',16);
ylabel('Log_{10}(Predicted/Experimental)','FontSize',12);
xlabel('Error distribution','FontSize',12);
text([1,2],[-1,-1],{['r=' num2str(round(err_metricFermentation,2))], ['r=' num2str(round(err_metricFermentationR,2))]}, 'HorizontalAlignment','center');

% Visulizing individual error statistics
figure()
scatter(1:length(errDistFermentationR),errDistFermentationR);
title('Fermentation','FontSize',16);
ylabel('Log_{10}(Predicted/Experimental)','FontSize',12);
xticks(1:length(errDistFermentationR));
xticklabels(fermentation.enzNames);
xtickangle(90);
disp('Following proteins are over- or underpredicted by the model with over 5 orders of magnitude:');
disp(fermentation(errDistFermentationR>5|errDistFermentationR<-5,[1,2,4,5]));
disp('Following proteins are over- or underpredicted by the model with over 1 order of magnitude:');
disp(fermentation(errDistFermentationR>1|errDistFermentationR<-1,[1,2,4,5]));

% GENERATING FIGURES FOR COMPARING PROTEINS
figure()
for i = 1:size(fermentation,1)
b=bar(categorical(fermentation{i,1}),fermentation{i,[2,3,5]});
b(1).FaceColor = [153, 153, 153]/255;
b(2).FaceColor  = [230, 159, 0]/255;
b(3).FaceColor = [204, 121, 167]/255;
set(gca,'FontSize',30);
set(gca,'FontWeight','bold');  
saveas(gcf,['figures/fermentation/' fermentation{i,1}{1,1} '.png']);
end
%% flux analysis
FluxUnregulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure1/Simulations/multiscaleUnregulated20200612/fluxDist_reducedYeast.txt');
FluxRegulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure1/Simulations/multiscaleRegulated20200612nr2/fluxDist_reducedYeast.txt');

FluxRellativeDiff=FluxUnregulated(:,1:5);
FluxRellativeDiff.Flux=(FluxRegulated{:,6:86}-FluxUnregulated{:,6:86});
