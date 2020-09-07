function [modelProteomics]=addProtData(geneExprG0,proteomics)

% Work with percentages
modelProteomics = geneExprG0(:,[6 1 7:8]);
proteomicsData=table();
empty={0,0,0,0,0,0,0,0,0,0};
for i = 1:size(modelProteomics,1)
    if modelProteomics.proteomicsIndex(i) ~= 0
        proteomicsData(i,:)=proteomics(i,4:13);
    else
        proteomicsData(i,:)=cell2table(empty);
    end
end
modelProteomics(:,5:14)=proteomicsData(:,1:10);
modelProteomics.Properties.VariableNames(5:14) = proteomicsData.Properties.VariableNames(1:10);
% Convert to percentage of usage
sumProtein=sum(table2array(modelProteomics(:,5:14)));
for i=1:length(sumProtein)
    modelProteomics{:,4+i}=modelProteomics{:,4+i}/sumProtein(i);
end

modelProteomics.TotalAve=mean(modelProteomics{:,[5:8 12:14]},2);%based on three last and three first timepoints
modelProteomics.respirationAve=mean(modelProteomics{:,12:14},2);%based on three last timepoints
modelProteomics.fermentationAve=mean(modelProteomics{:,5:8},2);%based on three first timepoints