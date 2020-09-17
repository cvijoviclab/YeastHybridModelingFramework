function enzUsages_analysis(resultsPath)
model = load('../../models/reduced_ecYeast_fermentation.mat');
model = model.ecModel_ferm;
enzTable = readtable([resultsPath '/enzUsages_reducedYeast.txt'],'delimiter','\t');
enzTable_reg = readtable([resultsPath '/enzUsages_reg_reducedYeast.txt'],'delimiter','\t');
%Get genes and shortNames
enzymes = enzTable.enzymes;
enzGenes   = [];
shortNames = [];
for i=1:length(enzymes)
    enzyme = enzymes{i};
    idx = find(strcmpi(model.enzymes,enzyme));
    enzGenes   = [enzGenes;model.enzGenes(idx)];
    shortNames = [shortNames;model.enzNames(idx)];
end    

variables = {'enzyme' 'enzNames' 'genes' 'subSystems' 'usage' 'usage_regulated' 'usage_ratio'};
ratios  = (enzTable_reg.D_1+1E-15)./(enzTable.D_1+1E-15);
results = table(enzTable.enzymes,shortNames,enzGenes,enzTable.subSystems,enzTable.D_1,enzTable_reg.D_1,ratios,'variableNames',variables);
writetable(results,[resultsPath '/enzUsages_respiration.txt'],'delimiter','\t','QuoteStrings',false);

variables = {'enzyme' 'enzNames' 'genes' 'subSystems' 'usage' 'usage_regulated' 'usage_ratio'};
ratios  = (enzTable_reg.D_4+1E-15)./(enzTable.D_4+1E-15);
results = table(enzTable.enzymes,shortNames,enzGenes,enzTable.subSystems,enzTable.D_4,enzTable_reg.D_4,ratios,'variableNames',variables);
writetable(results,[resultsPath '/enzUsages_fermentation.txt'],'delimiter','\t','QuoteStrings',false);
end