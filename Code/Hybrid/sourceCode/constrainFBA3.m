function [mutantModel] = constrainFBA3(GeneExpr, PTM, model_in, settings)

%written by: Julia M�nch
%date: 2019-12-04
%description: this function modifies the kcats of enzymes of interest via the getMutant function and returns 
%input: 1.-2. Matrix containing simulated gene expression and enzyme activity changes
%       3. WT modle
%       4. percentage of how much turn over number should be changed
%returns: 1. Mutant version of model
%contains functions:
% 1. getMutant

%% make Mutants on enzUsage!
mutantModel=model_in;
load('data/regulatedGenes.mat')

GeneID = GeneExpr{(GeneExpr{:,2} ~= 0), 3};
sel_ExpressionChanges = GeneExpr{(GeneExpr{:,2} ~= 0), 2};
geneIdx = find(contains(mutantModel.enzGenes, GeneID));
Enzymes = mutantModel.enzymes(geneIdx);
enzUsageList = enzymeUsage_FVA(mutantModel, Enzymes);

for i=1:length(GeneID)
    gene2modIndex   = geneIdx(i); 
    level = regulatedGenes.regulation(regulatedGenes.indx==gene2modIndex,:);
%    meanpU=regulatedGenes.pU(regulatedGenes.indx==gene2modIndex,:);
    enzyme  = mutantModel.enzymes(gene2modIndex);
    enzRxn = find(contains(mutantModel.rxnNames,enzyme));
    % if enzyme is down regulated upper bounds (enzUsageList(i,4) are adjusted down
    if sel_ExpressionChanges(i) < 0
        if level == "low"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))/settings.enzymeUse(1);
            mutantModel.ub(enzRxn) = table2array(enzUsageList(i,4)) - ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(1));
        elseif level == "medium"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))/settings.enzymeUse(2);
            mutantModel.ub(enzRxn) = table2array(enzUsageList(i,4)) - ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(2));
        elseif level == "high"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))/settings.enzymeUse(3);
            mutantModel.ub(enzRxn) = table2array(enzUsageList(i,4)) - ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(3));
        end
    % if enzyme is up regulated lower bounds (enzUsageList(i,3) are adjusted up
    elseif sel_ExpressionChanges(i) > 0
        if level == "low"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))*settings.enzymeUse(1);
            mutantModel.lb(enzRxn) = table2array(enzUsageList(i,3)) + ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(1));
        elseif level == "medium"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))*settings.enzymeUse(2);
            mutantModel.lb(enzRxn) = table2array(enzUsageList(i,3)) + ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(2));
        elseif level == "high"
%             mutantModel.ub(enzRxn) = table2array(enzUsageList(i,5))*settings.enzymeUse(3);
            mutantModel.lb(enzRxn) = table2array(enzUsageList(i,3)) + ...
                ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse(3));
        end
    end
    
end

%% Make mutants on kcats

GeneID = PTM{(PTM{:,2} ~= 0), 3};
geneIdx = find(contains(mutantModel.enzGenes, GeneID));
kcat = [settings.kcat{find(contains(settings.kcat{:,1}, GeneID)),2}];

for i=1:length(GeneID)
    gene2modIndex   = geneIdx(i);  
    %Setting enzyme to modify
    if ~isempty(gene2modIndex)
        enzyme  =mutantModel.enzymes(gene2modIndex);
    else
        enzyme = [];
    end
    
    if ~isempty(enzyme)
    enzName    = ['prot_' enzyme{1}];
    enzMetIndx = find(strcmpi(mutantModel.metNames,enzName));
    enzKcats   = find(mutantModel.S(enzMetIndx,:));
    enzKcats   = enzKcats(1:end-1);
    mutantModel.S(enzMetIndx,enzKcats) = cell2mat(kcat(i));
    end
end

%% Test
% sol = solveLP(mutantModel,1);
% if ~isempty(sol.x)
%     disp(sol.f)
%     disp('Constraints successfully applied')
% end
