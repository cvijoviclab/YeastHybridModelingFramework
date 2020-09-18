function [mutantModel] = constrainFBA2(GeneExpr, PTM, model_in, settings)

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


GeneID = GeneExpr{(GeneExpr{:,2} ~= 0), 3};
sel_ExpressionChanges = GeneExpr{(GeneExpr{:,2} ~= 0), 2};
geneIdx = find(contains(mutantModel.enzGenes, GeneID));
Enzymes = mutantModel.enzymes(geneIdx);
enzUsageList = enzymeUsage_FVA(mutantModel, Enzymes);

for i=1:length(GeneID)
    gene2modIndex   = geneIdx(i); 
    enzyme  = mutantModel.enzymes(gene2modIndex);
    enzRxn = find(contains(mutantModel.rxnNames,enzyme));
    % if enzyme is down regulated upper bounds (enzUsageList(i,4) are adjusted down
    if sel_ExpressionChanges(i) < 0
        mutantModel.ub(enzRxn) = table2array(enzUsageList(i,4)) - ...
            ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse);
        
    % if enzyme is up regulated lower bounds (enzUsageList(i,3) are adjusted up
    elseif sel_ExpressionChanges(i) > 0
        mutantModel.lb(enzRxn) = table2array(enzUsageList(i,3)) + ...
            ((table2array(enzUsageList(i,4))-table2array(enzUsageList(i,3)))* settings.enzymeUse);
    end
    
end


%% Test
% sol = solveLP(mutantModel,1);
% if ~isempty(sol.x)
%     disp(sol.f)
%     disp('Constraints successfully applied')
% end
