function [Mutant] = constrainFBA(GeneExpr, model_in, perc)

%written by: Julia M�nch
%date: 2019-12-04
%description: this function modifies the kcats of enzymes of interest via the getMutant function and returns 
%input: 1.-2. Matrix containing simulated gene expression and enzyme activity changes
%       3. WT modle
%       4. percentage of how much turn over number should be changed
%returns: 1. Mutant version of model
%contains functions:
% 1. getMutant

%% make Mutants!

GeneID = GeneExpr{(GeneExpr{:,2} ~= 0), 3};
sel_ExpressionChanges = GeneExpr{(GeneExpr{:,2} ~= 0), 2};
mod_action = num2cell(ones(length(GeneID),1));
mod_action(:,:) = num2cell(2); %as we want to change Kcat directly
factor = zeros(length(GeneID),1);
for i=1:length(factor)
    if sel_ExpressionChanges(i) < 0
        factor(i) = 1-perc; %0.97
    elseif sel_ExpressionChanges(i) > 0
        factor(i) = 1+perc; %1.03
    end
end

modifications = [GeneID, mod_action, num2cell(factor)];

model=model_in;
geneIdx = find(contains(model.enzGenes, GeneID));
Enzymes = model.enzymes(geneIdx);
rxnIndx = find(contains(model.rxnNames,Enzymes));

%if you want to change boundaries: use enzyme usage (enzUsage) or
%fluxes(enzUs) as third argument in getMutant
%enzUs= sol.x(rxnIndx);
%enzUsage = enzymeUsage_FVA(model, Enzymes);

model = setParam(model, 'ub', {'acOUT','glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn'}, [1.3, 1.7, 0,0,0,0]);
Mutant= getMutant(model, modifications, true);
