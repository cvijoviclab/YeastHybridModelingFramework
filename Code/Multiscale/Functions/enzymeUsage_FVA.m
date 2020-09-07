function FVAtable = enzymeUsage_FVA(model,enzymes)
% enzymeUsage_FVA
%   Perf
%
%   model           Original model (either GEM or ecGEM)
%   modifications   [nx3] Cell array with the desired modifications. Column
%                   (1) contains strings with each of the individual gene IDs
%                   to modify. Column (2) indicates the modification
%                   action; 0 for deletions, 1 for overexpressions acting on
%                   enzyme usages bounds, 2 for overexpressionbs acting on
%                   Kcat values, and 3 for modification of native enzymes
%                   activity.
%                   for enzyme activity modification (substitution
%                   endogenous<-heterologous enzyme). Column (3) Over-expression 
%                    factor, 0 for deletions, a given number for overexpressions 
%                    and the new Kcat for heterologous expression (1/s)
%   enzUsage        Basal usage level for the enzyme to modufy (for
%                   overexpresion [mmol/gDw h]
%   message         Flag, true if a success message should be displayed
%
%   mutantModel     Mutant model with new constraints
%
%   Usage: mutantModel = getMutant(model,modifications,enzUsage,message)
%
%   Ivan Domenzain.     Last edited 2019-11-13

%Get parsimonious protein usages
tempModel  = model;
pool_indx  = find(strcmpi(model.rxns,'prot_pool_exchange'));
prot_indxs = find(contains(model.rxnNames,'prot_'));
prot_indxs = prot_indxs(1:(end-1));
tempModel = setParam(tempModel, 'obj',pool_indx,-1);
sol       = solveLP(tempModel,1);
%initialize variables
ranges    = [];
minUsgs   = [];
maxUsgs   = [];
enz_pUsgs = [];
if ~isempty(sol.x)
    pUsgs = sol.x(prot_indxs);
    %Loop through all the provided enzymes
    for i=1:length(enzymes)
        if ~isempty(enzymes{i})
            rxnIndx = find(contains(model.rxnNames,enzymes{i}));
            enzIndx = strcmpi(model.enzymes,enzymes{i});
            model = setParam(model, 'obj', rxnIndx, -1);
            sol   = solveLP(model);
            if ~isempty(sol.f)
                minFlux = sol.x(rxnIndx); 
                model   = setParam(model, 'obj', rxnIndx, +1);
                sol     = solveLP(model);
                if ~isempty(sol.f)
                   disp(['Ready with rxn #' num2str(i)])
                   maxFlux = sol.x(rxnIndx); 
                else
                   maxFlux = nan; 
                end
            else
                minFlux = nan;
                maxFlux = nan; 
            end
            ranges    = [ranges; (maxFlux-minFlux)];
            minUsgs   = [minUsgs; minFlux]; 
            maxUsgs   = [maxUsgs; maxFlux];
            enz_pUsgs = [enz_pUsgs;pUsgs(enzIndx)];
        else
            ranges    = [ranges; 0];
            minUsgs   = [minUsgs; 0]; 
            maxUsgs   = [maxUsgs; 0];
            enz_pUsgs = [enz_pUsgs;0];
        end
    end
else
    disp('Model is not feasible')
end
varNamesT = {'enzymes' 'ranges' 'minU' 'maxU' 'pU'};
FVAtable  = table(enzymes,ranges,minUsgs,maxUsgs,enz_pUsgs,'VariableNames', varNamesT);
end