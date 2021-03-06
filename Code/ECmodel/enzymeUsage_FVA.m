function FVAtable = enzymeUsage_FVA(model,enzymes,verbose)
% enzymeUsage_FVA
%
%   Ivan Domenzain.     Last edited 2020-04-14
if nargin<3
    verbose = true;
    if nargin<2
        enzymes = model.enzymes;
    end
end
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
                    if verbose
                        disp(['Ready with rxn #' num2str(i)])
                    end
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