function analyze_NetFluxes(resultsPath)
model      = load('../../models/reduced_ecYeast_fermentation.mat');
model      = model.ecModel_ferm;
fluxTable  = readtable('../../results/fluxDist_reg_reducedYeast.txt','delimiter','\t');
[rxnIds,~] = get_original_rxnIDs(fluxTable.rxns);
resp_reg   = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_1,rxnIds);
ferm_reg   = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_4,rxnIds);

fluxTable  = readtable('../../results/fluxDist_reducedYeast.txt','delimiter','\t');
[rxnIds,~] = get_original_rxnIDs(fluxTable.rxns);
resp_Unreg = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_1,rxnIds);
ferm_Unreg = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_4,rxnIds);

variables = {'rxns' 'rxnNames' 'netFlux' 'ratio' 'futileFlux' 'netFlux_reg' 'ratio_reg' 'futileFlux_reg'};
results = table(resp_Unreg.rxns,resp_Unreg.rxnNames,resp_Unreg.netFluxes,resp_Unreg.ratios,resp_Unreg.futileFlux,resp_reg.netFluxes,resp_reg.ratios,resp_reg.futileFlux,'VariableNames',variables);
writetable(results,'../../results/netFlux_respiration.txt','delimiter','\t','QuoteStrings',false);
variables = {'rxns' 'rxnNames' 'netFlux' 'ratio' 'futileFlux' 'netFlux_reg' 'ratio_reg' 'futileFlux_reg'};
results = table(ferm_Unreg.rxns,ferm_Unreg.rxnNames,ferm_Unreg.netFluxes,ferm_Unreg.ratios,ferm_Unreg.futileFlux,ferm_reg.netFluxes,ferm_reg.ratios,ferm_reg.futileFlux,'VariableNames',variables);
writetable(results,'../../results/netFlux_fermentation.txt','delimiter','\t','QuoteStrings',false);
end