function met_turnovers(resultsPath)
%load model
model      = open('../../models/ecYeast_CCEM.mat');
model      = model.model;
%open flux distribution files (regulated)
fluxTable = readtable([resultsPath '/fluxDist_reg_reducedYeast.txt'],'delimiter','\t');
[original_IDs,original_rxnNames] = get_original_rxnIDs(fluxTable.rxns,fluxTable.rxnNames);
results = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_1,original_IDs);
original_IDs = strrep(results.rxns,'GAL10','GAL_X');
%Get original S matrix
[Smat,mets,metNames] = get_original_Smat(original_IDs,model);
resp_reg = get_met_turnover(Smat,results.netFluxes);
results  = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_4,original_IDs);
ferm_reg = get_met_turnover(Smat,results.netFluxes);
%Repeat for unregulated flux distributions
fluxTable = readtable([resultsPath '/fluxDist_reducedYeast.txt'],'delimiter','\t');
[original_IDs,original_rxnNames] = get_original_rxnIDs(fluxTable.rxns,fluxTable.rxnNames);
results = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_1,original_IDs);
original_IDs = strrep(results.rxns,'GAL10','GAL_X');
[Smat,mets,metNames] = get_original_Smat(original_IDs,model);
resp_unreg = get_met_turnover(Smat,results.netFluxes);
results    = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_4,original_IDs);
ferm_unreg = get_met_turnover(Smat,results.netFluxes);
%Get fold-changes
FC_resp = (resp_reg+1E-15)./(resp_unreg+1E-15);
FC_ferm = (ferm_reg+1E-15)./(ferm_unreg+1E-15);
FC_reg  = (ferm_reg+1E-15)./(resp_reg+1E-15);
FC_unreg  = (ferm_unreg+1E-15)./(resp_unreg+1E-15);
%Write results
turnovers = table(mets,metNames,resp_unreg,resp_reg,ferm_unreg,ferm_reg,FC_resp,FC_ferm,FC_reg,FC_unreg);
writetable(turnovers,[resultsPath '/met_turnovers.txt'],'delimiter','\t','QuoteStrings',false);
end
%compute met turnover numbers 
function turnoverNumbers = get_met_turnover(S,fluxDist)
[m,~] = size(S);
turnoverNumbers = zeros(m,1);
for i=1:m
    turnoverNumbers(i) = 0.5*(sum(abs(S(i,:)'.*fluxDist)));     
end
end