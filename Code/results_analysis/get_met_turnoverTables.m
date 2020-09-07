function get_met_turnoverTables(sourceFile,resultsFile,resultsPath)
%load model
model      = open('../../models/ecYeast_CCEM.mat');
model      = model.model;
%open flux distribution files (regulated)
fluxTable = readtable([resultsPath '/' sourceFile],'delimiter','\t');
[original_IDs,~] = get_original_rxnIDs(fluxTable.rxns,fluxTable.rxnNames);
netFlux = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.D_1,original_IDs);
original_IDs = strrep(netFlux.rxns,'GAL10','GAL_X');
%Get original S matrix
[Smat,mets,metNames] = get_original_Smat(original_IDs,model);
%Get all dilution rates from fluxTable 
Drates = fluxTable.Properties.VariableNames;
Drates = Drates(contains(Drates,'D_'));
turnoverNumbers = table(mets,metNames);
%compute met turnovers for each dilution rate
for Drate = Drates
    %Get net fluxes from flux distribution
    eval(['netFlux = getNetFluxes(fluxTable.rxns,fluxTable.rxnNames,fluxTable.' Drate{1} ',original_IDs);'])
    %Get met turnovers
    turnovers = get_met_turnover(Smat,netFlux.netFluxes);
    eval(['turnoverNumbers.' Drate{1} '=turnovers;'])
end
writetable(turnoverNumbers,[resultsPath '/' resultsFile],'delimiter','\t','QuoteStrings',false);
end

function turnoverNumbers = get_met_turnover(S,fluxDist)
[m,~] = size(S);
turnoverNumbers = zeros(m,1);
for i=1:m
    turnoverNumbers(i) = 0.5*(sum(abs(S(i,:)'.*fluxDist)));     
end
end