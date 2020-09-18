function results = getNetFluxes(rxnIDs,rxn_Names,fluxes,original_IDs)
rxnIDs = strrep(rxnIDs,'GAL10','GAL_X');
fluxes(fluxes==0) = 1E-15;
%Get net fluxes for all reactions
rxns      = [];
rxnNames  = [];
netFluxes = [];
ratios    = [];
futileFlux = [];
for i=1:length(original_IDs)
    ID        = original_IDs{i};
    fluxRatio = 0;
    % Check if i-th reaction is reversible
    revFlag = any(find(contains(rxnIDs,ID) & contains(rxnIDs,'_REV')));
    % Map rxn indexes in flux dist table
    rxnIdxs = rxnMapping(ID,rxnIDs,revFlag);
    %calculate flux ratios and net flux in case the reaction is reversible
    if ~isempty(rxnIdxs)
        if length(rxnIdxs) == 2
            fluxRatio = fluxes(rxnIdxs(2))/fluxes(rxnIdxs(1));
            netFlux   = fluxes(rxnIdxs(1)) -fluxes(rxnIdxs(2));
        elseif length(rxnIdxs) == 1
            netFlux   = fluxes(rxnIdxs(1));
        else
            disp(ID)
            rxnIDs(rxnIdxs)
        end
        rxnName = rxn_Names(rxnIdxs(1));
        rxnName = strrep(rxnName,'(arm)','');
        rxnName = strrep(rxnName,'No1','');
        rxnName = strtrim(strrep(rxnName,'()',''));
        rxns      = [rxns;original_IDs(i)];
        ratios    = [ratios;fluxRatio];
        netFluxes = [netFluxes; netFlux];
        rxnNames  = [rxnNames; rxnName];
        futileFlux = [futileFlux;abs(max(fluxes(rxnIdxs)))-abs(netFlux)];
    else
        disp(original_IDs{i})
    end
end
rxns    = strrep(rxns,'GAL_X','GAL10');
results = table(rxns,rxnNames,netFluxes,ratios,futileFlux);
end