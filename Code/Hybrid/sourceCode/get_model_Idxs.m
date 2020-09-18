function idxs = get_model_Idxs(model)
%Get GAM related indexes
bioRxn  = find(contains(model.rxns,'GROWTH'));
protIdx = find(strcmpi(model.rxns,'prot_pool_exchange'));
GAMmet  = find(strcmpi(model.metNames,'Maintainance for growth'));
GAMrxn  = find(model.S(GAMmet,:)>0);
%find relevant exchange reactions
cSource  = find(strcmpi(model.rxnNames,'uptake of glucose'));
oxyIndex = find(strcmpi(model.rxnNames,'Uptake of O2'));
CO2Index = find(strcmpi(model.rxnNames,'production of co2'));
ethIndex = find(strcmpi(model.rxnNames,'production of ethanol'));
aceIndex = find(strcmpi(model.rxnNames,'production of acetate'));
glyIndex = find(strcmpi(model.rxnNames,'production of glycerol'));
idxs = [bioRxn protIdx GAMmet GAMrxn cSource oxyIndex CO2Index ethIndex aceIndex glyIndex];
end