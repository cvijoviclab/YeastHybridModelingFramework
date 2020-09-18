function [rxnIds,rxn_Names] = get_original_rxnIDs(rxns,rxnNames)
%Get original rxn IDs
rxns = strrep(rxns,'GAL10','GAL_X');
rxnIds = unique(rxns);
rxnIds = strrep(rxnIds,'_REV','');
rxnIds = strrep(rxnIds,'arm_','');
rxnIds = regexprep(rxnIds,'No(\d)','');
rxnIds = unique(rxnIds(~contains(rxnIds,'prot_')));
rxnIds = strtrim(rxnIds);
rxn_Names = [];
if nargin>1
    %Get original rxn names
    rxn_Names = unique(rxnNames);
    rxn_Names = strrep(rxn_Names,'(reversible)','');
    rxn_Names = strrep(rxn_Names,'(arm)','');
    rxn_Names = regexprep(rxn_Names,'(No(\d))','');
    rxn_Names = unique(rxn_Names(~contains(rxn_Names,'prot_')));
    rxn_Names = strtrim(rxn_Names);
end
end