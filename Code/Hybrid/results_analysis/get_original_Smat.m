function [Smat,mets,metNames] = get_original_Smat(original_IDs,model)
model.rxns = strrep(model.rxns,'GAL10','GAL_X');
Smat = zeros(length(model.mets),length(original_IDs));
mets = model.mets;
for i=1:length(original_IDs)
    rxnID   = original_IDs{i};
    rxnIdxs = find(contains(model.rxns,rxnID) & ~contains(model.rxns,'_REV') & ~contains(model.rxns,[rxnID 'Rev']));
    if length(rxnIdxs) > 1
        rxnIdxs = rxnIdxs(contains(model.rxns(rxnIdxs),'arm_') | contains(model.rxns(rxnIdxs),'No1'));
        rxnVector = model.S(:,rxnIdxs(1)) + model.S(:,rxnIdxs(2)); 
    else
    	rxnVector = model.S(:,rxnIdxs);
    end
    Smat(:,i) = rxnVector;
end
toKeep   = (~contains(mets,'pmet_') & ~startsWith(mets,'prot_')); 
mets     = mets(toKeep);
metNames = model.metNames(toKeep);
Smat     = Smat(toKeep,:);  
end