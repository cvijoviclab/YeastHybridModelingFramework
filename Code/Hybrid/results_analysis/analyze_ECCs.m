function analyze_ECCs(condition,resultsPath)
% KEq_table = readtable('../../data/kEqTable.txt','delimiter','\t');
 model = load('../../../models/reduced_ecYeast_fermentation.mat');
 model = model.ecModel_ferm;
 ECC_table_reg = readtable(['../../../results/ECCs/ECC_glc_' condition '_reg.txt'],'delimiter','\t');
 ECC_table     = readtable(['../../../results/ECCs/ECC_glc_' condition '.txt'],'delimiter','\t');
% fluxTable = readtable('../../results/fluxDist_reducedYeast.txt','delimiter','\t');
% fluxes = fluxTable.D_1;
% rxnIDs = KEq_table.rxnID;
% ratios = [];
% ECCs = [];
% ECCs_reg = [];
% KEqs = [];
% reactions = [];
% rxnNames  = [];
% netFluxes = [];
 subSystems = [];
% %normalize ECCs
% % sum(ECC_table_reg.ECC)
% % ECC_table_reg.ECC = ECC_table_reg.ECC/sum(ECC_table_reg.ECC);
%  ECC_table_reg.rxnID = strrep(ECC_table_reg.rxnID,'GAL10','GAL_X');
%  fluxTable.rxns = strrep(fluxTable.rxns,'GAL10','GAL_X');
% 
%  % sum(ECC_table.ECC)
%  %ECC_table.ECC = ECC_table.ECC/sum(ECC_table.ECC);
%  ECC_table.rxnID = strrep(ECC_table.rxnID,'GAL10','GAL_X');
% %For each reaction in the Keq table, search its corresponding entries in
% %the fluxes and ECCs tables
% model.rxns = strrep(model.rxns,'GAL10','GAL_X');
% for i=1:length(rxnIDs)
%     ID = rxnIDs{i};
%     fluxRatio = 0;
%     if strcmpi(ID,'TPH')
%         ID = 'GLD';
%     end
%     if strcmpi(ID,'GAL10')
%         ID = 'GAL_X';
%     end
%     % Check if i-th reaction is reversible
%     revFlag = any(find(contains(fluxTable.rxns,ID) & contains(fluxTable.rxns,'_REV')));
%     % Map rxn indexes in flux dist table
%     rxnIdxs = rxnMapping(ID,model.rxns,revFlag);
%     %calculate flux ratios and net flux in case the reaction is reversible
%     if length(rxnIdxs) == 2
%         %disp(model.rxnNames(rxnIdxs(2)))
%         %disp(model.rxnNames(rxnIdxs(1)))
%         fluxRatio = fluxes(rxnIdxs(2))/fluxes(rxnIdxs(1));
%         netFlux   = fluxes(rxnIdxs(1)) -fluxes(rxnIdxs(2));
%     elseif length(rxnIdxs) == 1
%         %disp(model.rxnNames(rxnIdxs(1)))
%         netFlux   = fluxes(rxnIdxs(1));
%     else
%         disp(rxnIDs{i})
%     end
%     rxnName = model.rxnNames{rxnIdxs(1)};
%     subSystem = model.subSystems(rxnIdxs(1));
%     rxnName = strrep(rxnName,'arm','');
%     rxnName = strrep(rxnName,'No1','');
%     rxnName = strtrim(strrep(rxnName,'()',''));
%     %Search reaction in ECCs table
%     idxs = find(contains(ECC_table.rxnID,rxnIDs(i)));
%     idxs_reg = find(contains(ECC_table_reg.rxnID,rxnIDs(i)));
%     if ~isempty(idxs) & ~isempty(idxs_reg)
%         reactions = [reactions;rxnIDs(i)];
%         ratios = [ratios;fluxRatio];
%         netFluxes = [netFluxes; netFlux];
%         KEqs = [KEqs;KEq_table.Value(i)];
%         rxnNames  = [rxnNames; {rxnName}];
%         %Save the max ECCs associated to the reaction ID
%         ECCs = [ECCs;max(ECC_table.ECC(idxs))];
%         ECCs_reg = [ECCs_reg;max(ECC_table_reg.ECC(idxs_reg))];
%     end
% end
% results = table(reactions,rxnNames,netFluxes,ratios,ECCs,ECCs_reg);
% %Write table

rxnIDs = unique([ECC_table_reg.rxnID;ECC_table.rxnID]);
ECCs   = [];
ECCs_reg = [];
for i=1:length(rxnIDs)
    idx = find(contains(ECC_table_reg.rxnID,rxnIDs(i)));
    if ~isempty(idx)
        ECCs_reg = [ECCs_reg; ECC_table_reg.ECC(idx)]; 
    else
        ECCs_reg = [ECCs_reg; 0]; 
    end
    idx = find(contains(ECC_table.rxnID,rxnIDs(i)));
    if ~isempty(idx)
        ECCs = [ECCs; ECC_table.ECC(idx)]; 
    else
        ECCs = [ECCs; 0]; 
    end
    idx = find(contains(model.rxns,rxnIDs(i)),1);
    subSystems = [subSystems;model.subSystems(idx)];

end
%rxnIDs   = regexprep(rxnIDs,'No(\d)','');
results2 = table(rxnIDs,subSystems,ECCs,ECCs_reg);
results2 = sortrows(results2,'ECCs','descend');
writetable(results2,[resultsPath '/ECCs_' condition '.txt'],'delimiter','\t','QuoteStrings',false);
end