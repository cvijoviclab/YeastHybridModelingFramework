%load flux distributions for all dilution rates
fluxTable  = readtable('../../../results/fluxDist_reducedYeast.txt','delimiter','\t');
fluxTable_reg  = readtable('../../../results/fluxDist_reg_reducedYeast.txt','delimiter','\t');
%Load proteomics data for both conditions
data_resp = readtable('../../../data/respirationDataset1.txt','delimiter','\t','HeaderLines',4);
data_ferm = readtable('../../../data/fermentationDataset1.txt','delimiter','\t','HeaderLines',3);
%reformat datasets
data_ferm = data_ferm(:,[1 3]);
data_resp.NaN(isnan(data_resp.NaN))   = 0;
data_ferm.Var3(isnan(data_ferm.Var3)) = 0;
data_resp.presence = logical(data_resp.NaN);
data_ferm.presence = logical(data_ferm.Var3);
%Get all reactions with isoenzymes
arm_idxs = find(startsWith(fluxTable.rxns,'arm_'));
armIds   = strrep(fluxTable.rxns(arm_idxs),'arm_','');
armIds   = strrep(armIds,'_REV','');
armIds   = unique(armIds);
%extract flux distributions for D=0.1 and D = 0.4 from both ecModel and
%hybrid model results
fluxMAT  = [fluxTable.D_1,fluxTable_reg.D_1,fluxTable.D_4,fluxTable_reg.D_4];
resultsT = table();
%Iterate through all arm reactions
for i=1:length(armIds)
    armRxn = armIds(i);
    %Search all isoenzymes for i-th arm reaction
    iso_idxs = find(startsWith(fluxTable.rxns,armRxn) & ~contains(fluxTable.rxns,'_REVNo'));
    %get gene associations for all isoenzymes for i-th reaction
    iso_grRules = fluxTable.grRules(iso_idxs);
    %compare to proteomics data
    protPresence = zeros(numel(iso_idxs),2);
    %for each enzymatic complex catalizing each of the isoenzyme reactions
    %search for presence of its genes in the respiration and fermentation
    %proteomics datasets
    for j=1:length(iso_grRules)
        genes = iso_grRules{j};
        genes = strsplit(genes,' and ');
        for gene=genes
            try
                protPresence(j,1) = protPresence(j,1) + data_resp.presence(strcmpi(data_resp.Var1,gene{1}));
            catch
                protPresence(j,1) = protPresence(j,1) + 0;
            end
            
            try
            	protPresence(j,2) = protPresence(j,2) + data_ferm.presence(strcmpi(data_ferm.Var1,gene{1}));
            catch
                protPresence(j,2) = protPresence(j,2) + 0;
            end
        end
    end
    %Check which isoenzyme reactions are carrying flux
    carryFlux = logical(fluxMAT(iso_idxs,:));
    %in case that no flux is carried by forward reactions check also
    %reversible directions (they have an impact on the protein burden too).
    iso_idxs_rev = find(startsWith(fluxTable.rxns,armRxn) & contains(fluxTable.rxns,'_REVNo'));
    if ~isempty(iso_idxs_rev)
        revFlux = logical(fluxMAT(iso_idxs_rev,:));%num2cell(logical(fluxMAT(iso_idxs_rev,:)));
        carryFlux = logical(carryFlux+revFlux);
    end
    %append to results
    protPresence = num2cell(logical(protPresence));
    carryFlux    = num2cell(carryFlux);
    newBlock     = [repelem(armRxn,numel(iso_idxs),1),strrep(fluxTable.rxns(iso_idxs),armRxn,''),iso_grRules,carryFlux,protPresence];
    resultsT     = [resultsT;newBlock];
end
resultsT.Properties.VariableNames = {'rxn' 'isoenzymes' 'genes' 'ecModel_resp' 'hybrid_resp' 'ecModel_ferm' 'hybrid_ferm' 'protein_resp' 'protein_ferm'};
%Get confusion matrices and associated metrics
[confMat_ecM_R,ecM_R] = getConfusionMetrics(resultsT.protein_resp,resultsT.ecModel_resp);
[confMat_ecM_F,ecM_F] = getConfusionMetrics(resultsT.protein_ferm,resultsT.ecModel_ferm);
[confMat_hyb_R,hyb_R] = getConfusionMetrics(resultsT.protein_resp,resultsT.hybrid_resp);
[confMat_hyb_F,hyb_F] = getConfusionMetrics(resultsT.protein_ferm,resultsT.hybrid_ferm);
%save confusion matrices file
confMat = [confMat_ecM_R,confMat_hyb_R,confMat_ecM_F,confMat_hyb_F];
confMat = num2cell(confMat);
confMat = cell2table(confMat);
confMat.Properties.VariableNames = {'obs_T_ecM_R' 'obs_F_ecM_R' 'obs_T_hyb_R' 'obs_F_hyb_R' ...
                                    'obs_T_ecM_F' 'obs_F_ecM_F' 'obs_T_hyb_F' 'obs_F_hyb_F'};
confMat.Properties.RowNames = {'predictions_true' 'predictions_false'};  
writetable(confMat,'../../../results/isoenzymes_confusionMatrix.txt','delimiter','\t','QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true);

%arrange and save results from confusion metrics as tables
summ_Table = table(ecM_R,hyb_R,ecM_F,hyb_F);
summ_Table.Properties.RowNames      = {'sensitivity' 'specificity' 'precision' 'accuracy' 'F1_score' 'MCC' 'FMI' 'informedness' 'markedness'};
writetable(summ_Table,'../../../results/isoenzymes_utilization_summary.txt','delimiter','\t','QuoteStrings',false,'WriteRowNames',true,'WriteVariableNames',true);
%save overall results for isoenzymes_utilization as .txt file
resultsT = resultsT(:,[1 2 3 8 4 5 9 6 7]);
writetable(resultsT,'../../../results/isoenzymes_utilization.txt','delimiter','\t','QuoteStrings',false,'WriteRowNames',false,'WriteVariableNames',true);
