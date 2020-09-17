function model = prepare_ecModel(model)
%prepare_ecModel
% 
% Prepares the ecModel for integration with regulatory boolean network
%
idxs = get_model_Idxs(model);
protIdx  = idxs(2);
cSource  = idxs(5);
bioRxn   = idxs(1);
aceIndex = idxs(9);
glyIndex = idxs(10);
GAMmet   = idxs(3);
%Add exchange reaction for cytosolic pyruvate
rxnFormula = 'pyruvate[c] => ';
rxnsToAdd.equations = {rxnFormula};
rxnsToAdd.rxns      = {'PyrOut'};
rxnsToAdd.grRules   = {''};
rxnsToAdd.rxnNames  = {'Pyruvate production'};
rxnsToAdd.c  = 0;
rxnsToAdd.lb = 0;
rxnsToAdd.ub = 0.055;
% Introduce changes to the model
model = addRxns(model,rxnsToAdd,3);
% %Standardize gene related fields
[grRules, rxnGeneMat] = standardizeGrRules(model,true);
model.grRules    = grRules;
model.rxnGeneMat = rxnGeneMat;
%Prepare for simulations
%Set growth maximization as an objective function
model_pool = setParam(model,'obj',bioRxn,1);
%Set physiological levels of byproducts secretion as upper bounds
model_pool = setParam(model_pool,'ub',aceIndex,0.682);
model_pool = setParam(model_pool,'ub',glyIndex,0.165);
solution   = solveLP(model_pool,1);
if ~isempty(solution.f)
    maxGrowth = solution.x(bioRxn);
    disp(['The maximum growth rate for the model on glucose minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
end

%Ptot was fitted according to exp data on exchange fluxes
model_pool.ub(protIdx) = 0.038;
model_pool.S(GAMmet,bioRxn) = -18; % Start value for GAM
model_pool = setParam(model_pool,'lb',cSource,-1000); % Unconstrain glucose uptake
model_pool = setParam(model_pool,'ub',cSource,1000); % Unconstrain glucose uptake
%Save curated model
model = model_pool;
save('../../models/ecYeast_CCEM.mat','model');
end