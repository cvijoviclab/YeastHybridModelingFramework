% multiscaleModelSimulation2 
% This script is the next version of the multiscale simulation model.
% Created by Linnea, last edited: 2020-04-19

clc; 
close all;
clear;
%add path to functions
addpath ('sourceCode');
%load models
load('data/Bool.mat');
%load('data/tempModel.mat')
% model = tempModel;
load('../../models/reduced_ecYeast_fermentation.mat');
model=ecModel_ferm;
% Create folder and directery for results
[path] = createResultsFolder('multiscaleRegulatedDeletionTOR', yyyymmdd(datetime));
% Make settings
%manually change crosstalks (0 turns off, 1 turns on -> in order to work 
%as described in literature: all crosstalks must be turned on). 
settings.activeCrosstalk = [1 1 1 1 1]; 
% specify knockout(s)
settings.knockouts = {'Tor1','Tor2'};
%manually change initial nutrient levels
settings.gluc = [0]; %initilizing glucose setting
settings.nitr = [1 1]; %sequence of nitrogen concentration
% change in Kcat and enzyme usage
load('data/temporarykcats.mat');
settings.kcat = kcat; %put in value from the sensitivity analysis
%settings.enzymeUse = 0.0; 
settings.enzymeUse = (0.05); %put in value from the analysis in respiration the fermentation value is changes in "connectBool2FBA"
[Bool] = knockout(Bool, settings.knockouts);
% Load experimental data
fID  = fopen('data/chemostatData.txt');
data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1);
Drate = data{1};
qO2   = data{3};
qCO2  = data{4};
qGluc = data{5};
qEtOH = data{6};
%Get GAM related indexes
bioRxn = find(contains(model.rxns,'GROWTH'));
GAMmet = find(strcmpi(model.metNames,'Maintainance for growth'));
GAMrxn = find(model.S(GAMmet,:)>0);
%find relevant exchange reactions
cSource     = find(strcmpi(model.rxnNames,'uptake of glucose'));
oxyIndex    = find(strcmpi(model.rxnNames,'Uptake of O2'));
CO2Index    = find(strcmpi(model.rxnNames,'production of co2'));
ethIndex    = find(strcmpi(model.rxnNames,'production of ethanol'));
aceIndex    = find(strcmpi(model.rxnNames,'production of acetate'));
glyIndex    = find(strcmpi(model.rxnNames,'production of glycerol'));
protIdx     = find(strcmpi(model.rxns,'prot_pool_exchange'));
exchIndexes = [cSource;oxyIndex;CO2Index;ethIndex;aceIndex];
exchIndexes = [cSource;oxyIndex;CO2Index;ethIndex];
%Relevant pathways for further exploration
CC_paths = {'Glycolysis' 'TCA' 'pentose phosphate' 'Oxidative Phosphorylation' 'Anaerobic excretion'};
%Add exchange reaction for cytosolic pyruvate
pyr = find(strcmpi(model.metNames,'pyruvate'),1);
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
[grRules, rxnGeneMat]    = standardizeGrRules(model,true);
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
fermentation       = false;
SSE      = [];
results  = [];
error    = [];
pathW_protB = [];
%Initialize tables for saving all flux distributions and enzyme usage profiles
fluxDist = table();
fluxDist.rxns       = model_pool.rxns;
fluxDist.rxnNames   = model_pool.rxnNames;
fluxDist.formulas   = constructEquations(model_pool);
fluxDist.grRules    = model_pool.grRules;
fluxDist.subSystems = model_pool.subSystems;
Enz_usages          = table();
Enz_usages.enzymes  = model_pool.enzymes;
pathways            = mapEnzymeSubSystems(model_pool.enzymes,model);
Enz_usages.subSystems = pathways;
Enz_massUsages      = Enz_usages;

j=0;
i=0;
conBool= [];
RxnIdx  = [7 45 60 61 62 63 64 69 158]; % reactions relevant for Boolean model glucose signals
model_pool.S(GAMmet,bioRxn) = -18; % Start value for GAM
model_pool = setParam(model_pool,'lb',cSource,-1000); % Unconstrain glucose uptake
model_pool = setParam(model_pool,'ub',cSource,1000); % Unconstrain glucose uptake
% Simulations!!
for subopt_growth  = 0:(0.4/80):0.4
    tempModel      = model_pool;
    exchangeVector = zeros(length(exchIndexes),1);
    protDemand     = 0;
    %Fix suboptimal growth rate
    tempModel = setParam(tempModel,'lb',bioRxn,0.999999*subopt_growth);
    tempModel = setParam(tempModel,'lb',bioRxn,1.000001*subopt_growth);
    %Set GAM
    upper = 30;
    fermL = 25;
    respL = 18;
    tempModel = setParam(tempModel,'obj',cSource,-1); 
    if ~fermentation
        GAM = -respL+(respL-upper)*j/57;
    else
        GAM = -upper+(upper-fermL)*(j-58)/(80-58);
    end
    tempModel.S(GAMmet,bioRxn) = GAM;
    %Change Boolean
    [Bool,settings] = connectFBA2Bool(Bool,settings,conBool,i);
    [Bool, TFAct, EnzAct] = runBool(Bool,settings);
    [GeneExpr] = expression(TFAct.End, tempModel);  
    [PTM] = postTranscriptionalModifications2(EnzAct.End, tempModel);       
    %Solve FBA problem
    solution = solveLP(tempModel);
    j = j+1;    
    if ~isempty(solution.f)
        %Fix optimal glucose uptake and set protein usage as objective to 
        %minimize
        %tempModel = setParam(tempModel,'ub',cSource,1.25*solution.x(cSource));
        tempModel = setParam(tempModel,'ub',cSource,1.15*solution.x(cSource));
        %Constrain based on Boolean model
        tempModel = constrainFBA2(GeneExpr, PTM, tempModel, settings);
        %Solve!
        tempModel = setParam(tempModel,'lb',cSource,0);
        tempModel = setParam(tempModel,'ub',cSource,1000);
        solution  = solveLP(tempModel,1);
        if ~isempty(solution.x)
            exchangeVector = solution.x(exchIndexes);
            protDemand     = solution.x(protIdx);
        else
            solution.x = zeros(length(tempModel.rxns),1);
        end
    end
    
    if ~fermentation
        disp(['Dilution rate = ' num2str(subopt_growth) ': Respiration'])
        if exchangeVector(4)>1E-2 
            fermentation = true;
            disp(['The critical dilution rate is: ' num2str(subopt_growth) ' [1/h]'])
        end
    else
        disp(['Dilution rate = ' num2str(subopt_growth) ': Fermentation'])
    end
    
    exchangeVector(exchangeVector==0) = 1E-6;
    newRow  = [subopt_growth, exchangeVector'];
    results = [results; newRow];
    %Save the connection points with Boolean model
    if isempty(solution.f)
       solution.x = zeros(length(tempModel.rxns),1);
    end
    conBool= [conBool, solution.x(RxnIdx)];
    %Save flux distributions and enzyme usages (in percentage of protein mass)
    str = ['fluxDist.D_' strrep(num2str(subopt_growth),'0.','')  '=solution.x;'];
    eval(str)
    T = topUsedEnzymes(solution.x,tempModel,{''},length(tempModel.enzymes));
    str = ['Enz_usages.D_' strrep(num2str(subopt_growth),'0.','') '=T.Usages_;'];
    eval(str)
    % Save FVA
    FVAtable = enzymeUsage_FVA(tempModel,tempModel.enzymes,false);
    str = ['enzUsage_VarAnalysis_D_' strrep(num2str(subopt_growth),'0.','') '.txt'];
    writetable(FVAtable,[path str],'Delimiter','\t','QuoteStrings',false);
    %Get total protein burden per pathway
    subSystems = mapEnzymeSubSystems(tempModel.enzymes,tempModel);
    burden     = [];
    for k=1:length(CC_paths)
        path_idx = find(contains(subSystems,CC_paths{k}));
        burden   = [burden, sum(T.Usages_(path_idx))];
    end
    pathW_protB = [pathW_protB; burden];
    %compare with experimental data
    %Search subopt_growth in Drate exp data
    [~,dataIndex] = ismember(subopt_growth,Drate);
    if dataIndex>0
        Drate = data{1};
        %expExchFlux = [data{5}(dataIndex);data{3}(dataIndex);data{4}(dataIndex);data{6}(dataIndex);data{7}(dataIndex)];
        expExchFlux = [data{5}(dataIndex);data{3}(dataIndex);data{4}(dataIndex);data{6}(dataIndex)];
        expExchFlux(expExchFlux==0) = 1E-6;
        %Get median relative error for exchange fluxes prediction for each
        %Drate
        SSE = [SSE; mean(abs(expExchFlux-exchangeVector)./(expExchFlux))];
    end
    i = i+1;
end
% Write results
toRemove = find(contains(fluxDist.rxns,'prot_'));
fluxDist(toRemove,:) = [];
writetable(fluxDist,[path 'fluxDist_reducedYeast.txt'],'Delimiter','\t','QuoteStrings',false);
writetable(Enz_usages,[path 'enzUsages_reducedYeast.txt'],'Delimiter','\t','QuoteStrings',false);
%Display the median relative error for exchange fluxes prediction across
%dilution rates
SSE = mean(SSE);
disp(['The median error is: ' num2str(SSE)]);
error = [error;SSE];
%Plot results
figure()
names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol' 'Acetate'};
names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol'};

for i=1:(length(exchIndexes))
    plot(results(:,1),results(:,i+1),'LineWidth',3)
    hold on
end
%Add experimental data points
vector = [5 3 4 6];
for i=vector
    scatter(Drate,data{i})
    hold on
end
legend(names)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
xlim([0 max(results(:,1))])
ylim([0 25])

hold off
save([path,'settings.mat'],'-struct','settings');
% Make plot for Boolean input
figure()
Reactions  = model_pool.rxnNames(RxnIdx);
plot(0:(0.4/80):0.4,conBool','LineWidth',3)
legend(Reactions)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Fluxes [mmol/gDw h]','FontSize',18)
conBool2 = table(Reactions,conBool);
save([path,'conBool.mat'],'conBool2');

%Plot protein burden by pathways
figure()
for i=1:(length(CC_paths))
    plot(results(:,1),pathW_protB(:,i),'LineWidth',3)
    hold on
end
legend(CC_paths)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Protein burden [mmol/gDw]','FontSize',18)
xlim([0 max(results(:,1))])
ylim([0 0.03])