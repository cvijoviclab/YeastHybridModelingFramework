% multiscaleModelSimulation2 
% This script is the next version of the multiscale simulation model.
% Created by Linnea, last edited: 2020-04-19

clc; 
close all;
clear;
tolerance = 1E-12;
perturbation = 1+(1E-6);
%add path to functions
addpath ('sourceCode');
%load models
load('data/Bool.mat');
%load('data/tempModel.mat')
% model = tempModel;
load('../../models/reduced_ecYeast_fermentation.mat');
model = ecModel_ferm;
% Create folder and directery for results
[path] = createResultsFolder('multiscaleRegulated', yyyymmdd(datetime));
% Make settings
%manually change crosstalks (0 turns off, 1 turns on -> in order to work 
%as described in literature: all crosstalks must be turned on). 
settings.activeCrosstalk = [1 1 1 1 1]; 
% specify knockout(s)
settings.knockouts = {};
%manually change initial nutrient levels
settings.gluc = 0; %initilizing glucose setting
settings.nitr = [1 1]; %sequence of nitrogen concentration
% change in Kcat and enzyme usage
load('data/temporarykcats.mat');
settings.kcat = kcat; %put in value from the sensitivity analysis
%settings.enzymeUse = 0.0; %put in value from the analysis
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
%Get relevant rxn and met indexes for simulations and analysis
idxs = get_model_Idxs(model);
exchIndexes = [idxs(5);idxs(6);idxs(7);idxs(8)];
cSource = idxs(5);
GAMmet  = idxs(3);
GAMrxn  = idxs(4);
bioRxn  = idxs(1);
%Relevant pathways for further exploration
CC_paths = {'Glycolysis' 'TCA' 'pentose phosphate' 'Oxidative Phosphorylation' 'Anaerobic excretion'};
% reactions relevant for Boolean model glucose signals
signal_idxs  = [7 45 60 61 62 63 64 69 158]; 
%Prepare model for simulations
model_pool = prepare_ecModel(model);
%Initialize tables for saving all flux distributions and enzyme usage profiles
formulas       = constructEquations(model_pool);
variables      = {'rxns' 'rxnNames' 'formulas' 'grRules' 'subSystems'};
fluxDist_reg   = table(model_pool.rxns,model_pool.rxnNames,formulas,model_pool.grRules,model_pool.subSystems,'VariableNames',variables);
fluxDist       = fluxDist_reg;
pathways       = mapEnzymeSubSystems(model_pool.enzymes,model);
variables      = {'enzymes' 'enzNames' 'genes' 'MWs' 'subSystems'};
Enz_usages_reg = table(model_pool.enzymes,model_pool.enzNames,model_pool.enzGenes,model_pool.MWs,pathways);
Enz_usages     = Enz_usages_reg;
fermentation = false;
conBool      = [];
pathW_protB  = [];
SSE      = [];
results  = [];
error    = [];
j=0;
i=0;
% Simulations!!
mkdir('../../results/ECCs')
for subopt_growth  = 0:(0.4/80):0.4
    tempModel      = model_pool;
    exchangeVector = zeros(length(exchIndexes),1);
    %Fix suboptimal growth rate
    tempModel = setParam(tempModel,'lb',bioRxn,subopt_growth);
    tempModel = setParam(tempModel,'ub',bioRxn,(1+tolerance)*subopt_growth);
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
    solution = solveLP(tempModel,1);
    initSol  = solution.x;
    j = j+1;    
    if ~isempty(solution.f)
        if subopt_growth == 0.1 | subopt_growth == 0.4
            tempModel = setParam(tempModel,'lb',cSource,0);
            tempModel = setParam(tempModel,'ub',cSource,1000);
            [ECC_glc,~] = get_ECCs(tempModel,60,perturbation,tolerance,'descend');
            if subopt_growth == 0.1
                writetable(ECC_glc,'../../results/ECCs/ECC_glc_respiration.txt','Delimiter','\t','QuoteStrings',false);
            else
                writetable(ECC_glc,'../../results/ECCs/ECC_glc_fermentation.txt','Delimiter','\t','QuoteStrings',false);
            end
        end
        %Fix optimal glucose uptake and set protein usage as objective to
        %minimize
        tempModel = setParam(tempModel,'ub',cSource,1.15*solution.x(cSource));
        %Constrain based on Boolean model
        tempModel = constrainFBA2(GeneExpr, PTM, tempModel, settings);
        %Solve!
        tempModel = setParam(tempModel,'lb',cSource,0);
        tempModel = setParam(tempModel,'ub',cSource,1000);
        solution  = solveLP(tempModel,1);
        if ~isempty(solution.x)
            if subopt_growth == 0.1 | subopt_growth == 0.4
                [ECC_glc,~] = get_ECCs(tempModel,60,perturbation,tolerance,'descend');
                if subopt_growth == 0.1
                    writetable(ECC_glc,'../../results/ECCs/ECC_glc_respiration_reg.txt','Delimiter','\t','QuoteStrings',false);
                else
                    writetable(ECC_glc,'../../results/ECCs/ECC_glc_fermentation_reg.txt','Delimiter','\t','QuoteStrings',false);
                end
            end
            exchangeVector = solution.x(exchIndexes);         
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
    conBool= [conBool, solution.x(signal_idxs)];
    %Save flux distributions and enzyme usages (in percentage of protein mass)
    str = ['fluxDist_reg.D_' strrep(num2str(subopt_growth),'0.','')  '=solution.x;'];
    str2 = ['fluxDist.D_' strrep(num2str(subopt_growth),'0.','')  '= initSol;'];
    eval(str)
    eval(str2)
    T  = topUsedEnzymes(solution.x,tempModel,{''},length(tempModel.enzymes));
    T2 = topUsedEnzymes(initSol,tempModel,{''},length(tempModel.enzymes));
    str = ['Enz_usages_reg.D_' strrep(num2str(subopt_growth),'0.','') '=T.Usages_;'];
    str2 = ['Enz_usages.D_' strrep(num2str(subopt_growth),'0.','') '=T2.Usages_;'];
    eval(str)
    eval(str2)
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
    disp(' ')
end
% Write results
toRemove = find(contains(fluxDist_reg.rxns,'prot_'));
fluxDist_reg(toRemove,:) = [];
toRemove = find(contains(fluxDist.rxns,'prot_'));
fluxDist(toRemove,:) = [];
writetable(fluxDist_reg,'../../results/fluxDist_reg_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
writetable(Enz_usages_reg,'../../results/enzUsages_reg_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
writetable(fluxDist,'../../results/fluxDist_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
writetable(Enz_usages,'../../results/enzUsages_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
%Display the median relative error for exchange fluxes prediction across
%dilution rates
SSE = mean(SSE);
disp(['The median error is: ' num2str(SSE)]);
error = [error;SSE];
%Plot results
figure()
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
Reactions  = model_pool.rxnNames(signal_idxs);
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