function []=FBAmodel(regulatedGene, RankRespiration, RankFermentation, model)
%% Create folder and directery for results
[path] = createResultsFolder('regulationAnalysis', yyyymmdd(datetime),regulatedGene );

%% Make settings
%manually change crosstalks (0 turns off, 1 turns on -> in order to work 
%as described in literature: all crosstalks must be turned on). 
%settings.activeCrosstalk = [1 1 1 1 1]; 
% specify knockout(s)
%settings.knockouts = {};
%manually change initial nutrient levels
%settings.gluc = [0]; %initilizing glucose setting
%settings.nitr = [1 1]; %sequence of nitrogen concentration
% change in Kcat and enzyme usage
%load('../data/temporarykcats.mat');
%settings.kcat = kcat; %put in value from the sensitivity analysis
%settings.enzymeUse = 0.0; %put in value from the analysis
settings.enzymeUse = 0.05; %put in value from the analysis in respiration the fermentation value is changes in "connectBool2FBA"
%% Regulated enzyme
% geneIdx = find(contains(model.enzGenes, regulatedGene));
% Enzymes = model.enzymes(geneIdx);
% enzUsageList = enzymeUsage_FVA(model, Enzymes);
% enzyme  = model.enzymes(geneIdx);
% enzRxn = find(contains(model.rxnNames,enzyme));
% % if enzyme is down regulated upper bounds (enzUsageList(i,4) are adjusted down
% if RankRespiration < 0
%     model.ub(enzRxn) = table2array(enzUsageList(:,4)) - ...
%         ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
% % if enzyme is up regulated lower bounds (enzUsageList(i,3) are adjusted up
% elseif RankRespiration > 0
%     model.lb(enzRxn) = table2array(enzUsageList(:,3)) + ...
%         ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
% end

%% Load experimental data
fID  = fopen('data/chemostatData.txt');
data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1);
Drate = data{1};
qO2   = data{3};
qCO2  = data{4};
qGluc = data{5};
qEtOH = data{6};

%% Preper for simulations
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
exchIndexes = [cSource;oxyIndex;CO2Index;ethIndex];

%Set growth maximization as an objective function
model_pool = setParam(model,'obj',bioRxn,1);

solution   = solveLP(model_pool,1);
if ~isempty(solution.f)
    maxGrowth = solution.x(bioRxn);
    disp(['The maximum growth rate for the model on glucose minimal media is: ' num2str(maxGrowth) ' [g/gDw h]'])
end

%Ptot was fitted according to exp data on exchange fluxes
model_pool.ub(end) = 0.038;
fermentation       = false;
j=0;
i=0;
SSE      = [];
results  = [];
error    = [];

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
%% Simulations!!
for subopt_growth = 0:(0.4/80):0.4
    tempModel = model_pool;
    %Fix suboptimal growth rate
    tempModel = setParam(tempModel,'lb',bioRxn,0.9999*subopt_growth);
    %Set minimization of glucose uptake as objective function
    tempModel = setParam(tempModel,'obj',cSource,-1);
    %Set GAM
    upper = 30;
    fermL = 26;
    respL = 18;

    if ~fermentation
        GAM = -respL+(respL-upper)*j/57;
    else
        GAM = -upper+(upper-fermL)*(j-58)/(80-58);
    end
    tempModel.S(GAMmet,bioRxn) = GAM;
    
    %Change Boolean
    %% Regulated enzyme
    if ~fermentation
        geneIdx = find(contains(tempModel.enzGenes, regulatedGene));
        Enzymes = tempModel.enzymes(geneIdx);
        enzUsageList = enzymeUsage_FVA(tempModel, Enzymes);
        enzyme  = tempModel.enzymes(geneIdx);
        enzRxn = find(contains(tempModel.rxnNames,enzyme));
        % if enzyme is down regulated upper bounds (enzUsageList(i,4) are adjusted down
        if RankRespiration < 0
            tempModel.ub(enzRxn) = table2array(enzUsageList(:,4)) - ...
                ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
            % if enzyme is up regulated lower bounds (enzUsageList(i,3) are adjusted up
        elseif RankRespiration > 0
            tempModel.lb(enzRxn) = table2array(enzUsageList(:,3)) + ...
                ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
        end
    else
        geneIdx = find(contains(tempModel.enzGenes, regulatedGene));
        Enzymes = tempModel.enzymes(geneIdx);
        enzUsageList = enzymeUsage_FVA(tempModel, Enzymes);
        enzyme  = tempModel.enzymes(geneIdx);
        enzRxn = find(contains(tempModel.rxnNames,enzyme));
        % if enzyme is down regulated upper bounds (enzUsageList(i,4) are adjusted down
        if RankFermentation < 0
            tempModel.ub(enzRxn) = table2array(enzUsageList(:,4)) - ...
                ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
            % if enzyme is up regulated lower bounds (enzUsageList(i,3) are adjusted up
        elseif RankFermentation > 0
            tempModel.lb(enzRxn) = table2array(enzUsageList(:,3)) + ...
                ((table2array(enzUsageList(:,4))-table2array(enzUsageList(:,3)))* settings.enzymeUse);
        end
    end
%     [Bool,settings] = connectFBA2Bool(Bool,settings,conBool,i);
%     [Bool, TFAct, EnzAct] = runBool(Bool,settings);
%     [GeneExpr] = expression(TFAct.End, tempModel);  
%     [PTM] = postTranscriptionalModifications2(EnzAct.End, tempModel);       
    
    %Solve FBA problem
    solution  = solveLP(tempModel);
    j   = j+1;
    
    if ~isempty(solution.f)
            %Fix optimal glucose uptake and set protein usage as objective to 
            %minimize
            tempModel = setParam(tempModel,'lb',cSource,0.9999*solution.x(cSource));
            tempModel = setParam(tempModel,'ub',cSource,1.0001*solution.x(cSource));
            %Set minimization of total protein usage as objective function
            tempModel = setParam(tempModel,'obj',{'prot_pool_exchange'},-1);
            %Constrain based on Boolean model
            %tempModel = constrainFBA2(GeneExpr, PTM, tempModel, settings);
            %tempModel = constrainFBA2(GeneExpr, PTM, tempModel, settings);
            %Solve!
            solution  = solveLP(tempModel,1);
        if ~isempty(solution.f)
            exchangeVector = solution.x(exchIndexes);
            FVAtable = enzymeUsage_FVA(tempModel,tempModel.enzymes,false);
            str = ['enzUsage_VarAnalysis_D_' strrep(num2str(subopt_growth),'0.','') '.txt'];
            writetable(FVAtable,[path str],'Delimiter','\t','QuoteStrings',false);
        end
    else
        exchangeVector = zeros(length(exchIndexes),1);
    end
    if ~fermentation
        disp(['Dilution rate = ' num2str(subopt_growth) ': Respiration'])
        if exchangeVector(4)>1E-2 
            disp(['The critical dilution rate is: ' num2str(subopt_growth) ' [1/h]'])
            fermentation = true;
        end
    else
        disp(['Dilution rate = ' num2str(subopt_growth) ': Fermentation'])
    end
    exchangeVector(exchangeVector==0) = 1E-6;
    newRow  = [subopt_growth, exchangeVector'];
    results = [results; newRow];
    %Save the connection points with Boolean model
    if ~isempty(solution.f)
        conBool= [conBool, solution.x(RxnIdx)];
    else
        conBool= [conBool, zeros(length(conBool),1)];
    end
    %Save flux distributions and enzyme usages (in percentage of protein
    %mass)
    str = ['fluxDist.D_' strrep(num2str(subopt_growth),'0.','')  '=solution.x;'];
    eval(str)
    T = topUsedEnzymes(solution.x,tempModel,{''},length(tempModel.enzymes));
    str = ['Enz_usages.D_' strrep(num2str(subopt_growth),'0.','') '=T.Usages_;'];
    eval(str)
    %compare with experimental data
    %Search subopt_growth in Drate exp data
    [~,dataIndex] = ismember(subopt_growth,Drate);
    if dataIndex>0
        Drate = data{1};
        expExchFlux = [data{5}(dataIndex);data{3}(dataIndex);data{4}(dataIndex);data{6}(dataIndex)];
        expExchFlux(expExchFlux==0) = 1E-6;
        %Get median relative error for exchange fluxes prediction for each
        %Drate
        SSE = [SSE; mean(abs(expExchFlux-exchangeVector)./(expExchFlux))];
    end
    i = i+1;
end
%% Write results
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
% figure()
% names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol'};
% for i=1:(length(exchIndexes))
%     plot(results(:,1),results(:,i+1),'LineWidth',3)
%     hold on
% end
% %Add experimental data points
% vector = [5 3 4 6];
% for i=vector
%     scatter(Drate,data{i})
%     hold on
% end
% legend(names)
% xlabel('Dilution rate [1/h]','FontSize',18)
% ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
% xlim([0 max(results(:,1))])
% hold off
save([path,'settings.mat'],'-struct','settings');
% % Make plot for Boolean input
% figure()
% Reactions  = model_pool.rxnNames(RxnIdx);
% plot(0:(0.4/80):0.4,conBool','LineWidth',3)
% legend(Reactions)
% xlabel('Dilution rate [1/h]','FontSize',18)
% ylabel('Fluxes [mmol/gDw h]','FontSize',18)
% conBool2 = table(Reactions,conBool);
% save([path,'conBool.mat'],'conBool2');