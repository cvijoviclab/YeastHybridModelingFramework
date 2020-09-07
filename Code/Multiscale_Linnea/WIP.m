%WIP

%% Load model
load('data/reduced_ecYeast_batch.mat');
model = ecModel_batch;

%% Load experimental data
fID  = fopen('data/chemostatData.txt');
data = textscan(fID,'%f %f %f %f %f %f %f %f %f %f','Delimiter','\t','HeaderLines',1);
Drate = data{1};
qO2   = data{3};
qCO2  = data{4};
qGluc = data{5};
qEtOH = data{6};

%% print biomass pseudoreaction formulation
bioRxn = find(contains(model.rxns,'GROWTH'));
printModel(model,bioRxn)
%It seems that the growth-associated ATP maintenance is encompassed 
%in the pseudo metabolite "Maintainance for growth", save the index of the
%GAM reaction/mets for further fitting of GAM
GAMmet = find(strcmpi(model.metNames,'Maintainance for growth'));
GAMrxn = find(model.S(GAMmet,:)>0);
printModel(model,GAMrxn)
%Prevent uncoupling (from Avlant's code)
model = setParam(model, 'ub', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
model = setParam(model, 'lb', {'ShuttleXRev', 'ATPTransportRev', 'OAC1Rev', 'PyrTransRev', 'CAT2Rev'}, 0);
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
%Close all uptakes except for medium components
model_pool.ub(find(model_pool.ub==Inf)) = 1000;
uptakes    = find(contains(model_pool.rxnNames,'Uptake'));
model_pool = setParam(model_pool,'ub',uptakes,0);
medium = {'glucose' 'O2' 'Phosphate' 'H2O'};
for component=medium
    rxn   = ['Uptake of ' component{1}];
    index = find(strcmpi(model.rxnNames,rxn));
    model_pool = setParam(model_pool,'ub',index,1000);
end
%Block active uptake of glucose
model_pool = setParam(model_pool,'ub',{'gluActiveIn'},0);
%Limit acetate production
model_pool   = setParam(model_pool,'ub',aceIndex,0.6);
%Set NGAM
model_pool.lb((contains(model.rxns,'ATPX'))) = 0.7;
%Incorporate Avlant's manually curated data:
[model_pool,~] = manualModifications(model_pool);
%Ptot was fitted according to exp data on exchange fluxes
model_pool.ub(end) = 0.099;
fermentation       = false;
j=0;
i=0;
SSE      = [];
results  = [];
error    = [];
%% Initialize tables for saving all flux distributions and enzyme usage
%profiles
fluxDist = table();
fluxDist.rxns       = model_pool.rxns;
fluxDist.formulas   = constructEquations(model_pool);
fluxDist.grRules    = model_pool.grRules;
fluxDist.subSystems = model_pool.subSystems;
Enz_usages          = table();
Enz_usages.enzymes  = model_pool.enzymes;
%Simulations!!
for subopt_growth = 0:(0.4/80):0.4
    tempModel = model_pool;
    %Fix suboptimal growth rate
    tempModel = setParam(tempModel,'lb',bioRxn,0.9999*subopt_growth);
    %Set minimization of glucose uptake as objective function
    tempModel = setParam(tempModel,'obj',cSource,-1);
    %Solve FBA problem
    solution  = solveLP(tempModel);
    %Set GAM
    if ~fermentation
        GAM = -45+(45-60)*j/56;
    else
    	GAM = -60+(60-30)*(j-56)/(80-56);
    end
    tempModel.S(GAMmet,bioRxn) = GAM;
    j   = j+1;
    if ~isempty(solution.f)
            %Fix optimal glucose uptake and set protein usage as objective to 
            %minimize
            tempModel = setParam(tempModel,'lb',cSource,0.9999*solution.x(cSource));
            tempModel = setParam(tempModel,'ub',cSource,1.0001*solution.x(cSource));
            %Set minimization of total protein usage as objective function
            tempModel = setParam(tempModel,'obj',{'prot_pool_exchange'},-1);
            %Solve!
            solution  = solveLP(tempModel,1);
        if ~isempty(solution.f)
            exchangeVector = solution.x(exchIndexes);
            FVAtable = enzymeUsage_FVA(tempModel,tempModel.enzymes,false);
            str = ['enzUsage_VarAnalysis_D_' strrep(num2str(subopt_growth),'0.','') '.txt'];
            writetable(FVAtable,['results/WIP/' str],'Delimiter','\t','QuoteStrings',false);
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
writetable(fluxDist,'results/WIP/fluxDist_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
writetable(Enz_usages,'results/WIP/enzUsages_reducedYeast.txt','Delimiter','\t','QuoteStrings',false);
%Display the median relative error for exchange fluxes prediction across
%dilution rates
SSE = mean(SSE);
disp(['The median error is: ' num2str(SSE)]);
error = [error;SSE];
%Plot results
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
hold off