%Ivan Domenzain. Last modified 2020-04-22
system('git clone https://github.com/SysBioChalmers/GECKO');
load('../../models/reduced_ecYeast_fermentation.mat')
model = ecModel_ferm;
gIndex = find(contains(model.rxnNames,'Production of biomass'));
Pindex = find(contains(model.rxnNames,'prot_pool_'));
%Set growth rate as a maximization objective
model = setParam(model,'obj',gIndex(1),1);
%Block "reversible" formation of biomass
model = setParam(model,'ub',gIndex(2),0);
%Get enzyme control coefficients on max growth
[limitations,~] = find_All_Limitations(model);
%Create a .txt file as an output
writetable(limitations,'../../results/growthRate_ECC.txt','Delimiter','\t','QuoteStrings',false)
%Get enzyme fluxes (in terms of protein mass percentage)
cd GECKO/geckomat/kcat_sensitivity_analysis
%bi-level optimization
solution = solveLP(model);
model = setParam(model,'lb',gIndex(1),0.999999*solution.x(gIndex(1)));
model = setParam(model,'obj',Pindex,-1);
solution = solveLP(model,1);
enzTable = topUsedEnzymes(solution.x,model,{'max_gRate'},'reducedYeast_fermentation',false,length(model.enzymes));
cd ../../..
writetable(enzTable,'../../results/fermentation_enzUsages.txt','Delimiter','\t','QuoteStrings',false)
%Analyze fluxes and enzyme usage distributions across growth rates
fluxTable = readtable('../../results/fluxDist_reducedYeast.txt');
enzTable = readtable('../../results/enzUsages_reducedYeast.txt');
enzMassTable = readtable('../../results/enzMassUsages_reducedYeast.txt');
pathways = {'Glycolysis' 'Oxidative Phosphorylation' 'TCA' 'pentose phosphate' 'Anaerobic excretion'};

%Get flux distributions metrics
variables = {'D_rate' 'nZeros' 'maxflux' 'minflux' 'meanflux' 'medianflux' 'stddev' ...
             'fermentation' 'glyc_meanFlux' 'oxPhos_meanFlux' 'TCA_meanFlux' ...
             'PPP_meanFlux'};
dist_metrics = getDistributionMetrics(fluxTable,variables,pathways);
writetable(dist_metrics,'../../results/fluxDist_metrics.txt','Delimiter','\t','QuoteStrings',false);
%Get enzyme usage distributions metrics (percentage)
variables = {'D_rate' 'nZeros' 'maxUsage' 'minUsage' 'meanUsage' 'medianUsage' 'stddev' ...
             'fermentation' 'glyc_meanUsage' 'oxPhos_meanUsage' 'TCA_meanUsage' ...
             'PPP_meanUsage' 'ana_meanUsage' 'glyc_burden' 'oxPhos_burden' 'TCA_burden' 'PPP_burden' 'ana_burden'};
dist_metrics1 = getDistributionMetrics(enzTable,variables,pathways,model);
writetable(dist_metrics1,'../../results/enzUsage_metrics.txt','Delimiter','\t','QuoteStrings',false);
%Get enzyme usage distributions metrics (mass)
variables = {'D_rate' 'nZeros' 'maxUsage' 'minUsage' 'meanUsage' 'medianUsage' 'stddev' ...
             'fermentation' 'glyc_meanUsage' 'oxPhos_meanUsage' 'TCA_meanUsage' ...
             'PPP_meanUsage' 'ana_meanUsage' 'glyc_burden' 'oxPhos_burden' 'TCA_burden' 'PPP_burden' 'ana_burden'};
dist_metrics = getDistributionMetrics(enzMassTable,variables,pathways,model);
writetable(dist_metrics,'../../results/enzMassUsage_metrics.txt','Delimiter','\t','QuoteStrings',false);
