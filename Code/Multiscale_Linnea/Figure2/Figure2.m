% Figure2
% Desctiption: This scripts analyzes the enzyme usage from the multiscale 
% simulations of the deletion strains.

%add path to functions
addpath ('functions');
%% PROTEIN ALLOCATION
% Load unregulated and regulated simulation data and relative abundance.
respiration.regulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulated20200612nr2/enzUsage_VarAnalysis_D_1.txt');
respiration.TORdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedTORdeletion20200616/enzUsage_VarAnalysis_D_1.txt');
respiration.SNF1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedSNF1deletion20200616/enzUsage_VarAnalysis_D_1.txt');
respiration.PKAdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedPKAdeletion20200616/enzUsage_VarAnalysis_D_1.txt');
respiration.REG1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedREG1deletion20200616/enzUsage_VarAnalysis_D_1.txt');

fermentation.regulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulated20200612nr2/enzUsage_VarAnalysis_D_4.txt');
fermentation.TORdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedTORdeletion20200616/enzUsage_VarAnalysis_D_4.txt');
fermentation.SNF1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedSNF1deletion20200616/enzUsage_VarAnalysis_D_4.txt');
fermentation.PKAdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedPKAdeletion20200616/enzUsage_VarAnalysis_D_395.txt');
fermentation.REG1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedREG1deletion20200616/enzUsage_VarAnalysis_D_4.txt');


%% FLUX ANALYSIS
flux.regulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulated20200612nr2/fluxDist_reducedYeast.txt');
flux.TORdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedTORdeletion20200616/fluxDist_reducedYeast.txt');
flux.SNF1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedSNF1deletion20200616/fluxDist_reducedYeast.txt');
flux.PKAdeletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedPKAdeletion20200616/fluxDist_reducedYeast.txt');
flux.REG1deletion=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/Figure2/Simulations/multiscaleRegulatedREG1deletion20200616/fluxDist_reducedYeast.txt');

