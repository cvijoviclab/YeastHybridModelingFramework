% regulationSensitivityAnalisis
% Created by Linnea, last edited: 2020-04-19

clc; 
close all;
clear;

%add path to functions
addpath ('sourceCode');
%load models
load('../data/Bool.mat');
% load('data/tempModel.mat')
% model = tempModel;
load('../../../models/reduced_ecYeast_fermentation.mat');
model=ecModel_ferm;

%% RunBoolean
settings.activeCrosstalk = [1 1 1 1 1]; 
%specify knockout(s)
settings.knockouts = {};
%manually change initial nutrient levels
settings.gluc = [0]; %initilizing glucose setting
settings.nitr = [1 1]; %sequence of nitrogen concentration
%change in Kcat and enzyme usage
%load('../data/temporarykcats.mat');
%settings.kcat = kcat; %put in value from the sensitivity analysis
%settings.enzymeUse = 0.0; %put in value from the analysis
%settings.enzymeUse = 0.05;

[Bool, TFAct, EnzAct] = runBool(Bool,settings);
[GeneExpr1] = expression(TFAct.End, model);  
settings.gluc = [1];
[Bool, TFAct, EnzAct] = runBool(Bool,settings);
[GeneExpr2] = expression(TFAct.End, model);  

regulatedGenes = GeneExpr1.enzGenes(GeneExpr1.rank~=0 | GeneExpr2.rank~=0);

for i = 1:length(regulatedGenes)
    FBAmodel(regulatedGenes(i), GeneExpr1.rank(i), GeneExpr1.rank(i),model);
end
regulatedEnzNames = GeneExpr1.enzNames(GeneExpr1.rank~=0 | GeneExpr2.rank~=0);