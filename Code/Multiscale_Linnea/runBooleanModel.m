%runBooleanModel

%written by: Linnea Österberg and Julia Münch
%date: 2020-04-02
%description: 
%contains functions:
% 1. createResultsFolder

clc; 
close all;
clear;

%add path to functions
addpath ('sourceCode');
%load model
load('data/Bool.mat');
% Create folder and directery for results
[path] = createResultsFolder('booleanSimulations', yyyymmdd(datetime));
%% Make settings
%manually change crosstalks (0 turns off, 1 turns on -> in order to work 
%as described in literature: all crosstalks must be turned on). 
settings.activeCrosstalk = [1 1 1 1 1]; 
% specify knockout(s)
settings.knockouts = {};
%manually change initial nutrient levels
settings.gluc = [0 1]; %sequence of glucose concentration
settings.nitr = [0 1]; %sequence of nitrogen concentration
%% Run first itteration of Bool 0 glucose 1 nitrogen. 
[Bool, TFAct, EnzAct] = runBool2(Bool,settings,path);
pathway = {dir([path, 'Activity/Transitions/' '*.txt']).name};

%loop over pathways
for i = 1:length(pathway)
    
    %create folder for figures
    mkdir([path, 'Figures/'], 'Transitions/')
    mkdir([path, 'Figures/'], 'Heatmaps/')
    mkdir([path, 'Figures/'], 'Change/')
    
    plotActivity(pathway{i}, path, nutr)
    plotTransitions(pathway{i}, path)
    barPlot(pathway{i}, path)
        
end