%% Boolean Model

%written by: Julia M�nch
%date: 2019-10-28
%description: 
%   script to produce data for different nutrient conditions that are provided 
%   in vector form for glucose and nitrogen:
%       no glucose/nitrogen present (0)
%       glucose/nitrogen present (1)
%   -crosstalk activities between pathways that can
%   manually be turned off and on
%   -knockouts can be simulated by providing name of component to turn off

%contains functions:
%   1. initialization
%   2. knockouts
%   3. reachSteadyState 
%       contains 4. activityConverter
%                5. saveTransitions
%                6. crosstalk
%                7. TFtargets
%                8. getRanks
%   9. plotActivity
%   10. plotTransitions
%   11. barPlot

%clear command window, close figures, clear variables
clc; 
close all;
clear;

%add path to functions
addpath ('Functions');

%% Initialize Model %%%

[Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = initialization ();

%% Produce data for different nutrient levels %%%

%manually change nutrient levels
gluc = [0 1]; %sequence of glucose concentration
nitr = [0 1]; %sequence of nitrogen concentration

%determine knockout
knockouts = {};

foldername = ['glu', num2str(gluc(1)), 'nit', num2str(nitr(1)), '-glu', num2str(gluc(2)), 'nit', num2str(nitr(2)), '/'];
path = ['../../results/Boolean/', foldername];
%path = ['Data/Validation/', foldername]; for validation knockouts

%manually change crosstalks (0 turns off, 1 turns on -> in order to work as described in literature: all crosstalks must be turned on)
activeCrosstalk = [1 1 1 1 1]; 

nutr = {[num2str(gluc(1)), '|', num2str(nitr(1))], [num2str(gluc(2)), '|', num2str(nitr(2))]};

%loop over glucose levels
for i = 1:length(gluc)
    
    glucLevel = gluc(i);
    nitrLevel = nitr(i);
       
    %create folders if they don't already exist
    mkdir(path)
    mkdir(path, 'Transitions/')
    mkdir(path, 'Figures/')
    mkdir(path, 'Activity/')
    mkdir([path, 'Activity/'], 'Transitions/')
    
    disp(['Glucose level: ', num2str(gluc(i)), ', Nitrogen level: ', num2str(nitr(i))]);
    
    %run the boolean model to create txt files
    
    %create knockout (first after system has reached steady state for
    %inital conditions)
    if i > 1
        [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = knockout(Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, knockouts);        
    end
    
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, transTFAct] = ...
        reachSteadyState(Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, glucLevel, nitrLevel, path, activeCrosstalk);

end


%% Create Figures

% create array containing all pathways of interest

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