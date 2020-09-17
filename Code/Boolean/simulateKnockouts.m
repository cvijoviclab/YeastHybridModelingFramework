%% Boolean Model running knockouts

%written by: Linnea and Julia
%date: 2020-04-16
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
gluc = [1 0]; %sequence of glucose concentration
nitr = [1 0]; %sequence of nitrogen concentration

%determine knockout sequence
knockouts = {"WT","Snf1","Reg1"',["Tpk1","Tpk2","Tpk3"], ["Tor1","Tor2"]};

foldername = ['KO_glu', num2str(gluc(2)), 'nit', num2str(nitr(2)), '/'];
path = ['../../results/Boolean/', foldername];
%path = ['Data/Validation/', foldername]; %for validation knockouts

%manually change crosstalks (0 turns off, 1 turns on -> in order to work as described in literature: all crosstalks must be turned on)
activeCrosstalk = [1 1 1 1 1]; 

nutr = {[num2str(gluc(1)), '|', num2str(nitr(1))]};

%loop over glucose levels and knockouts
for i = 1:length(knockouts)
    
    glucLevel = gluc(2);
    nitrLevel = nitr(2);
       
    %create folders if they don't already exist
    mkdir(path,[num2str(i),'/'])
    mkdir(path, ['Transitions/',num2str(i),'/'])
    mkdir(path, 'Figures/')
    mkdir(path, 'Activity/')
    mkdir([path, 'Activity/'], ['Transitions/',num2str(i),'/'])
    
    disp(['Glucose level: ', num2str(gluc(1)), ', Nitrogen level: ', num2str(nitr(1))]);
    
    %run the boolean model to create txt files
    
    %create knockout (first after system has reached steady state for
    %inital conditions)
 
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = initialization ();   
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, transTFAct] = ...
        reachSteadyState2(i,Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, gluc(1), nitr(1), path, activeCrosstalk);
    if i > 1
        [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = knockout2(Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, knockouts{1,i});        
    end
    
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, transTFAct] = ...
        reachSteadyState2(i,Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, glucLevel, nitrLevel, path, activeCrosstalk);

end


%% Create Figures

% create array containing all pathways of interest

pathway = {dir([path, 'Activity/Transitions/1/' '*.txt']).name};

%loop over pathways
for i = 1:length(pathway)
    
    %create folder for figures
    mkdir([path, 'Figures/'], 'Transitions/')
    mkdir([path, 'Figures/'], 'Heatmaps/')
    mkdir([path, 'Figures/'], 'Change/')
    
    plotActivity2(pathway{i}, path, nutr, knockouts)
    %plotTransitions(pathway{i}, path)
    %barPlot(pathway{i}, path)
        
end