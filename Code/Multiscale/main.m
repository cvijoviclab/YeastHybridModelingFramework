%% main script to run multiscale model
%written by: Julia M�nch
%date: 2019-11-25
%description: 
%   script to simulate nutrient flux dependency on growthrate while 
%   modifying turnover number of enzymes contained in the metabolic model 
%   to simulated changed gene expression due to variation in nutrient
%   availability. Simulation of knockouts of signalling pathway components
%   and variation in crosstalk activity can be manually applied; saves gene
%   Expression and Enzyme Activity changes as txt files
%
%contains functions:
%   1. initialization
%   2. runMultiscale 
%       contains 3. knockouts
%                4. Expression
%                5. reachSteadyState
%                   contains 6. crosstalk
%                            7. activityConverter
%                            8. TFtargets                    
%                9. BoolToFBA2
%                   contains 10. getMutant
%                11. plotfluxdistr
%                12. plotfluxes
%                

%clear command window, close figures, clear variables
clc; 
close all;
clear;

%add path to functions
addpath ('Functions');
load('RequiredDocs/reduced_ecYeast_batch.mat');
model_in = ecModel_batch;

%% Initialize Model %%%
[Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = initialization ();

%% manually change crosstalks (0 turns off, 1 turns on -> in order to work as described in literature: all crosstalks must be turned on)
activeCrosstalk = [1 1 1 1 1]; 

%% specify knockout
knockouts = {};

%% Produce data for different nutrient levels %%%

%manually change initial nutrient levels
gluc = [0 1]; %sequence of glucose concentration
nitr = [0 1]; %sequence of nitrogen concentration

%% percentage of change in Kcat
perc = 0.03;
foldername = ['changeKcat', strrep(num2str(perc),'.', ','), '/'];
path = ['Data/', foldername];
       
%% create folders if they don't already exist
mkdir Data
mkdir('Data', foldername)

%% run multiscale model
[ExpressionChange, EnzActChange, o2, co2, eth, ac, glcin, growth] = runMultiscale(ecModel_batch, knockouts, activeCrosstalk, gluc, nitr, Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, perc, path);

%% save
writetable(ExpressionChange, [path, 'ExpressionChange.txt'], 'Delimiter', '\t');
writetable(EnzActChange, [path, 'EnzActChange.txt'], 'Delimiter', '\t');
