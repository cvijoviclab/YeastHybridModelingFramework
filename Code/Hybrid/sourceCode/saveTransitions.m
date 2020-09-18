function[transitions] = ...
    saveTransitions(iteration, MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, transitions)

%written by: Julia Münch
%date: 2019-11-07
%description: this function returns the states' transitions when changing
%the nutrient conditions
%arguments: 
%   1. number of iteration i
%   2.-7. matrix representing state of each pathway before iteration i
%   8. vector containing matrices of for each pathway that are extended
%   after each iteration
%returns: vector containing matrices of each pathway describing components' states after each iteration   

%% rename and save initial conditions
if (iteration == 1)
    for i = 1:length(transitions)
        transitions{i}.Properties.VariableNames{2} = ['pres ', num2str(iteration-1)];
    end
    for i = 2:5
        transitions{i}.Properties.VariableNames{3} = ['phos ', num2str(iteration-1)];
    end
    transitions{2}.Properties.VariableNames{4} = ['act ', num2str(iteration-1)];
    
% rename and append new states after iteration
else
    y = {MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld};
    for i = 1:length(y)
        y{i}.Properties.VariableNames{2} = ['pres ', num2str(iteration-1)];
    end
    for i = 2:5
        y{i}.Properties.VariableNames{3} = ['phos ', num2str(iteration-1)];
    end
    y{2}.Properties.VariableNames{4} = ['act ', num2str(iteration-1)];
    
    %append new states
    for i = 1:length(transitions)
        transitions{i} = join(transitions{i}, y{i}, 'keys', 'Name');
    end
    
end