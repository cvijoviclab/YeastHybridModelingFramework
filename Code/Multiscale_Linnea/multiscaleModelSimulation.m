%% main script to run multiscaleModelSimulation
%written by: Linnea Österberg and Julia Münch
%date: 2020-04-02
%description: 
%   script to simulate nutrient flux dependency on growthrate while 
%   modifying turnover number of enzymes contained in the metabolic model 
%   to simulated changed gene expression due to variation in nutrient
%   availability. Simulation of knockouts of signalling pathway components
%   and variation in crosstalk activity can be manually applied
%
%contains functions:
%   1. createResultFolder
%   2. runBool
%       contains 3. activityConverter
%                4. boolRules
%                5. crosstalk
%                6. TFtargets
%   7. knockouts
%   8. expression
%   9. postTranscriptionalModifications
%   10. constrainFBA
%        contains 11. getMutant
%   12. plotfluxdistr
%   13. plotfluxes
%                

%clear command window, close figures, clear variables
clc; 
close all;
clear;

%add path to functions
addpath ('sourceCode');
%load models
load('data/Bool.mat');
load('data/reduced_ecYeast_batch.mat');
model = ecModel_batch;
%% Create folder and directery for results
[path] = createResultsFolder('changeKcat', yyyymmdd(datetime));

%% Make settings
%manually change crosstalks (0 turns off, 1 turns on -> in order to work 
%as described in literature: all crosstalks must be turned on). 
settings.activeCrosstalk = [1 1 1 1 1]; 
% specify knockout(s)
settings.knockouts = {};
%manually change initial nutrient levels
settings.gluc = [0 1]; %sequence of glucose concentration
settings.nitr = [1 1]; %sequence of nitrogen concentration
% percentage of change in Kcat
settings.perc = 0.03; %at some point try perc[0.03; 0.06;0.09] for low mediuma and high effect.

%save([path,'settings.mat'],'-struct','settings');

%% Run first itteration of Bool 0 glucose 1 nitrogen. 
[Bool, TFAct, EnzAct] = runBool(Bool,settings);

if numel(settings.knockouts) ~= 0
    [Bool] = knockout(BoolIn, settings.knockouts);
    [Bool, TFAct, EnzAct] = runBool(Bool.settings);
end
[GeneExpr] = expression(TFAct.End, model); % values are between -4 and 4 if deviding into differnt levels later. 
[GeneExpr] = postTranscriptionalModifications(EnzAct.End, GeneExpr);% values are between -4 and 5 if deviding into differnt levels later.This is when PTM gives 2 in rank.  
%% Run EC-FBA

mutant=constrainFBA(GeneExpr, model, settings.perc);


%% Notes

x_val = 0:0.005:5;
parameters.eth=[];
parameters.o2=[];
parameters.ac=[];
parameters.co2=[];
parameters.glcin=[];
parameters.growth=[];

i=1 ;

while true
    %% iteration over growth without bool activity
    if x_val(i) < 0.26 %3.3
        mutant = setParam(mutant, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn', 'GROWTH'}, [1.3, 1.7, 0, 0, 0, 0, x_val(i)]);
        mutant = setParam(mutant, 'obj', {'GROWTH'}, 1);
        sol = solveLP(mutant,1);
        fluxmodel = mutant;
    % turn on bool activity (modelled as (0|1) -> (1|1))
    elseif x_val(i) == 0.26 %3.3
        %Change glucose to 1
        settings.gluc(1)=1;
        [Bool, TFAct, EnzAct] = runBool(Bool,settings);
        [GeneExpr] = expression(TFAct.End, model); % values are between -4 and 4 if deviding into differnt levels later.
        [GeneExpr] = postTranscriptionalModifications(EnzAct.End, GeneExpr);% values are between -4 and 5 if deviding into differnt levels later.This is when PTM gives 2 in rank.
        [mutant] = constrainFBA(GeneExpr, model, settings.perc);
    elseif x_val(i) >= 0.26 %3.5
        mutant = setParam(mutant, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn', 'GROWTH'}, [1.3, 1.7, 0, 0, 0, 0, x_val(i)]);
        mutant = setParam(mutant, 'obj', {'GROWTH'}, 1);
        sol = solveLP(mutant,1);
        fluxmodel=mutant;      
    end
    parameters.eth(i) = sol.x(5);
    parameters.o2(i) = sol.x(8);
    parameters.ac(i)=sol.x(1);
    parameters.co2(i)=sol.x(4);
    parameters.growth(i) = sol.x(40);
    parameters.glcin(i) = sol.x(7);

    if i > 2
        if parameters.growth(i) <= parameters.growth(i-2)
            break
        else
            disp(['Iteration ', num2str(i), ' completed!']);
        end
    end
    plotfluxdistr(x_val(i), sol, fluxmodel, path)
    i=i+1; 

end

plotfluxes (parameters.o2, parameters.co2, parameters.eth, parameters.ac, ...
    parameters.glcin, parameters.growth, path)

%% Save
save([path,'settings.mat'],'-struct','settings');
save([path,'parameters.mat'],'-struct','parameters');