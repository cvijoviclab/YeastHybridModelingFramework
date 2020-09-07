function[PKApwAct, Snf1pwAct, TORpwAct, EnzymesAct, TargetsAct, MetabolitesAct] = ...
    activityConverter(PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, MetabolitesOld)

%written by: Julia M�nch
%date: 2019-11-12
%description: this function converts each pathway component's vector
%   describing (1) presence, (2) phosphorylation and (3) specific activation
%   into a binary value indicating the component's activity (0 = inactive,
%   1 = active) to enable better depiction via 'plotTransitions'
%arguments: 
%   1.-4. input values for each iteration in 'reachSteadyState' (PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld)
%   5. number of iteration
%   6.-10. matrix containing converted activities for each time step that is
%   extended after each iteration (transPKApwAct, transSnf1pwAct, transTORpwAct, transEnzymesAct, transTargetsAct)
%   11. path of folder containing data for specific nutrient conditions
%returns: matrix for each pathway depicting the states of the components'
%   activity transitions; saves matrices as txt files in folder specified
%   under path

 
%% PKA %%
activity = zeros(height(PKApwOld), 1);
if PKApwOld{1,4} == 1 && PKApwOld{1,2} == 1
    activity(1) = 1;
end
if PKApwOld{2,4} == 1 && PKApwOld{2,2} == 1
    activity(2) = 1;
end
if PKApwOld{3,4} == 1 && PKApwOld{3,2} == 1
    activity(3) = 1;
end
if PKApwOld{4,4} == 1 && PKApwOld{4,2} == 1
    activity(4) = 1;
end
if PKApwOld{5,4} == 1 && PKApwOld{5,2} == 1
    activity(5) = 1;
end
if PKApwOld{6,2} == 1 && PKApwOld{6,4} == 1
    activity(6) = 1;
end
if PKApwOld{7,4} == 1 && PKApwOld{7,2} == 1
    activity(7) = 1;
end
if PKApwOld{8,3} == 1 && PKApwOld{8,2} == 1
    activity(8) = 1;
end
if PKApwOld{9,2} == 1
    activity(9) = 1;
end
if PKApwOld{10,2} == 1
    activity(10) = 1;
end
if PKApwOld{11,2} == 1
    activity(11) = 1;
end
if PKApwOld{12,2} == 1
    activity(12) = 1;
end
if PKApwOld{13,2} == 1 && PKApwOld{13,4} == 1
    activity(13) = 1;
end
if PKApwOld{14,3} == 0 && PKApwOld{14,2} == 1
    activity(14) = 1;
end
if PKApwOld{15,3} == 1 && PKApwOld{15,2} == 1
    activity(15) = 1;
end
if PKApwOld{16,3} == 1 && PKApwOld{15,2} == 1
    activity(16) = 1;
end

PKApwAct = [table(PKApwOld{:,1}),table(activity)];

%% Snf1pw %%
activity = zeros(height(Snf1pwOld), 1);

if Snf1pwOld{1,2} == 1
    activity(1) = 1;
end
if Snf1pwOld{2,2} == 1
    activity(2) = 1;
end
if Snf1pwOld{3,2} == 1
    activity(3) = 1;
end
if Snf1pwOld{4,2} == 1 && Snf1pwOld{4,4} == 1
    activity(4) = 1;
end
if Snf1pwOld{5,2} == 1 && Snf1pwOld{5,3} == 1
    activity(5) = 1;
end
if Snf1pwOld{6,2} == 1 
    activity(6) = 1;
end
if Snf1pwOld{7,2} == 1
    activity(7) = 1;
end
if Snf1pwOld{8,2} == 1
    activity(8) = 1;
end
if Snf1pwOld{9,2} == 1 && Snf1pwOld{9,3} == 1
    activity(9) = 1;
end
if Snf1pwOld{10,2} == 1
    activity(10) = 1;
end
if Snf1pwOld{11,2} == 1 && Snf1pwOld{11,3} == 1
    activity(11) = 1;
end
if Snf1pwOld{12,2} == 1 && Snf1pwOld{12,3} == 1
    activity(12) = 1;
end
if Snf1pwOld{13,2} == 1 && Snf1pwOld{13,3} == 1 && PKApwOld{13,4} == 0
    activity(13) = 1;
end
if Snf1pwOld{14,2} == 1 && Snf1pwOld{14,3} == 0
    activity(14) = 1;
end

Snf1pwAct = [table(Snf1pwOld{:,1}),table(activity)];

%% TORpw %%
activity = zeros(height(TORpwOld), 1);

if TORpwOld{1,2} == 1 && TORpwOld{1,4} == 1
    activity(1) = 1;
end
if TORpwOld{2,2} == 1
    activity(2) = 1;
end
if TORpwOld{3,2} == 1
    activity(3) = 1;
end
if TORpwOld{4,2} == 1
    activity(4) = 1;
end
if TORpwOld{5,2} == 1
    activity(5) = 1;
end
if TORpwOld{6,2} == 1
    activity(6) = 1;
end
if TORpwOld{7,2} == 1 && TORpwOld{7,4} == 1
    activity(7) = 1;
end
if TORpwOld{8,2} == 1 && TORpwOld{8,3} == 1
    activity(8) = 1;
end
if TORpwOld{9,2} == 1 && TORpwOld{9,3} == 1
    activity(9) = 1;
end
if TORpwOld{10,2} == 1 && TORpwOld{10,3} == 1
    activity(10) = 1;
end
if TORpwOld{11,2} == 1 && TORpwOld{11,4} == 1
    activity(11) = 1;
end
% Mks: not sure if active when phosphorylated (binds Bmh1/2) or
% unphosphorylated (binds Rgt2)
if TORpwOld{12,2} == 1 && TORpwOld{12,3} == 1
    activity(12) = 1;
end
if TORpwOld{13,2} == 1
    activity(13) = 1;
end
if TORpwOld{14,2} == 1 && TORpwOld{14,3} == 0
    activity(14) = 1;
end
if TORpwOld{15,2} == 1 && TORpwOld{15,3} == 0
    activity(15) = 1;
end
if TORpwOld{16,2} == 1 && TORpwOld{16,3} == 0
    activity(16) = 1;
end

TORpwAct = [table(TORpwOld{:,1}),table(activity)];

%% Enzymes %%
activity = zeros(height(EnzymesOld), 1);

if EnzymesOld{1,2} == 1 && EnzymesOld{1,3} == 0 %if pfk is phosphorylated -> acts as F-2,6-BP
    activity(1) = 1;
end
if EnzymesOld{2,2} == 1 && EnzymesOld{2,3} == 1
    activity(2) =1;
end
if EnzymesOld{3,2} == 1 && EnzymesOld{3,3} == 1 %(Portela, 2006), increases activity in the absence of fructose-2,6-biphosphate
    activity(3) = 1;
end
if EnzymesOld{4,2} == 1 && EnzymesOld{4,3} == 0 %(Horn, Holzer, 1986)
    activity(4) = 1;
end
if EnzymesOld{5,2} == 1 && EnzymesOld{5,3} == 0
    activity(5) = 1;
end
if EnzymesOld{6,2} == 1 && EnzymesOld{6,3} == 1
    activity(6) = 1;
end
if EnzymesOld{7,2} == 1 && EnzymesOld{7,3} == 1
    activity(7) = 1;
end

EnzymesAct = [table(EnzymesOld{:,1}),table(activity)];

%% Targets %% 
%no conversion necessary

TargetsAct = TargetsOld;

%% Metabolites%%

MetabolitesAct = MetabolitesOld;
MetabolitesOld.Properties.VariableNames{2} = 'activity';

end