function [BoolOut, TFAct, EnzAct] = runBool(BoolIn,settings)
%written by: Linnea Österberg and Julia Münch
%date: 2020-04-02
%description: 
%input:
%returns:
%contains functions:1. activityConverter
%                   2. boolRules
%                   3. crosstalk
%                   4. TFtargets

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set starting values %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Metabolites = BoolIn.Metabolites;
PKApw = BoolIn.PKApw;
Snf1pw = BoolIn.Snf1pw;
TORpw = BoolIn.TORpw;
Enzymes = BoolIn.Enzymes;
Targets = BoolIn.Targets;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Glucose concentration%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Metabolites{1,2} = settings.gluc(1); %determine glucose concentration (0) no, (1) high
Metabolites{5,2} = settings.nitr(1);
disp(['Glucose level: ', num2str(Metabolites{1,2}), ', Nitrogen level: ', num2str(Metabolites{5,2})]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% run until steady state is reached %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iteration = 1;

while true
        
    %save old values for comparison
    MetabolitesOld = Metabolites;
    PKApwOld = PKApw;
    Snf1pwOld = Snf1pw;
    TORpwOld = TORpw;
    EnzymesOld = Enzymes;
    TargetsOld = Targets;
    
    if iteration == 1
        [PKApwAct, Snf1pwAct, TORpwAct, EnzymesAct, TargetsAct, MetabolitesAct] = ...
            activityConverter(PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, MetabolitesOld);
        
        %create file with transitions of TF factors
        [TFAct.Start] = TFtargets(PKApwAct, Snf1pwAct, TORpwAct);
        [EnzAct.Start] = EnzymesAct;
    end
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = ...
        boolRules(MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld,...
        Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets);
    %% include crosstalk
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = ...
        crosstalk(MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld,...
        Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, settings.activeCrosstalk);
      
           
    %% if states do not change anymore, steady state is reached -> break the loop
    
    if isequal(MetabolitesOld, Metabolites) && isequal(PKApwOld, PKApw) && ...
            isequal(Snf1pwOld, Snf1pw) && isequal(TORpwOld, TORpw) &&...
            isequal(EnzymesOld, Enzymes) && isequal(TargetsOld, Targets)
        disp('You have reached a steady state! :)');
        break;
%     else
%         disp(['Iteration ', num2str(iteration), ' completed!']);
    end
    
    %% stop simulation if steady state is not reached after 100 iterations
    
    if (iteration == 100)
        disp('Steady state could not be reached. Process aborted.')
        break;
    end
    
    iteration = iteration + 1;
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% readout %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%save states in new SS
% writetable(Metabolites, [path, 'Metabolites.txt'], 'Delimiter', '\t');
% writetable(PKApw, [path, 'PKApw.txt'], 'Delimiter', '\t');
% writetable(Snf1pw, [path, 'Snf1pw.txt'], 'Delimiter', '\t');
% writetable(TORpw, [path, 'TORpw.txt'], 'Delimiter', '\t');
% writetable(Enzymes, [path, 'Enzymes.txt'], 'Delimiter', '\t');
% writetable(Targets, [path, 'Targets.txt'], 'Delimiter', '\t');
BoolOut.Metabolites = Metabolites;
BoolOut.PKApw = PKApw;
BoolOut.Snf1pw = Snf1pw;
BoolOut.TORpw = TORpw;
BoolOut.Enzymes = Enzymes;
BoolOut.Targets = Targets;

% convert modifications into activity
[PKApwAct, Snf1pwAct, TORpwAct, EnzymesAct, TargetsAct, MetabolitesAct] = ...
    activityConverter(PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, MetabolitesOld);

%create file with transitions of TF factors
[TFAct.End] = TFtargets(PKApwAct, Snf1pwAct, TORpwAct);
[EnzAct.End] = EnzymesAct;
      
%save results in SS for new nutrient conditions
% writetable(MetabolitesAct, [path, 'Activity/', 'Metabolites.txt'], 'Delimiter', '\t');
% writetable(PKApwAct, [path, 'Activity/', 'PKApw.txt'], 'Delimiter', '\t');
% writetable(Snf1pwAct, [path, 'Activity/', 'Snf1pw.txt'], 'Delimiter', '\t');
% writetable(TORpwAct, [path, 'Activity/', 'TORpw.txt'], 'Delimiter', '\t');
% writetable(EnzymesAct, [path, 'Activity/', 'Enzymes.txt'], 'Delimiter', '\t');
% writetable(TargetsAct, [path, 'Activity/', 'Targets.txt'], 'Delimiter', '\t');
% writetable(TFActEnd, [path, 'Activity/', 'TF.txt'], 'Delimiter', '\t');
end