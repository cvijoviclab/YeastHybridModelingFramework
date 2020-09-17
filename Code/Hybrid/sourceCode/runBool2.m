function [BoolOut, TFAct, EnzAct] = runBool2(BoolIn,settings,path)
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
        %% convert modifications into activity
    if iteration == 1
        transPKApwAct = [];
        transSnf1pwAct = [];
        transTORpwAct = [];
        transEnzymesAct = [];
        transTargetsAct = [];
    end
    
    [transPKApwAct, transSnf1pwAct, transTORpwAct, transEnzymesAct, transTargetsAct] = ...
        activityConverter(PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, iteration,...
        transPKApwAct, transSnf1pwAct, transTORpwAct, transEnzymesAct, transTargetsAct);
      
    %% save initial conditions and...
    if (iteration == 1)
        writetable(MetabolitesOld, [path, 'Metabolites_init.txt'], 'Delimiter', '\t');
        writetable(PKApwOld, [path, 'PKApw_init.txt'], 'Delimiter', '\t');
        writetable(Snf1pwOld, [path, 'Snf1pw_init.txt'], 'Delimiter', '\t');
        writetable(TORpwOld, [path, 'TORpw_init.txt'], 'Delimiter', '\t');
        writetable(EnzymesOld, [path, 'Enzymes_init.txt'], 'Delimiter', '\t');
        writetable(TargetsOld, [path, 'Targets_init.txt'], 'Delimiter', '\t');
        
    % ...create file that contains transitions (is extended after each iteration)
        transMetabolites = MetabolitesOld;
        transPKApw = PKApwOld;
        transSnf1pw = Snf1pwOld;
        transTORpw = TORpwOld;
        transEnzymes = EnzymesOld;
        transTargets = TargetsOld;
        transitions = {transMetabolites, transPKApw, transSnf1pw, transTORpw, transEnzymes, transTargets};
        transNames = {'transMetabolites', 'transPKApw', 'transSnf1pw', 'transTORpw', 'transEnzymes', 'transTargets'};
    end
    
    [transitions] = saveTransitions(iteration, MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld, transitions);
   
   
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
%writetable(Metabolites, [path, 'Metabolites.txt'], 'Delimiter', '\t');
%writetable(PKApw, [path, 'PKApw.txt'], 'Delimiter', '\t');
%writetable(Snf1pw, [path, 'Snf1pw.txt'], 'Delimiter', '\t');
%writetable(TORpw, [path, 'TORpw.txt'], 'Delimiter', '\t');
%writetable(Enzymes, [path, 'Enzymes.txt'], 'Delimiter', '\t');
%writetable(Targets, [path, 'Targets.txt'], 'Delimiter', '\t');

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
 writetable(MetabolitesAct, [path, 'Metabolites.txt'], 'Delimiter', '\t');
 writetable(PKApwAct, [path,'PKApw.txt'], 'Delimiter', '\t');
 writetable(Snf1pwAct, [path,'Snf1pw.txt'], 'Delimiter', '\t');
 writetable(TORpwAct, [path,'TORpw.txt'], 'Delimiter', '\t');
 writetable(EnzymesAct, [path,'Enzymes.txt'], 'Delimiter', '\t');
 writetable(TargetsAct, [path,'Targets.txt'], 'Delimiter', '\t');
 writetable(TFAct.End, [path,'TF.txt'], 'Delimiter', '\t');
 
 
%save files that contain transitions
for i = 1:length(transitions)
    writetable(transitions{i}, [path, transNames{i}, '.txt'], 'Delimiter', '\t');
end

%save files that contain activity transitions
writetable(transPKApwAct, [path, 'transPKApw.txt'], 'Delimiter', '\t');
writetable(transSnf1pwAct, [path, 'transSnf1pw.txt'], 'Delimiter', '\t');
writetable(transTORpwAct, [path, 'transTORpw.txt'], 'Delimiter', '\t');
writetable(transEnzymesAct, [path, 'transEnzymes.txt'], 'Delimiter', '\t');
writetable(transTargetsAct, [path, 'transTargets.txt'], 'Delimiter', '\t');


end

    
    




