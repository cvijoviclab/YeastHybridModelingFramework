function [BoolOut] = knockout(BoolIn, knockouts)

%taken from group folder in dropbox, slightly modified
%date: 2019-10-30
%description: this function modifies each pathway's components by setting
%   the genes/proteins in knockouts to 0
%arguments:
%   1.-6. pathway tables to be changed (Enzymes, PKApw, Snf1pw, TORpw, Enzymes, Targets)
%   7. knockouts = cell of strings with components that should be knocked
%   out
%returns: modified tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create knockouts %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Metabolites = BoolIn.Metabolites;
PKApw = BoolIn.PKApw;
Snf1pw = BoolIn.Snf1pw;
TORpw = BoolIn.TORpw;
Enzymes = BoolIn.Enzymes;
Targets = BoolIn.Targets;

nKnockouts = length(knockouts);
 
for i = 1:nKnockouts
    proteinName = knockouts{i};
     
    % find position of protein in data files
    for j = {Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets}
        players = j{1}.Name;
        pos = find(proteinName == players, 1);
        if ~isempty(pos)
            break;
        end
    end
     
    % change presence/activity in corresponding table to 0
    if isequal(j{1},Metabolites)
        Metabolites{pos,2} = 0;
    elseif isequal(j{1},PKApw)
        PKApw{pos,2} = 0;
    elseif isequal(j{1},Snf1pw)
        Snf1pw{pos,2} = 0;
    elseif isequal(j{1},TORpw)
        TORpw{pos,2} = 0; 
    elseif isequal(j{1},Enzymes)
        Enzymes{pos,2} = 0; 
    elseif isequal(j{1},Targets)
        Targets{pos,2} = 0;
    end
end
BoolOut.Metabolites = Metabolites;
BoolOut.PKApw = PKApw;
BoolOut.Snf1pw = Snf1pw;
BoolOut.TORpw = TORpw;
BoolOut.Enzymes = Enzymes;
BoolOut.Targets = Targets;
end