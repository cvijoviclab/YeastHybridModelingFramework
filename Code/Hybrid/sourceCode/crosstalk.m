function [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = ... 
crosstalk(MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld,...
Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, activeCrosstalk)

%written by: Julia Münch
%date: 2019-11-18
%description: Function that defines crosstalks between pathway; depending
%   on state of crosstalks (0 off and 1 on), states of input's components are
%   changed and returned
%arguments: 
%   1.-6. Old pathway matrices (MetabolitesOld, PKApwOld, Snf1pwOld,
%   TORpwOld, EnzymesOld, TargetsOld) to refer to
%   7.-12. Actual pathway matrices to be changed (Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets)
%   13. vector containing which crosstalks should be turned on/off
%%returns: States of each pathways components that were modified according
%%to crosstalk activity

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Define crosstalks %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% #1 %%
    %phosphorylated Sch9 can phosphorylate Rim15
    %(actually not relevant as Sch9 only active if glucose (and nitrogen)
    %there, but then also PKA is active that also phosphorylates Sch9 -> may be
    %relevant when PKA is knocked out)
    if activeCrosstalk(1) == 1 && TORpwOld{8,3} == 1 && PKApwOld{14,2} == 1
        PKApw{14,3} = 1;    
    end

    %% #2 %%
    %active PKA prevents expression of ADH2 through Adr1 inhibition
    %(also redundant, PKA active when glucose abundant whereas Snf1p that activates ADH2
    %expression is active when no glucose abundant --> effect would be only visible if dephosphorylation of Adr1 not able (but here gap, 
    %so else statement and PKA would be only way to turn off ADH2 expression or if PKA constitutively active))
    if activeCrosstalk(2) == 1 && PKApwOld{13,4} == 1 && Snf1pwOld{13,3} == 1 %&& Snf1pwOld{13,2} == 1 && PKApwOld{13,4} == 1
        Targets{6,2} = 0;  
    end

    %% #3 %%
    %Snf1p inhibits activation of TOR
    % -> only active when glucose is abundant 
    %(needed, otherwise model behaves differently from literature)

    if activeCrosstalk(3) == 1 && Snf1pwOld{9,3} == 1 && TORpwOld{1,4} == 1 && TORpwOld{7,2} == 1
        TORpw{7,4} = 0;
    end

    %Snf1p phosphorylates Tap42 (only for model, not proven)
    if activeCrosstalk(3) == 1 && Snf1pwOld{9,3} == 1 && TORpwOld{10,2} == 1 
            TORpw{10,3} = 1;
    end

    if activeCrosstalk(3) == 1 && TORpwOld{10,2} == 1 && TORpwOld{7,4} == 0 && Snf1pwOld{9,3} == 0
            TORpw{10,3} = 0;
    end

    %% Glucose activation of Reg1-Glc7 via PKA

    if activeCrosstalk(4) == 1 && Metabolites{1,2} == 1 && PKApwOld{13,4} == 1 && Snf1pwOld{4,2} == 1 && Snf1pwOld{4,4} == 1 && Snf1pwOld{5,2} == 1 && Snf1pwOld{5,3} == 1 && Snf1pwOld{9,3} == 1 
        Snf1pw{9,3} = 0;
    end

    %% SNF1 deactivates AC/CYR1

    % if activeCrosstalk(5) == 1 && Metabolites{1,2} == 0 && PKApwOld{7,2} == 1 && Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1)
    %         PKApw{7,4} = 0; 
    % end
    if activeCrosstalk(5) == 1 && PKApwOld{2,4} == 0 && PKApwOld{5,4} == 0 && PKApwOld{7,2} == 1 && Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1)
            PKApw{7,4} = 0; 
    end
end