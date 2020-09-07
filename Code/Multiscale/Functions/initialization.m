function [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = initialization ()
%written by: Julia Münch
%date: 2019-10-31
%description: function to define the pathways and respective components 
%and set the initial condition for a Boolean Network including
%   (1) Metabolites
%   (2) Target Genes
%   (3) Enzymes that are regulated and serve as input for stochastic model
%   (4) PKA pathway
%   (5) Snf1 pathway
%   (6) TOR pathway
%returns: matrix representation of each pathway 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% components and initial values %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Metabolites%%
Name = [string('GLUex'), string('ATP'), string('cAMP'), string('F16BP'), string('NH3')]';
presence = ones(length(Name), 1);
Metabolites = table(Name, presence);
Metabolites{3,2} = 0;
Metabolites{4,2} = 0;
Metabolites{5,2} = 0;

%%% Target genes %%
Name = [string('PDS'),string('STRE'), string('RTG'), string('NCR'), string('SUC2'), string('ADH2'), string('CSRE'), string('PCK1'), string('FBP1'), string('Ribosomal')]';
activity = zeros(length(Name), 1);
Targets = table(Name, activity);

%%% Connection to metabolic model/regualted enzymes%%
Name = [string('PFK2'), string('TREH'), string('PK'), string('FBP1'), string('ACC'), string('HXK2'), string('FBP2')]';
presence = ones(length(Name), 1);
phosphorylation = zeros(length(Name),1);
Enzymes = table(Name, presence, phosphorylation);
Enzymes{6,3} = 1;

%%% PKApw %%%
Name = [string('Gpr1'), string('Gpa2'), string('Krh'), string('Cdc25'), string('Ras'), string('Ira'), string('AC'), string('Pde'), string('Tpk1'), string('Tpk2'), string('Tpk3'), string('Bcy1'), string('PKA'), string('Rim15'), string('Gis1'), string('Msn2,4')]';
presence = ones(length(Name), 1);
phosphorylation = zeros(length(Name),1);
%specific activation for RAS (GDP (0) or GTP(1) bound), GPA2(GDP (0) or GTP(1) bound) BCY (no cAMP(0) or
%cAMP bound(1), CDC25 (no (0) or direct (1) interaction with F-1,6-BP)
spec_activation = zeros(length(Name),1);
PKApw = table(Name, presence, phosphorylation, spec_activation);
PKApw{5,4} = 0; %RAS bound to GDP
PKApw{4,4} = 0; %CDC25 does not interact with F-1,6-BP
PKApw{12,4} = 0; %BCY not bound to cAMP
PKApw{13,2} = 0; %PKA-complex -> not abundant
 
%%% Snf1pw %%%
Name = [string('Tos3'), string('Sak1'), string('Elm1'), string('Glc7'), string('Reg1'), string('Sip2'), string('Sip1'), string('Gal83'), string('Snf1'), string('Snf4'), string('Cat8'), string('Sip4'), string('Adr1'), string('Mig1')]';
presence = ones(length(Name), 1);
phosphorylation = zeros(length(Name),1);
spec_activation = zeros(length(Name),1);
Snf1pw = table(Name, presence, phosphorylation, spec_activation);
Snf1pw{5,3} = 1; %Reg1 is phosphorylated
Snf1pw{9,3} = 1; %Snf1 is phosphorylated

%%% TORpw %%% mainly regulated by nitrogen availability
Name = [string('EGO'), string('Tor1'), string('Tor2'), string('Kog1'), string('Tco89'), string('Lst8'), string('TORC1')...
    string('Sch9'), string('Sfp1'), string('Tap42'), string('PP2A'), string('Mks1'), string('Rtg2'), string('Rtg1,3')...
    string('Gln3'), string('Gat1')]';
presence = ones(length(Name), 1);
phosphorylation = zeros(length(Name), 1);
spec_activation = zeros(length(Name),1);
TORpw = table(Name, presence, phosphorylation,spec_activation);
TORpw{7,2} = 0; %TORC not present
end