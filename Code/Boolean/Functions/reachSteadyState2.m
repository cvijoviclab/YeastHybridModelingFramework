function [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, transTFAct] = ...
    reachSteadyState2(count,MetabolitesIn, PKApwIn, Snf1pwIn, TORpwIn, EnzymesIn, TargetsIn, glucoseLevel, nitrogenLevel, path, activeCrosstalk)

%written by: Julia M�nch
%data: 2019-10-31
%description: function to define the Boolean rules; loops over discrete time
%   steps until states do not change anymore; stop looping if steady state is
%   not reached after 100 interations; then find SS for next nutrient
%   condition; states for inital conditions before changing nutrient conditions, SS
%   for specific nutrient conditions and transitions from init to SS are
%   saved as txt files, activity transitions for each pathway tables are
%   saved as txt files, simulated gene ranks are saved as txt file
%input: 
%   1.-6. matrix representations of pathways at time = t
%   7. glucose level, changed when SS is reached for respective input
%   8. nitrogen level, changed when SS is reached for repective input
%   9. path of folder containing data for specific nutrient conditions 
%returns: matrix representations of pathways at time = t+1
%contains functions:
%   1. activityConverter (converts states in vector form into integer representing activity )
%   2. saveTransitions (saves transitions)
%   3. TFtargets (isolates transitions of TF activity from pahtway tables)
%   4. getRanks (creates ranks for target genes)


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% set starting values %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Metabolites = MetabolitesIn;
PKApw = PKApwIn;
Snf1pw = Snf1pwIn;
TORpw = TORpwIn;
Enzymes = EnzymesIn;
Targets = TargetsIn;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Glucose concentration%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Metabolites{1,2} = glucoseLevel; %determine glucose concentration (0) no, (1) high
%HXT transporter could be added to model?
Metabolites{5,2} = nitrogenLevel;

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
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% PKA pathway %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% GPR1 part %%%

    %Glucose Sensing via GPR1 induces activation of GPA2 (Lorenz et al., 2000)
    if MetabolitesOld{1,2} == 1 && PKApwOld{1,2} == 1
        PKApw{1,4} = 1;
    else
        PKApw{1,4} = 0;
    end

    if PKApwOld{1,4} == 1 && PKApwOld{2,2} == 1
        PKApw{2,4} = 1;
    end
    if (PKApwOld{1,4} == 0 && PKApwOld{2,2} == 1) || PKApwOld{2,2} == 0
        PKApw{2,4} = 0;   
    end

    %active GPA2 inhibits Krh activity (Peeters et al., 2006)
    if (PKApwOld{2,4} == 1 && PKApwOld{3,2} == 1) || PKApwOld{3,2} == 0
        PKApw{3,4} = 0;
    else
        PKApw{3,4} = 1;
    end

    %%% AC part %%%

    %F16BP is present when glucose is present (Peeters et al., 2017)
    if MetabolitesOld{1,2} == 1
        Metabolites{4,2} = 1;
    else
        Metabolites{4,2} = 0;
    end 

    %CDC25 is activated by F16BP (Peeters et al., 2017), exact mechanism
    %unknown
    if MetabolitesOld{4,2} == 1 && PKApwOld{4,2} == 1
        PKApw{4,4} = 1; 
    else
        PKApw{4,4} = 0;
    end

    %IRA is active during glucose depletion (Bond, 2008)
    if MetabolitesOld{1,2} == 0 && PKApwOld{6,2} == 1
        PKApw{6,4} = 1;
    else
        PKApw{6,4} = 0;
    end

    %active CDC25 stimulates guanine nucleotide exchange in RAS from GDP to
    %GTP (Conrad et al., 2014) 
    %include IRA?
    if PKApwOld{5,2} == 1 && PKApwOld{5,4} == 0 && PKApwOld{4,4} == 1 && PKApwOld{6,4} == 0 
        PKApw{5,4} = 1;
    end

    %IRA stimulates guanine nucleotide exchange in RAS from GTP to GDP (Conrad et al., 2014)
    % or in case of RAS knockout -> inactive
    if (PKApwOld{5,2} == 1 && PKApwOld{5,4} == 1 && PKApwOld{4,4} == 0 && PKApwOld{6,4} == 1) || PKApwOld{5,2} == 0
        PKApw{5,4} = 0;
    end

    %AC is activated when RAS is activated, or RAS and GPR2 is active (Conrad et al., 2014) 
    if PKApwOld{5,4} == 1 && PKApwOld{7,2} == 1 || (PKApwOld{5,4} == 1 && PKApwOld{2,4} == 1 && PKApwOld{7,2} == 1)
        PKApw{7,4} = 1;
    end
    if PKApwOld{7,2} == 0
        PKApw{7,4} = 0; 
    end

    %cAMP present when ATP is present and AC is activated ()
    if MetabolitesOld{2,2} == 1 && PKApwOld{7,4} == 1
        Metabolites{3,2} = 1;
    end

    %no cAMP when AC is not activated and PDE is activated
    if PKApwOld{7,4} == 0 && PKApwOld{8,2} == 1 && PKApwOld{8,3} == 1
        Metabolites{3,2} = 0;
    end

    %%% PKA part %%%

    %PKA is present when BCY1 and TPK1/2/3 are present
    if (PKApwOld{9,2} == 1 || PKApwOld{10,2} == 1 || PKApwOld{11,2} == 1)
        PKApw{13,2} = 1;
    else
        PKApw{13,2} = 0;
    end

    %PKA is active when cAMP is there and Krh is inactive 
    if (PKApwOld{13,2} == 1 && MetabolitesOld{3,2} == 1 && PKApwOld{3,4} == 0 && PKApwOld{12,2} == 1) || PKApwOld{12,2} == 0
        PKApw{13,4} = 1;
    end

    %PKA is inactive if KRH is active or no cAMP is there, PKA not there
    if (PKApwOld{13,2} == 1 && PKApwOld{12,2} == 1 && (PKApwOld{3,4} == 1 || MetabolitesOld{3,2} == 0)) || PKApwOld{13,2} == 0
        PKApw{13,4} = 0;
    end

    %PDE acitvated through PKA mediated phosphorylation (Hu et al., 2010)
    if PKApwOld{13,4} == 1 && PKApwOld{8,2} == 1
        PKApw{8,3} = 1;
    else
        PKApw{8,3} = 0;
    end 

    %Downstream of PKA

    %active PKA phosphorylates TREH, PFK2/F-2,6-BP, PK, F-1,6-BP
    if PKApwOld{13,4} == 1 && EnzymesOld{1,2} == 1 
        Enzymes{1,3} = 1;
    else
        Enzymes{1,3} = 0;
    end
    if PKApwOld{13,4} == 1 && EnzymesOld{2,2} == 1 
        Enzymes{2,3} = 1;
    else
        Enzymes{2,3} = 0;
    end
    if PKApwOld{13,4} == 1 && EnzymesOld{3,2} == 1 
        Enzymes{3,3} = 1;
    else
        Enzymes{3,3} = 0;
    end
    if PKApwOld{13,4} == 1 && EnzymesOld{4,2} == 1 
        Enzymes{4,3} = 1;
    else
        Enzymes{4,3} = 0;
    end


    %active PKA as well as phosphorylted Sch9 can phosphorylate Rim15
    %if (PKApwOld{13,4} == 1 || TORpwOld{8,3} == 1) && PKApwOld{14,2} == 1
    if PKApwOld{13,4} == 1 && PKApwOld{14,2} == 1
        PKApw{14,3} = 1; 
    else
        PKApw{14,3} = 0;
    end

    %phosphorylated Rim gets into cytoplasm -> no activation of Gis and Msn activity
    if PKApwOld{14,3} == 0 && PKApwOld{15,2} == 1
        PKApw{15,3} = 1;
    else
        PKApw{15,3} = 0;
    end
    if PKApwOld{14,3} == 0 && PKApwOld{16,2} == 1
        PKApw{16,3} = 1;
    else
       PKApw{16,3} = 0;
    end

    %Gis1 induces transcription of PDS genes (Conrad et al.,2014)
    if PKApwOld{15,3} == 1
        Targets{1,2} = 1;
    else
        Targets{1,2} = 0;
    end

    %Msn induces transcription of STRE genes (Contrad et al., 2014)
    if PKApwOld{16,3} == 1
        Targets{2,2} = 1;
    else
        Targets{2,2} =0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% Snf1 pathway %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% 3 states from Sanz et al., 1999

    %glucose deprivation induces phosphorylation of Hxk2     
    if MetabolitesOld{1,2} == 0 && Enzymes{6,2} == 1
        Enzymes{6,3} = 1;
    else
        Enzymes{6,3} = 0;
    end

    % HIGH GLUCOSE TO LOW GLUCOSE1:
    % Tos3, Sak1 and Elm1 kinases are constitutively active and
    % phosphorylate Snf1 (Busti et al., 2010) 
    if MetabolitesOld{1,2} == 0 && Snf1pwOld{9,2} == 1 && (Snf1pwOld{1,2} == 1 || Snf1pwOld{2,2} == 1 || Snf1pwOld{3,2} == 1) && Snf1pwOld{9,3} == 0 && Snf1pwOld{5,3} == 0
        Snf1pw{9,3} = 1;
    end

    % LOW GLUCOSE1 TO LOW GLUCOSE2:
    % phosphorylation of Snf1 (Snf1 complex there) leads to phosphorylation of Reg1,
    % phosphorylated Hxk2 promotes Reg1 phosphorylation
    if Snf1pwOld{4,2} == 1 && Snf1pwOld{5,2} == 1 && Snf1pwOld{5,3} == 0 && Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) || (EnzymesOld{6,3} == 1)
        Snf1pw{5,3} = 1;
    end

    % LOW GLUCOSE2 TO HIGH GLUCOSE:
    %Reg 1 phosphorylation stimulates Glc7 activity -> dephosphorylation of
    %Reg1 (first after Glc7 has dephosphorylated Snf1, see in crosstalk), dephosphorylation inhibited by Hxk2p and PKA
    if Snf1pwOld{5,2} == 1 && Snf1pwOld{5,3} == 1 && Snf1pwOld{4,2} == 1 
        Snf1pw{4,4} = 1;
    end

    if (Snf1pwOld{5,3} == 0 && Snf1pwOld{4,2} == 1) || Snf1pwOld{4,2} == 0
        Snf1pw{4,4} = 0;
    end

    if EnzymesOld{6,3} == 0 && Snf1pwOld{4,2} == 1 && Snf1pwOld{4,4} == 1 && Snf1pwOld{5,2} == 1 && Snf1pwOld{5,3} == 1 && Snf1pwOld{9,3} == 0
        Snf1pw{5,3} = 0;
    end

    %Glc7 dephoshorylates Snf1 and Hxk2 if glucose is there (include in crosstalk als active PKA required)
    if Snf1pwOld{9,2} == 0
        Snf1pw{9,3} = 0;
    end

    %if Snf1 is phosphorylated and Snf4 and Gal83/Sip1/Sip2 are there, Mig
    %is deactivated and Cat8, Sip4 and Adr1 are activated via
    %phosphorylation (Rahner et al., 1999)
    %ACC is phosphorylated and therefore inactivated
    if Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) && Snf1pwOld{14,2} == 1
        Snf1pw{14,3} = 1;
    else
        Snf1pw{14,3} = 0;
    end

    if Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) && Snf1pwOld{12,2} == 1
        Snf1pw{12,3} = 1;
    else
        Snf1pw{12,3} = 0;
    end

    if Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) && Snf1pwOld{11,2} == 1
        Snf1pw{11,3} = 1;
    else
        Snf1pw{11,3} = 0;
    end

    if Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) && EnzymesOld{5,2} == 1
        Enzymes{5,3} = 1;
    else
        Enzymes{5,3} = 0;
    end

    %Adr phosphorylated by pSnf1 (activity is inhibited by active PKA, now included in crosstalk)
    if Snf1pwOld{9,3} == 1 && Snf1pwOld{10,2} == 1 && (Snf1pwOld{6,2} == 1 || Snf1pwOld{7,2}  == 1 || Snf1pwOld{8,2} == 1) && Snf1pwOld{13,2} == 1 %&& PKApwOld{13,4} == 0
        Snf1pw{13,3} = 1;
    else
        Snf1pw{13,3} = 0;
    end

    %crosstalk with PKA inactivates Snf1pw

    % Mig represses the transcription of SUC2
    if Snf1pwOld{14,3} == 1
        Targets{5,2} = 1;
    else
        Targets{5,2} = 0;
    end

    % Adr1 induces the expression of ADH2
    if Snf1pwOld{13,3} == 1
        Targets{6,2} = 1;
    else
        Targets{6,2} = 0;
    end

    %Sip4 and Cat8 induce expression of CSRE genes (FBP1, PCK1...)
    if Snf1pwOld{12,3} == 1 || Snf1pwOld{11,3} == 1
        Targets{7,2} = 1;
        Targets{8,2} = 1;
        Targets{9,2} = 1;
    else
        Targets{7,2} = 0;
        Targets{8,2} = 0;
        Targets{9,2} = 0;  
    end    

    if TargetsOld{9,2} == 1
        Enzymes{4,2} = 1;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %% TOR pathway %%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

    %%% depending on nitrogen levels = Metabolites{5,2} %%
    %if glucose 1 and nitro 0: TORC inactive -> Sch9 path inacitve and PP2A
    %   branch active
    %if glucose 1 and nitro 1: TORC active -> Sch9 path active and PP2A
    %   branch inactive
    %if glucose 0: Sch9 and PP2A inactive regardless of nitrogen
    %   availability!
    %(Hallett, Luo, Capaldi, 2014)
    %corsstalk: TORC inhibited by pSnf1 

    %EGO is activated when nitrogen/AA is available (Conrad et al., 2014)
    if MetabolitesOld{5,2} == 1 && TORpwOld{1,2} == 1
        TORpw{1,4} = 1;
    else
        TORpw{1,4} = 0;
    end

    %TORC1 consists of Tor1/2, Kog1, Tco89 and Lst8 (Inoue and Nomura, 2017)
    if (TORpwOld{2,2} == 1 || TORpwOld{3,2} == 1) && TORpwOld{4,2} == 1 && TORpwOld{5,2} == 1 && TORpwOld{6,2} == 1
        TORpw{7,2} = 1;
    else
        TORpw{7,2} = 0;
    end

    %active EGO activates TORC (Conrad et al., 2014) (phosphorylated SNF1
    %inhibits activation (Conrad, Junkel et al., 2019), is now included in crosstalks)
    if TORpwOld{1,4} == 1 && TORpwOld{7,2} == 1 %&& Snf1pwOld{9,3} == 0
        TORpw{7,4} = 1;
    else
        TORpw{7,4} = 0;
    end

    %active TOR phosphorylates Sch9, Tap42 and Sfp1 (Conrad et al., 2014; )
    %Sch9
    if TORpwOld{7,4} == 1 && TORpwOld{8,2} == 1
        TORpw{8,3} = 1;
    else
        TORpw{8,3} = 0; %% rapidly dephosphorylated upon starvation (carbon, nitrogen...) (Conrad)
    end

    %Tap42 is phosphorylated if TORC is active or if TORC is inactive and  Snf1 is phosphorylated
    %(glucose is not abundant) --> not known!! but have to find way to
    %phosphorylate Tap42 in glucose and nitrogen depletion
    %include OR (TORpwOld{1,4} == 1 && Snf1pwOld{9,3} == 0 && TORpwOld{7,4} == 0 && TORpwOld{10,2} == 1)
    %if you want to skip step inbetween where TORC is not activated yet and
    %Tap42 is therefore shortly dephosphorylated
    if (TORpwOld{7,4} == 1 && TORpwOld{10,2} == 1) %|| (Snf1pwOld{9,3} == 1 && TORpwOld{10,2} == 1) 
        TORpw{10,3} = 1;
    end

    %Tap42 is dephosphorylated (by PP2A, cannot be included) if TOR is not active 
    %(and Snf1 is unphosphorylated -> now included in crosstalk)
    %(no nitorgen but glucose)
    if TORpwOld{10,2} == 1 && TORpwOld{7,4} == 0 || TORpwOld{10,2} == 0 %&& Snf1pwOld{9,3} == 0 
        TORpw{10,3} = 0;
    end

    %Sfp1 (https://www.yeastgenome.org/locus/SFP1/regulation)
    if TORpwOld{7,4} == 1 && TORpwOld{9,2} == 1
        TORpw{9,3} = 1;
    else
        TORpw{9,3} = 0;
    end

    %active Sfp1p provides negative feedpack on Sch9p -> if feedback
    %activated, crosstalk to Rim15 cannot compensate loss of PKA as Sch9
    %only shortly active
    %     if (TORpwOld{9,3} == 1 && TORpwOld{7,4} == 1 && TORpwOld{8,2} == 1) || TORpwOld{8,2} == 0
    %         TORpw{8,3} = 0;
    %     end

    %active Sfp1p activates expression of genes for ribosomal biogenesis
    if TORpwOld{9,3} == 1
        Targets{10,2} = 1;
    else
        Targets{10,2} = 0;
    end

    %active TOR phosphorylates Mks, also when no glucose (TOR inactive)
    if (TORpwOld{7,4} == 1 && TORpwOld{12,2} == 1) || (Metabolites{1,2} == 0 && TORpwOld{12,2} == 1)
        TORpw{12,3} = 1;
    end

    %active PP2A dephoshporylates Mks
    if (TORpwOld{12,2} == 1 && TORpwOld{11,2} == 1 && TORpwOld{11,4} == 1) || TORpwOld{12,2} == 0
        TORpw{12,3} = 0;
    end

    %unphosphorylated Mks complexes with Rtg2 allowing dephosphorylation of
    %Rtg1,3 (Hill, Van Remmen, 2014) by PP2A --> include PP2A?
    if TORpwOld{12,2} == 1 && TORpwOld{12,3} == 0 && TORpwOld{13,2} == 1 && TORpwOld{14,2} == 1
        TORpw{14,3} = 0;
    else
        TORpw{14,3} = 1;
    end

    %unphosphorylated Rtg1,3 can enter nucleus and induce expression of RTG
    %genes
    if TORpwOld{14,3} == 0 && TORpwOld{14,2} == 1
        Targets{3,2} = 1;
    else
        Targets{3,2} = 0;
    end      

    %%% Sch9 path %%%
    %active (phosphorylated Sch9) inhibits Rim15 via phosphorylation -> included in
    %PKA pathway

    %%% Tap42 pathway %%%
    % phoshporylated Tap42 binds PP2A therefore inhibiting phosphatase
    % activity
    if (TORpwOld{10,3} == 1 && TORpwOld{10,2} == 1 && TORpwOld{11,2} == 1) || TORpwOld{11,2} == 0
        TORpw{11,4} = 0;
    else
        TORpw{11,4} = 1;
    end

    %active PP2A can dephosphorylate Gln3 and Gat1
    if TORpwOld{11,2} == 1 && TORpwOld{11,4} == 1 && TORpwOld{15,2} == 1
        TORpw{15,3} = 0;
    else
        TORpw{15,3} = 1;
    end

    if TORpwOld{11,2} == 1 && TORpwOld{11,4} == 1 && TORpwOld{16,2} == 1
        TORpw{16,3} = 0;
    else
        TORpw{16,3} = 1;
    end

    %dephosphorylated Gln3 and Gat1 can get into nucleus and activate
    %expression of NCR genes
    if (TORpwOld{16,3} == 0 && TORpwOld{16,2} == 1) || (TORpwOld{15,3} == 0 && TORpwOld{15,2} == 1)
        Targets{4,2} = 1;
    else
        Targets{4,2} = 0;
    end
    
    %% include crosstalk
    [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = ...
        crosstalk(MetabolitesOld, PKApwOld, Snf1pwOld, TORpwOld, EnzymesOld, TargetsOld,...
        Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, activeCrosstalk);
    
    %% create table for transcription factors
    
    
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
    else
        disp(['Iteration ', num2str(iteration), ' completed!']);
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

%save results in SS for new nutrient conditions
writetable(Metabolites, [path,num2str(count),'/','Metabolites_SS.txt'], 'Delimiter', '\t');
writetable(PKApw, [path,num2str(count),'/', 'PKApw_SS.txt'], 'Delimiter', '\t');
writetable(Snf1pw, [path,num2str(count),'/', 'Snf1pw_SS.txt'], 'Delimiter', '\t');
writetable(TORpw, [path,num2str(count),'/', 'TORpw_SS.txt'], 'Delimiter', '\t');
writetable(Enzymes, [path,num2str(count),'/', 'Enzymes_SS.txt'], 'Delimiter', '\t');
writetable(Targets, [path,num2str(count),'/', 'Targets_SS.txt'], 'Delimiter', '\t');

%save files that contain transitions
for i = 1:length(transitions)
    writetable(transitions{i}, [path, 'Transitions/',num2str(count),'/', transNames{i}, '.txt'], 'Delimiter', '\t');
end

%save files that contain activity transitions
writetable(transPKApwAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transPKApw.txt'], 'Delimiter', '\t');
writetable(transSnf1pwAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transSnf1pw.txt'], 'Delimiter', '\t');
writetable(transTORpwAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transTORpw.txt'], 'Delimiter', '\t');
writetable(transEnzymesAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transEnzymes.txt'], 'Delimiter', '\t');
writetable(transTargetsAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transTargets.txt'], 'Delimiter', '\t');

%create file with transitions of TF factors
[transTFAct] = TFtargets(transPKApwAct, transSnf1pwAct, transTORpwAct);
writetable(transTFAct, [path, 'Activity/Transitions/',num2str(count),'/', 'transTF.txt'], 'Delimiter', '\t');

%create files with regulated genes
%load('../Combined/reduced_ecYeast_batch.mat');
%model = ecModel_batch;
[ranks] = getRanks(transTFAct);
writetable(ranks, [path, 'enzRanks.txt'], 'Delimiter', '\t');
end
