function[TFAct] = TFtargets(PKApwAct, Snf1pwAct, TORpwAct)

%written by: Julia Münch
%date: 02-12-07
%description: this function returns a single matrix containing the activity
%transitions of all transcription factors included in the boolean model
%arguments: 
%   1-3. matrices including the activity transitions of the PKA, Snf1 and
%   Tor pathays
%returns: matrix with transitions of TF activity

TF = [string('Gis'), string('Msn2'), string('Msn4'), string('Adr1'), string('Mig1'), string('Sip4'), string('Cat8'), string('Gln3'), string('Gat1'), string('Rtg1'), string('Rtg3'), string('Sfp1')]';

Gis = PKApwAct{15,2};
Msn2 = PKApwAct{16,2};
Msn4 = PKApwAct{16,2};
Adr1 = Snf1pwAct{13,2};
Mig1 = Snf1pwAct{14,2};
Sip4 = Snf1pwAct{12,2};
Cat8 = Snf1pwAct{11,2};
Gln3 = TORpwAct{15,2};
Gat1 = TORpwAct{16,2};
Rtg1 = TORpwAct{14,2};
Rtg3 = TORpwAct{14,2};
Sfp1 = TORpwAct{9,2};

TFAct = [Gis; Msn2; Msn4; Adr1; Mig1; Sip4; Cat8; Gln3; Gat1; Rtg1; Rtg3; Sfp1];
TFAct = table(TF, TFAct);
%TFAct = array2table([TF, TFAct]);
TFAct.Properties.VariableNames{1} = 'Name';
TFAct.Properties.VariableNames{2} = 'activity';
end