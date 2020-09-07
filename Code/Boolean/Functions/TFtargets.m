function[transTFAct] = TFtargets(transPKApwAct, transSnf1pwAct, transTORpwAct)

%written by: Julia Münch
%date: 02-12-07
%description: this function returns a single matrix containing the activity
%transitions of all transcription factors included in the boolean model
%arguments: 
%   1-3. matrices including the activity transitions of the PKA, Snf1 and
%   Tor pathays
%returns: matrix with transitions of TF activity

TF = [string('Gis'), string('Msn2'), string('Msn4'), string('Adr1'), string('Mig1'), string('Sip4'), string('Cat8'), string('Gln3'), string('Gat1'), string('Rtg1'), string('Rtg3'), string('Sfp1')]';

Gis = transPKApwAct{15,(2:size(transPKApwAct,2))};
Msn2 = transPKApwAct{16,(2:size(transPKApwAct,2))};
Msn4 = transPKApwAct{16,(2:size(transPKApwAct,2))};
Adr1 = transSnf1pwAct{13,(2:size(transSnf1pwAct,2))};
Mig1 = transSnf1pwAct{14,(2:size(transSnf1pwAct,2))};
Sip4 = transSnf1pwAct{12,(2:size(transSnf1pwAct,2))};
Cat8 = transSnf1pwAct{11,(2:size(transSnf1pwAct,2))};
Gln3 = transTORpwAct{15,(2:size(transTORpwAct,2))};
Gat1 = transTORpwAct{16,(2:size(transTORpwAct,2))};
Rtg1 = transTORpwAct{14,(2:size(transTORpwAct,2))};
Rtg3 = transTORpwAct{14,(2:size(transTORpwAct,2))};
Sfp1 = transTORpwAct{9, (2:size(transTORpwAct,2))};

transTFAct = [Gis; Msn2; Msn4; Adr1; Mig1; Sip4; Cat8; Gln3; Gat1; Rtg1; Rtg3; Sfp1];
transTFAct = array2table([TF, transTFAct]);
transTFAct.Properties.VariableNames{1} = 'Name';
for i = 1:(width(transTFAct)-1)
    transTFAct.Properties.VariableNames{i+1} = ['act', num2str(i)];
end
end