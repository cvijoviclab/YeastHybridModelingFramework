function[rank] = getRanks(TFAct)

%written by: Julia M�nch
%date: 02-12-07
%description: this function returns simulated ranks for gene expression by
%adding up all positive and negative regulations (derived from Yeastract database) that the specific gene
%receives by the transcription factors included in the boolean model
%arguments: 
%   1. matrix of transitions of transcription factor activity
%returns: vector containing matrices of each pathway describing components' states after each iteration   

for i=2:width(TFAct)
    TFAct.(i) = str2double(TFAct{:,i});
end
%enzNames = model.enzNames;
enzNames = table2array(readtable('RequiredDocs/POI.txt', 'ReadVariableNames', false));
rank = zeros(length(enzNames),1);
rank = table(enzNames, rank);
var = size(TFAct, 2);

expression = '\,';
TF = readtable('RequiredDocs/TF-Targets-import.xlsx');
GisUp = regexp(cell2mat(TF{1,2}), expression, 'split');
GisDown = regexp(cell2mat(TF{1,3}), expression, 'split');
Msn2Up = regexp(cell2mat(TF{2,2}), expression, 'split');
Msn2Down = regexp(cell2mat(TF{2,3}), expression, 'split');
Msn4Up = regexp(cell2mat(TF{3,2}), expression, 'split');
Msn4Down = regexp(cell2mat(TF{3,3}), expression, 'split');
AdrUp = regexp(cell2mat(TF{4,2}), expression, 'split');
AdrDown = regexp(cell2mat(TF{4,3}), expression, 'split');
MigUp = regexp(cell2mat(TF{5,2}), expression, 'split');
MigDown = regexp(cell2mat(TF{5,3}), expression, 'split');
SipUp = regexp(cell2mat(TF{6,2}), expression, 'split');
SipDown = regexp(cell2mat(TF{6,3}), expression, 'split');
CatUp = regexp(cell2mat(TF{7,2}), expression, 'split');
CatDown = regexp(cell2mat(TF{7,3}), expression, 'split');
GlnUp = regexp(cell2mat(TF{8,2}), expression, 'split');
GlnDown = regexp(cell2mat(TF{8,3}), expression, 'split');
GatUp = regexp(cell2mat(TF{9,2}), expression, 'split');
GatDown = regexp(cell2mat(TF{9,3}), expression, 'split');
Rtg1Up = regexp(cell2mat(TF{10,2}), expression, 'split');
Rtg1Down = regexp(cell2mat(TF{10,3}), expression, 'split');
Rtg3Up = regexp(cell2mat(TF{11,2}), expression, 'split');
Rtg3Down = regexp(cell2mat(TF{11,3}), expression, 'split');
Sfp1Up = regexp(cell2mat(TF{12,2}), expression, 'split');
Sfp1Down = regexp(cell2mat(TF{12,3}), expression, 'split');

%% part 1: sums up TF influence on each enzyme -> rank for each enzyme
if TFAct{1,var} == 1
    up = find(ismember(enzNames, GisUp));
    down = find(ismember(enzNames, GisDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{2,var} == 1
    up = find(ismember(enzNames, Msn2Up));
    down = find(ismember(enzNames, Msn2Down));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{3,var} == 1
    up = find(ismember(enzNames, Msn4Up));
    down = find(ismember(enzNames, Msn4Down));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{4,var} == 1
    up = find(ismember(enzNames, AdrUp));
    down = find(ismember(enzNames, AdrDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2}+ 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{5,var} == 1
    up = find(ismember(enzNames, MigUp));
    down = find(ismember(enzNames, MigDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2}+ 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{6,var} == 1
    up = find(ismember(enzNames, SipUp));
    down = find(ismember(enzNames, SipDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{7,var} == 1
    up = find(ismember(enzNames, CatUp));
    down = find(ismember(enzNames, CatDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{8,var} == 1
    up = find(ismember(enzNames, GlnUp));
    down = find(ismember(enzNames, GlnDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{9,var} == 1
    up = find(ismember(enzNames, GatUp));
    down = find(ismember(enzNames, GatDown));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{10,var} == 1
    up = find(ismember(enzNames, Rtg1Up));
    down = find(ismember(enzNames, Rtg1Down));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{11,var} == 1
    up = find(ismember(enzNames, Rtg3Up));
    down = find(ismember(enzNames, Rtg3Down));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end
if TFAct{12,var} == 1
    up = find(ismember(enzNames, Sfp1Up));
    down = find(ismember(enzNames, Sfp1Down));
    for i = 1:length(up)
        rank{up(i),2} = rank{up(i),2} + 1;
    end
    for i = 1:length(down)
        rank{down(i),2} = rank{down(i),2} - 1;
    end
end