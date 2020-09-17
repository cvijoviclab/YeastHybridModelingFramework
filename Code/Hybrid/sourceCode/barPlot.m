function barPlot(pathway, path)

%written by: Julia Münch
%date: 2019-11-07
%description: this function depicts and saves the pathway components'
%activity changes when changing the nutrient conditions
%arguments:
%   1. pathway to depict (Enzymes, PKA, Snf1, Targets, TOR)
%   2. path to save images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process simulated data %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import = readtable([path, 'Activity/Transitions/', pathway], 'ReadVariableNames', true, 'ReadRowNames', true);
Name = categorical([import.Properties.RowNames]);
start= import{:,1};
stop = import{:, size(import,2)};
data = stop - start;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plot data %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

b = bar(Name, data);
b.FaceColor = [0.6,0.6,0.6];
ylabel('Change in Activity');
if isequal(pathway, 'transEnzymes.txt') || isequal(pathway, 'transTargets.txt')
    title(strrep(strrep(pathway, 'trans', ''), '.txt', ''));
elseif isequal(pathway, 'transTF.txt')
    title('Transcriptional Regulators')
else
    title([strrep(strrep(pathway, 'trans', ''), 'pw.txt', ''), ' Pathway']);
end    
b = gcf;
%b.CurrentAxes.DataAspectRatio = [1.5 1 1];
b.CurrentAxes.PlotBoxAspectRatio = [2.8 1 1];
b.CurrentAxes.FontSize = 16;
b.CurrentAxes.YTick = [-1 1];
b.CurrentAxes.YTickLabel = {'off', 'on'};
b.CurrentAxes.TickLength = [0 0];
b.CurrentAxes.YLim = [-1.1 1.1];
b.CurrentAxes.Box = 'off';
b.CurrentAxes.XColor = [0 0 0];
b.CurrentAxes.YColor = [0 0 0];

print([path,'Figures/Change/' strrep(strrep(pathway, 'trans', ''), '.txt', '')],'-dpng','-r0')
close(b)



