function plotTransitions(pathway, path)

%written by: Julia M�nch
%date: 2019-11-07
%description: this function depicts and saves the states' transitions when changing
%the nutrient conditions
%arguments:
%   1. pathway to depict (Enzymes, PKA, Snf1, Targets, TOR)
%   2. path to save images

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process simulated data %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
import = readtable([path, 'Activity/Transitions/', pathway], 'ReadVariableNames', true, 'ReadRowNames', true);
Names = [import.Properties.RowNames];
data = table2array(import);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Create Heatmap %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h = heatmap(data);

%create title depending on input
if isequal(pathway, 'transTargets.txt')
    h.Title = 'Target Gene Expression';
elseif isequal(pathway, 'transEnzymes.txt')
    h.Title = 'Enzyme Activity';
elseif isequal(pathway, 'transTF.txt')
    h.Title = 'Activity of Transcriptional Regulators';    
else
    h.Title = ['Activity of ', strrep(strrep(pathway, 'trans', ''), 'pw.txt', ''), ' Pathway Components'];
end

%remove y-label
h.YLabel = '';

%present gene names in same order as input (don't sort alphabetically)
h.YData = Names;

%determine x-label
h.XLabel = 'Iteration';

%change x-ticks
h.XData = num2cell(0:(size(data,2)-1));

%determine font size
h.FontSize = 18;

%select colors (white & blue)
h.Colormap = [[1 1 1]; [122/255 146/255 164/255]];

%remove colorbar
h.ColorbarVisible = 'off';

%remove cell label
h.CellLabelColor = 'none';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save images %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=gcf;
print([path, 'Figures/Transitions/', strrep(pathway, '.txt', '')],'-dpng','-r0')
close(fig)

end