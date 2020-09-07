function plotActivity(pathway, path, nutr)

%written by: Julia M�nch
%date: 2019-10-30
%description: this function depicts and saves the activity of pathway
%components in a heatmap comparing the Steady States under different nutrient conditions (glucose|nitrogen)
%arguments:
%   1. pathway to depict (Enzymes, PKA, Snf1, Targets, TOR)
%   2. path to save images
%   3. nutrient conditions to display

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Process simulated data %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

import = readtable([path, 'Activity/Transitions/', pathway], 'ReadVariableNames', true, 'ReadRowNames', true);
Name = [import.Properties.RowNames];
data = [import{:,1}, import{:,size(import,2)}];


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
h.YData = Name;

%determine x-label
h.XLabel = 'Nutrient presence (Glucose|Nitrogen)';

%determine x-ticks
h.XData = nutr;

%determine font size
h.FontSize = 18;

%select colors (red & green)
h.Colormap = [[0.9 0.2 0.1]; [0.1 0.7 0.2]];

%remove colorbar
h.ColorbarVisible = 'off';

%remove cell label
h.CellLabelColor = 'none';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% save images %%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fig=gcf;
print([path,'Figures/Heatmaps/' strrep(strrep(pathway, 'trans', ''), '.txt', '')],'-dpng','-r0')
close(fig)

end