function plotActivity(pathway, path, nutr, knockouts)

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
data=[];
for i = 1:length(knockouts)
    import = readtable([path, 'Activity/Transitions/',num2str(i),'/', pathway], 'ReadVariableNames', true, 'ReadRowNames', true);
    data = [data,import{:,length(import{1,:})}];
end
Name = [import.Properties.RowNames];

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

for i = 1:length(knockouts)
    knockouts{i}=join(knockouts{i});
end
h.XLabel = ['Nutrient presence (Glucose|Nitrogen)',nutr];


%determine x-ticks
h.XData = knockouts;

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
print([path,'Figures/Heatmaps/' strrep(strrep(pathway, 'trans', ''), '.txt', '')],'-dpng','-r0')
saveas(fig,[path,'Figures/Heatmaps/' strrep(strrep(pathway, 'trans', ''), '.txt', '.fig')])
close(fig)

end