function []=plotfigure1(G0regulatedGenes,G1regulatedGenes)
colorsG0 = zeros(length(G0regulatedGenes.rank),3);
colorsG0(G0regulatedGenes.rank>0,2)= 1;
colorsG0(G0regulatedGenes.rank<0,1)= 1;
colorsG1 = zeros(length(G1regulatedGenes.rank),3);
colorsG1(G1regulatedGenes.rank>0,2)= 1;
colorsG1(G1regulatedGenes.rank<0,1)= 1;
figure(1)
tiledlayout(2,2)
%Ranges including all timepoints
nexttile
boxplot(sort(G0regulatedGenes.ranges(:,:))','Colors',colorsG0);
set(gca,'YScale','log');
%axes('ColorOrder',colors);
nexttile
boxplot(sort(G1regulatedGenes.ranges(:,60:81))','Colors',colorsG1);
set(gca,'YScale','log');
%Ranges including only regulated timepoints
nexttile
boxplot(sort(G0regulatedGenes.ranges(:,1:50))','Colors',colorsG0);
set(gca,'YScale','log');
nexttile
boxplot(sort(G1regulatedGenes.ranges(:,60:81))','Colors',colorsG1);
set(gca,'YScale','log');