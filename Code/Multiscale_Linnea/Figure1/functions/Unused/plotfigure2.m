function plotfigure2(geneExprG0,G0regulatedGenes,G1regulatedGenes)
colorsG0 = zeros(length(G0regulatedGenes.rank),3);
colorsG0(G0regulatedGenes.rank>0,2)= 1;
colorsG0(G0regulatedGenes.rank<0,1)= 1;
colorsG1 = zeros(length(G1regulatedGenes.rank),3);
colorsG1(G1regulatedGenes.rank>0,2)= 1;
colorsG1(G1regulatedGenes.rank<0,1)= 1;
figure(2)
tiledlayout(3,2)
%statistics on pU for all genes all dilution rates
nexttile([1 2])
boxplot(sort(geneExprG0.pU)');
%statistics on pU for regulated genes spanning all dilutionrates
nexttile
boxplot(sort(G0regulatedGenes.pU(:,:))','Colors',colorsG0);
ylim([0 2*10^-4])
nexttile
boxplot(sort(G1regulatedGenes.pU(:,:))','Colors',colorsG1);
ylim([0 2*10^-4])
%statistics on pU for regulated genes spanning G0 or G1
nexttile
boxplot(sort(G0regulatedGenes.pU(:,1:50))','Colors',colorsG0);
ylim([0 2*10^-4])
nexttile
boxplot(sort(G1regulatedGenes.pU(:,60:81))','Colors',colorsG1);
ylim([0 2*10^-4])