function plotfigure4(G0regulatedGenes,G1regulatedGenes)
colorsG0 = zeros(length(G0regulatedGenes.rank),3);
colorsG0(G0regulatedGenes.rank>0,2)= 1;
colorsG0(G0regulatedGenes.rank<0,1)= 1;
colorsG1 = zeros(length(G1regulatedGenes.rank),3);
colorsG1(G1regulatedGenes.rank>0,2)= 1;
colorsG1(G1regulatedGenes.rank<0,1)= 1;

figure(4)
tiledlayout(2,2)
nexttile
for i=1:50;
    scatter([1:length(G0regulatedGenes.pU(:,i))],(G0regulatedGenes.pU(:,i)-G0regulatedGenes.minU(:,i))./G0regulatedGenes.ranges(:,i)*.100,[],colorsG0,'filled')
    hold on
end
hold off

nexttile
for i=1:50;
    scatter([1:length(G1regulatedGenes.pU(:,i))],(G1regulatedGenes.pU(:,i)-G1regulatedGenes.minU(:,i))./G1regulatedGenes.ranges(:,i)*.100,[],colorsG1,'filled')
    hold on
end
hold off

nexttile
for i=60:81;
    scatter([1:length(G0regulatedGenes.pU(:,i))],(G0regulatedGenes.pU(:,i)-G0regulatedGenes.minU(:,i))./G0regulatedGenes.ranges(:,i).*100,[],colorsG0,'filled')
    hold on
end
hold off

nexttile
for i=60:81;
    scatter([1:length(G1regulatedGenes.pU(:,i))],(G1regulatedGenes.pU(:,i)-G1regulatedGenes.minU(:,i))./G1regulatedGenes.ranges(:,i).*100,[],colorsG1,'filled')
    hold on
end
hold off