function plotfigure3(G0regulatedGenes,G1regulatedGenes)
colorsG0 = zeros(length(G0regulatedGenes.rank),3);
colorsG0(G0regulatedGenes.rank>0,2)= 1;
colorsG0(G0regulatedGenes.rank<0,1)= 1;
colorsG1 = zeros(length(G1regulatedGenes.rank),3);
colorsG1(G1regulatedGenes.rank>0,2)= 1;
colorsG1(G1regulatedGenes.rank<0,1)= 1;
figure(3)
tiledlayout(2,2)
nexttile
scatter([1:length(G0regulatedGenes.pU(:,20))],(G0regulatedGenes.pU(:,20)-G0regulatedGenes.minU(:,20))./G0regulatedGenes.ranges(:,20)*100,[],colorsG0,'filled')
nexttile
scatter([1:length(G1regulatedGenes.pU(:,20))],(G1regulatedGenes.pU(:,20)-G1regulatedGenes.minU(:,20))./G1regulatedGenes.ranges(:,20)*100,[],colorsG1,'filled')
nexttile
scatter([1:length(G0regulatedGenes.pU(:,71))],(G0regulatedGenes.pU(:,71)-G0regulatedGenes.minU(:,71))./G0regulatedGenes.ranges(:,71)*100,[],colorsG0,'filled')
nexttile
scatter([1:length(G1regulatedGenes.pU(:,71))],(G1regulatedGenes.pU(:,71)-G1regulatedGenes.minU(:,71))./G1regulatedGenes.ranges(:,71)*100,[],colorsG1,'filled')
