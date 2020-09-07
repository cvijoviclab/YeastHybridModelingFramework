function plotfigure5(modelProteomics)
colorsG0 = zeros(length(modelProteomics.G0regulation),3);
colorsG0(modelProteomics.G0regulation>0,2)= 1;
colorsG0(modelProteomics.G0regulation<0,1)= 1;
colorsG1 = zeros(length(modelProteomics.G1regulation),3);
colorsG1(modelProteomics.G1regulation>0,2)= 1;
colorsG1(modelProteomics.G1regulation<0,1)= 1;
comb=colorsG0+colorsG1;
comb(comb>1)=1/2;
figure()
tiledlayout(3,2)
nexttile([2 2])
scatter(modelProteomics.respirationAve,modelProteomics.fermentationAve,[],comb,'filled')
ylabel('Fermentation','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off
nexttile
scatter(modelProteomics.respirationAve(modelProteomics.G0regulation~=0,:),modelProteomics.fermentationAve(modelProteomics.G0regulation~=0,:),[],colorsG0(modelProteomics.G0regulation~=0,:),'filled')
xlabel('Respiration','FontSize',18)
ylabel('Fermentation','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off
nexttile
scatter(modelProteomics.respirationAve(modelProteomics.G1regulation~=0,:),modelProteomics.fermentationAve(modelProteomics.G1regulation~=0,:),[],colorsG1(modelProteomics.G1regulation~=0,:),'filled')
xlabel('Respiration','FontSize',18)
ylabel('Fermentation','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off