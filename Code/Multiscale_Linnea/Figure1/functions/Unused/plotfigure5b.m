function plotfigure5b(modelProteomics,modelSimulated)
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
scatter(modelSimulated.TotalAve,modelProteomics.TotalAve,[],comb,'filled')
ylabel('Proteomics Data','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off
nexttile
scatter(modelSimulated.respitrationAve(modelProteomics.G0regulation~=0,:),modelProteomics.h33(modelProteomics.G0regulation~=0,:),[],colorsG0(modelProteomics.G0regulation~=0,:),'filled')
xlabel('Simulated Data','FontSize',18)
ylabel('Proteomics Data','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off
nexttile
scatter(modelSimulated.fermentationAve(modelProteomics.G1regulation~=0,:),modelProteomics.h5(modelProteomics.G1regulation~=0,:),[],colorsG1(modelProteomics.G1regulation~=0,:),'filled')
xlabel('Simulated Data','FontSize',18)
ylabel('Proteomics Data','FontSize',18)
hold on
plot(linspace(0,0.03,50),linspace(0,0.03,50),'LineWidth',3)
hold off

% ylim([0 0.07])
% scatter(modelSimulated.respitrationAve(modelProteomics.G0regulation~=0,:),modelProteomics.respirationAve(modelProteomics.G0regulation~=0,:))
% scatter(modelSimulated.fermentationAve(modelProteomics.G1regulation~=0,:),modelProteomics.fermentationAve(modelProteomics.G1regulation~=0,:))