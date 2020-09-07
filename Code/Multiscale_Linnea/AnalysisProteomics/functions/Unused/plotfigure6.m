function plotfigure6(classifyer)
figure(6)
tiledlayout(1,2)
nexttile
plot(classifyer.G0diff(classifyer.G0class==1,:),'LineWidth',3);
hold on
plot(classifyer.G0diff(classifyer.G0class==2,:),'LineWidth',3);
hold off

nexttile
plot(classifyer.G1diff(classifyer.G1class==1,:),'LineWidth',3);
hold on
plot(classifyer.G1diff(classifyer.G1class==2,:),'LineWidth',3);
hold off