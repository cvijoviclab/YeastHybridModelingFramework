function plot_multiscale_results(model,results,data,Drate,signal_idxs,conBool,pathW_protB,path)
%Get relevant rxn and met indexes for simulations and analysis
idxs = get_model_Idxs(model);
exchIndexes = [idxs(5);idxs(6);idxs(7);idxs(8)];
%Relevant pathways for further exploration
CC_paths = {'Glycolysis' 'TCA' 'pentose phosphate' 'Oxidative Phosphorylation' 'Anaerobic excretion'};
%Plot results
figure()
names = {'Glucose' 'Oxygen' 'CO2' 'Ethanol'};
for i=1:(length(exchIndexes))
    plot(results(:,1),results(:,i+1),'LineWidth',3)
    hold on
end
%Add experimental data points
vector = [5 3 4 6];
for i=vector
    scatter(Drate,data{i})
    hold on
end
legend(names)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Exchange fluxes [mmol/gDw h]','FontSize',18)
xlim([0 max(results(:,1))])
ylim([0 25])
hold off
% Make plot for Boolean input
figure()
Reactions  = model.rxnNames(signal_idxs);
plot(0:(0.4/80):0.4,conBool','LineWidth',3)
legend(Reactions)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Fluxes [mmol/gDw h]','FontSize',18)
conBool2 = table(Reactions,conBool);
save([path,'conBool.mat'],'conBool2');
%Plot protein burden by pathways
figure()
for i=1:(length(CC_paths))
    plot(results(:,1),pathW_protB(:,i),'LineWidth',3)
    hold on
end
legend(CC_paths)
xlabel('Dilution rate [1/h]','FontSize',18)
ylabel('Protein burden [mmol/gDw]','FontSize',18)
xlim([0 max(results(:,1))])
ylim([0 0.03])
end