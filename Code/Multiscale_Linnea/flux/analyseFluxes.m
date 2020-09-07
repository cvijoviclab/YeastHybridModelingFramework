%analyse fluxes

Regulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/flux/Simulations/multiscale20200601/fluxDist_reducedYeast.txt');
Unregulated=readtable('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/flux/Simulations/multiscaleUnregulated20200601/fluxDist_reducedYeast.txt');

Difference=Unregulated{:,6:86}-Regulated{:,6:86};
figure()
plot(Difference)

bar(max(Difference,[],2))
min(Difference,[],2)

differneceOA=table();
differneceOA=Regulated(:,[1 2 5]);
differneceOA.sumofdiff=sum(Difference,2);

differnece01=table();
differnece01=Regulated(:,[1 2 5]);
differnece01.sumofdiff=Difference(:,21);

differnece029=table();
differnece029=Regulated(:,[1 2 5]);
differnece029.sumofdiff=Difference(:,59);

differnece04=table();
differnece04=Regulated(:,[1 2 5]);
differnece04.sumofdiff=Difference(:,81);