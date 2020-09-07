function [geneExpr]= addSimData(geneExpr,D) 
ranges =[];
minU=[];
maxU=[];
pU=[];
for i = 1:length(D)
    ranges(:,i)=table2array(D{i}(:,2));
    minU(:,i)=table2array(D{i}(:,3));
    maxU(:,i)=table2array(D{i}(:,4));
    pU(:,i)=table2array(D{i}(:,5));
end
geneExpr.ranges = ranges;
geneExpr.minU = minU;
geneExpr.maxU = maxU;
geneExpr.pU = pU;
