% make table for kcats
load('../../models/reduced_ecYeast_fermentation.mat');
model = ecModel_ferm;
enzGenes = model.enzGenes;
kcat={}

GeneID = enzGenes;
for i=1:length(GeneID)
    %Setting enzyme to find
    if ~isempty(i)
        enzyme  =model.enzymes(i);
    else
        enzyme = [];
    end
    
    if ~isempty(enzyme)
    enzName    = ['prot_' enzyme{1}];
    enzMetIndx = find(strcmpi(model.metNames,enzName));
    enzKcats   = find(model.S(enzMetIndx,:));
    enzKcats   = enzKcats(1:end-1);
    kcat{i,1} = [model.S(enzMetIndx,enzKcats)];
    end
end
kcat=table(enzGenes,kcat)
%kcat.kcat{41}(1)=kcat.kcat{41}(1)*0.5