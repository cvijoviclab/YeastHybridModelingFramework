function [proteomics, D, enzUsages_reducedYeast, geneExprG0, geneExprG1]=preProcess(proteomics, D, enzUsages_reducedYeast, geneExprG0, geneExprG1)

%% GeneExpression table
%Add a index variable to table
geneExprG0.indx=[1:127]';
geneExprG1.indx=[1:127]';

% Add pathway information
load('/Users/linoste/Documents/GitHub/Crabtree/Code/Multiscale_Linnea/data/tempModel.mat');
pathways = mapEnzymeSubSystems(tempModel.enzymes,tempModel);
geneExprG0.pathways = pathways;
geneExprG1.pathways = pathways;

% Add proteomics index
% Couples all genes in model to a indecs in the proteomics file.
proteomicsIndex = zeros(size(geneExprG0,1),1);
for i=1:size(geneExprG0,1)
    j = (find(ismember(proteomics.ProteinAccession, geneExprG0.enzGenes(i))));
    if ~isempty(j)
        proteomicsIndex(i) = j;
    end
end
geneExprG0.proteomicsIndex = proteomicsIndex;
geneExprG1.proteomicsIndex = proteomicsIndex;
%% enzUsages


%% proteomics

