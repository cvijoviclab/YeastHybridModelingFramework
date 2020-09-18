function[rank] = postTranscriptionalModifications2(EnzAct, model)
%written by: Julia and Linnea
%date: 2019-12-04
%description: modifies the rank in the geneExpression after post
%transcriptional modifications
%input: 1. matrix of enzyme activity
%       2. GeneExpr table
%returns: table of gene expression, and mapping of from gene to protein
%names
enzNames = model.enzNames;
enzGenes = model.enzGenes;
rank = zeros(length(enzNames),1);
pos = zeros(length(enzNames),1);
rank = table(enzNames, rank, enzGenes, pos);
%The PTM in Hxk2 is not known to have an influence on metabolic function 
%https://www.ncbi.nlm.nih.gov/pmc/articles/PMC110673/
EnzAct{6,2} = 0;

% Adding isoenzyme
if EnzAct{3,2} == 1
    EnzAct(7,:) = cell2table({"PYK2" 1});
end
for i=1:height(EnzAct)
    if contains(EnzAct{i,1}, rank{:,1})
        idx=find(contains(rank{:,1}, EnzAct{i,1}));
        if EnzAct{i,2} > 0
            rank{idx,2} = rank{idx,2} +1;
        elseif EnzAct{i,2} < 0
            rank{idx,2} = rank{idx,2} -1;
        end
    end
end
