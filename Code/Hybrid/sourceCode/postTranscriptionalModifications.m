function[GeneExpr] = postTranscriptionalModifications(EnzAct, GeneExpr)
%written by: Julia M�nch
%date: 2019-12-04
%description: modifies the rank in the geneExpression after post
%transcriptional modifications
%input: 1. matrix of enzyme activity
%       2. GeneExpr table
%returns: table of gene expression, and mapping of from gene to protein
%names
for i=1:height(EnzAct)
    if contains(EnzAct{i,1}, GeneExpr{:,1})
        idx=find(contains(GeneExpr{:,1}, EnzAct{i,1}));
        if EnzAct{i,2} > 0
            GeneExpr{idx,2} = GeneExpr{idx,2} +2;
        elseif EnzAct{i,2} < 0
            GeneExpr{idx,2} = GeneExpr{idx,2} -2;
        end
    end
end