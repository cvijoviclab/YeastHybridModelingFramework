function []=analyze_CLSTargets()

CLS_target1=readtable('../../data/aging-07-0177-s001.xlsx','Range','3:217');
CLS_target2=readtable('../../data/aging-07-0177-s002.xlsx','Range','3:341');
enzyme_usageMutant=readtable('../../results/enzUsages_mutants.xlsx');
CLS_target1.Properties.VariableNames{'p_value'} = 'pValue_CR';
CLS_target1.Properties.VariableNames{'FoldChange'} = 'fc_CR';
CLS_target1.Properties.VariableNames{'p_value_1'} = 'pValue_ade4D';
CLS_target1.Properties.VariableNames{'FoldChange_1'} = 'fc_ade4D';
CLS_target2.Properties.VariableNames{'p_value'} = 'pValue_CRc';
CLS_target2.Properties.VariableNames{'FoldChange'} = 'fc_CRc';
CLS_target2.Properties.VariableNames{'p_value_1'} = 'pValue_INAM';
CLS_target2.Properties.VariableNames{'FoldChange_1'} = 'fc_INAM';

%Intersect lists of targets
A=string(CLS_target1{:,1});
B=string(enzyme_usageMutant{:,3});
[C,ia,ib]=intersect(A,B);
CLS_T1=[enzyme_usageMutant(ib,1:3) CLS_target1(ia,3:6) enzyme_usageMutant(ib,5:10)];

A=string(CLS_target2{:,1});
B=string(enzyme_usageMutant{:,3});
[C,ia,ib]=intersect(A,B);
CLS_T2=[enzyme_usageMutant(ib,1:3) CLS_target2(ia,3:6) enzyme_usageMutant(ib,5:10)];
writetable(CLS_T1 ,'../../results/CLS_Targets1.txt','delimiter','\t','QuoteStrings',false);
writetable(CLS_T2 ,'../../results/CLS_Targets2.txt','delimiter','\t','QuoteStrings',false);
end