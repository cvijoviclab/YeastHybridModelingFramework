function [results,breakFlag] = get_ECCs(model,target,perturbation,tolerance,direction)
%get_ECCs
%
%
% Ivan Domenzain    Last edited. 2020-08-13

if nargin <5
    direction = 'descend';
end
%Get an initial simulation optimizing for biomass production
base_sol = solveLP(model);
%Extract the metabolite indexes corresponding to the model enzymes
%(avoiding the protein pool)
results   = [];
breakFlag = false;
%If the model is feasible, analyse the individual kcat values of the flux
%carrier enzymes
if ~isempty(base_sol.f)
    %First, find the limiting Kcats in the model
    ECCs_table = computeECCs(model,base_sol,target,perturbation,tolerance);
    if ~isempty(ECCs_table.uniprot)
        %Sort the limiting Kcats cell according to their control coefficient
        ECCs_table = sortrows(ECCs_table,'ECC',direction);
    end
    
    if ~isempty(ECCs_table.uniprot)
        results = ECCs_table;
    else
        breakFlag = true;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function ECCstable = computeECCs(model,base_sol,target,perturbation,tolerance)
ECCstable = cell(1,10);
%Extract the enzyme usage pseudoreactions indexes
enzUsageIndxs = find(~cellfun(@isempty,strfind(model.rxns,'draw_prot_')));
%Get the original rxn indexes of the flux carrier enzymes (enzyme usages)
enzUsageIndxs = enzUsageIndxs(base_sol.x(enzUsageIndxs)>0);
enzUsageFlux  = base_sol.x(enzUsageIndxs);
%Fix all enzyme usages
%model.lb(enzUsageIndxs) = (1-tolerance)*enzUsageFlux;
%model.ub(enzUsageIndxs) = (1+tolerance)*enzUsageFlux;
%Extract the metabolite index for the flux carrying enzymes
usedEnzIndxs = zeros(size(enzUsageIndxs));
for i=1:length(usedEnzIndxs)
    usedEnzIndxs(i) = find(model.S(:,enzUsageIndxs(i))==1);
end
%For each flux carrying enzyme calculate the objective control
%coefficient (ObjCC) on each of its kinetic coefficients
for i = 1:length(usedEnzIndxs)
    option = 1;
    enzPos  = usedEnzIndxs(i);
    enzName = model.mets(enzPos);
    %Get the nonzero kinetic coefficients
    enz_rxnsCoeffs = find(model.S(enzPos,:));
    %For each enzymatic reaction (avoid enzyme usage pseudoRxn)
    if ~isempty(enz_rxnsCoeffs)
        for j=1:length(enz_rxnsCoeffs)-1
            temp_model = model;
            %Allow perturbation of the i-th enzyme
            temp_model.lb(enzUsageIndxs(i)) = 0;
            temp_model.ub(enzUsageIndxs(i)) = 1000;
            coeffPos = enz_rxnsCoeffs(j);
            kcoeff   = model.S(enzPos,coeffPos);
            if option == 1
                temp_model = setParam(temp_model,'lb',enzUsageIndxs(i),(1-tolerance)*perturbation*enzUsageFlux(i));
                temp_model = setParam(temp_model,'ub',enzUsageIndxs(i),(1+tolerance)*perturbation*enzUsageFlux(i));
                new_sol = solveLP(temp_model);
                factor = 1;
                if isempty(new_sol.f)
                    option = 2;
                end
            end
            if option == 2
                temp_model.S(enzPos,coeffPos) = kcoeff/perturbation;
                factor = perturbation;
            end
            a_i = kcoeff*enzUsageFlux(i);
            %Get a new solution with flexibilized j-th Kcat for the i-th
            %enzyme
            new_sol = solveLP(temp_model);
            if ~isempty(new_sol.f)
                a_i_new = kcoeff*factor*new_sol.x(enzUsageIndxs(i));
                delta_V  = (new_sol.x(target) - base_sol.x(target))/base_sol.x(target);
                delta_Ea = (a_i_new-a_i)/a_i;
                %if (a_i_new*delta_Ea*delta_V)~=0 %& 
                if delta_Ea>0
                    %If solution was feasible then calculate the control
                    %coefficient for the j-th Kcat
                    objCC = delta_V/delta_Ea;
                    if abs(objCC) > (tolerance/(perturbation-1))
                        Kcat         = (-1/model.S(enzPos,coeffPos))/3600;
                        ECCstable{1} = [ECCstable{1}; enzName];
                        ECCstable{2} = [ECCstable{2}; enzPos];
                        ECCstable{3} = [ECCstable{3}; coeffPos];
                        ECCstable{4} = [ECCstable{4}; Kcat];
                        ECCstable{5} = [ECCstable{5}; objCC];
                        ECCstable{6} = [ECCstable{6}; model.rxnNames(coeffPos)];
                        ECCstable{7} = [ECCstable{7}; model.rxns(coeffPos)];
                        ECCstable{8} = [ECCstable{8}; delta_V];
                        ECCstable{9} = [ECCstable{9}; delta_Ea];
                        ECCstable{10}= [ECCstable{10}; option];
                    end
                end
            end
        end
    end
end
vars = {'uniprot' 'rxnName' 'Kcat' 'ECC' 'rxnID' 'delta_V' 'delta_Ea' 'option'};
ECCstable = table(ECCstable{1},ECCstable{6},ECCstable{4},ECCstable{5},ECCstable{7},ECCstable{8},ECCstable{9},ECCstable{10},'VariableNames',vars);
%Clean ECCs table
rxns = unique(ECCstable.rxnID);
for i=1:length(rxns)
    rxn      = rxns(i);
    idxs     = find(strcmpi(ECCstable.rxnID,rxn));
    [~,MI]   = max(ECCstable.ECC(idxs));
    idxs(MI) = [];
    ECCstable(idxs,:) = [];
end
end
