function[ExpressionChange, EnzActChange, o2, co2, eth, ac, glcin, growth] = runMultiscale(WT_model, knockouts, activeCrosstalk, gluc, nitr, Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, perc, path)

%written by: Julia Münch
%date: 2019-12-11
%description: this function iterates over increasing growthrates and turns
%             on the boolean signalling model if growthrate > 0.26 (start
%             short-term Crabtree effect); simulated changes in gene
%             expression are used to modify turnover number of enzymes in
%             FBA model
%input: 1. WT model
%       2. knockouts
%       3. activeCrosstalk
%       4.-5. vectors containing glucose and nitrogen availability
%       6.-11. matrices with states of pathway components
%       12. percentage of how much of its value kcat should be in/decreased 
%       13. path to save results
%returns: 1.-2. change in gene expression and enzyme activity
%         3.-8. respective fluxes for o2in, co2out, ethout, acout, glcin
%         and growth
%contains functions:
%   1. knockouts  
%   2. reachSteadyState
%      contains 3. crosstalk
%               4. activityConverter
%               5. TFtargets   
%      6. Expression
%      7. BoolToFBA2
%         contains 8. getMutant
%                  9. plotfluxdistr
%                  10. plotfluxes

%% iterate over growth rates and run boolmodel at each iteration after critical growthrate (0.26) has been reached
x_val = 0:0.005:5;
eth=[];
o2=[];
ac=[];
co2=[];
glcin=[];
growth=[];

i=1;

while true
    
    %% iteration over growth without bool activity
    if x_val(i) < 0.26 %3.3
        model = WT_model;
        model = setParam(model, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn', 'GROWTH'}, [1.3, 1.7, 0, 0, 0, 0, x_val(i)]);
        model = setParam(model, 'obj', {'GROWTH'}, 1);
        sol = solveLP(model,1);
        fluxmodel = model;
        
    % turn on bool activity (modelled as (0|0) -> (1|1))
    elseif x_val(i) > 0.26 %3.3

        for k = 1:length(gluc)
            disp(['Glucose level: ', num2str(gluc(k)), ', Nitrogen level: ', num2str(nitr(k))]);

            %run the boolean model to create txt files
            if k > 1
                knockouts = knockouts;
                [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets] = knockout(Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, knockouts);        
            end

            [Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, TFActStart, TFActEnd, EnzActStart, EnzActEnd] = ...
                reachSteadyState(Metabolites, PKApw, Snf1pw, TORpw, Enzymes, Targets, gluc(k), nitr(k), path, activeCrosstalk);
        end

        %% convert TF activities into gene expression changes
        [ExpressionStart] = Expression(TFActStart, WT_model);
        [ExpressionEnd] = Expression(TFActEnd, WT_model);

        % calculate difference in Gene Expression and Enzyme Activity
        ExpressionChange = ExpressionStart;
        ExpressionChange{:,2} = ExpressionEnd{:,2} - ExpressionStart{:,2};
        EnzActChange = EnzActStart;
        EnzActChange{:,2} = EnzActEnd{:,2} - EnzActStart{:,2};

        
        % for first time boolean is acitve, mutant is based on WT model
        if x_val(i) < 0.27 %3.5
        mut = BoolToFBA2(ExpressionChange, EnzActChange, model, perc);
        
        % for further bool modification steps, mutant is based on mutant at
        % iteration before
        elseif x_val(i) >= 0.27 %3.5
            mut = BoolToFBA2(ExpressionChange,EnzActChange, mut, perc);
        end
        mut = setParam(mut, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn', 'GROWTH'}, [1.3, 1.7, 0, 0, 0, 0, x_val(i)]);
        mut = setParam(mut, 'obj', {'GROWTH'}, 1);
        sol = solveLP(mut,1);
        fluxmodel=mut;
               
    end
    
    plotfluxdistr(x_val(i), sol, fluxmodel, path)
    
    eth(i) = sol.x(5);
    o2(i) = sol.x(8);
    ac(i)=sol.x(1);
    co2(i)=sol.x(4);
    growth(i) = sol.x(40);
    glcin(i) = sol.x(7);
        
    
    if i > 2
        if growth(i) <= growth(i-2)
            break
        else
            disp(['Iteration ', num2str(i), ' completed!']);
        end
    end
    
    i=i+1;    
end

plotfluxes (o2, co2, eth, ac, glcin, growth, path)

end