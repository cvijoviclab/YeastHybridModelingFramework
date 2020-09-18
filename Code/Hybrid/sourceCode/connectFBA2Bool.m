function [Bool,settings] = connectFBA2Bool(Bool,settings,conBool,i)

% based on glucose uptake rate
if ~isempty(conBool)
%    if table2array(conBool(1,i)) > 3.2914
    if conBool(1,i) > 3.2914
        settings.gluc(1)=1;
        %settings.enzymeUse = (0.055);
    else
        settings.gluc(1)=0;
    end
end

%based on hexokinase

%based on FBP1

%based on ATP productionrate