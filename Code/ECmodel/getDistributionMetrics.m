function dist_metrics = getDistributionMetrics(fluxTable,variables,pathways,model)
if nargin<4
    model = [];
    if nargin<3
        pathways = '';
    end
end
D_rate  = [];
nZeros  = [];
maxVal  = [];
minVal  = [];
meanVal = [];
medianVal = [];
stddev    = [];
fermentation = [];
glyc_meanFlux   = [];
oxPhos_meanFlux = [];
TCA_meanFlux    = [];
PPP_meanFlux    = [];
ana_meanFlux    = [];

glyc_burden   = [];
oxPhos_burden = [];
TCA_burden    = [];
PPP_burden    = [];
ana_burden    = [];

for i=1:width(fluxTable)
    colName = fluxTable.Properties.VariableNames(i);
    if contains(colName,'D_')
        gRate = str2double(strrep(colName,'D_','0.'));
        eval(['distribution = fluxTable.' colName{1} ';'])
        D_rate  = [D_rate;gRate];
        nZeros  = [nZeros; sum(distribution==0)];
        %Get non-zero fluxes
        nonZeros   = distribution~=0;
        subSystems = fluxTable.subSystems(nonZeros);
        newDist  = distribution(nonZeros);
        maxVal   = [maxVal; max(newDist)];
        minVal   = [minVal; min(newDist)];
        meanVal  = [meanVal; mean(newDist)];
        medianVal = [medianVal; median(newDist)];
        stddev    = [stddev; std(newDist)];
        fermentation = [fermentation; gRate>=0.285];
        if ~isempty(pathways)
            for j=1:length(pathways)
                path = pathways(j);
                p_indexes = find(contains(subSystems,path));
                if isempty(p_indexes)
                    newDist = 0;
                else
                   dist = newDist(p_indexes); 
                end
                switch j
                    case 1
                        glyc_meanFlux = [glyc_meanFlux;mean(dist)];
                        glyc_burden = [glyc_burden;sum(dist)];
                    case 2
                        oxPhos_meanFlux = [oxPhos_meanFlux;mean(dist)];
                        oxPhos_burden = [oxPhos_burden;sum(dist)];
                    case 3
                        TCA_meanFlux = [TCA_meanFlux;mean(dist)];
                        TCA_burden = [TCA_burden;sum(dist)];
                    case 4
                        PPP_meanFlux = [PPP_meanFlux; mean(dist)];
                        PPP_burden = [PPP_burden; sum(dist)];
                    case 5
                        ana_meanFlux = [ana_meanFlux; mean(dist)];
                        ana_burden = [ana_burden; sum(dist)];
                end
            end
        end
    end
end
if isempty(model)
    dist_metrics = table(D_rate,nZeros,maxVal,minVal,meanVal,medianVal,stddev,...
                     fermentation,glyc_meanFlux,oxPhos_meanFlux,TCA_meanFlux,...
                     PPP_meanFlux,'VariableNames',variables);
else
    dist_metrics = table(D_rate,nZeros,maxVal,minVal,meanVal,medianVal,stddev,...
                     fermentation,glyc_meanFlux,oxPhos_meanFlux,TCA_meanFlux,...
                     PPP_meanFlux,ana_meanFlux,glyc_burden,oxPhos_burden,TCA_burden,PPP_burden,ana_burden,...
                     'VariableNames',variables);
end
end
