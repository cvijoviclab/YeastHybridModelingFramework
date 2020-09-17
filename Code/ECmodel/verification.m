%% model verification of metabolic model -> reproduce fig 2b from Nilsson paper
%written by: Julia M�nch
%date: 2019-12-04
%description: this script reproduces the results by Nilsson (nutrient flux dependency on growthrate/dilution) using the
%             converted ECmodel; in addition to showing dependency on
%             growthrate, dependency on glucose influx is also depicted and
%             the resulting flux distributions are depicted as bar plots
%%             
load('reduced_ecYeast_batch.mat');
model = ecModel_batch;

for j=1:2
    x = {'glcIN', 'GROWTH'};

    eth=[];
    o2=[];
    ac=[];
    co2=[];
    glcin=[];
    growth=[];
    
    %max glucose uptake rate approx. 20
    if j == 1
        x_val=0:0.5:25;

    %max growth rate approx. 0.4
    elseif j == 2
        x_val = 0:0.01:1;
    end

    %iterate over fluxes and extract fluxes of interest
    %upper bound set for acetate (1.3) and glycerol (1.7) (see Nilsson paper)
    i=1;
    while true
 
        model = ecModel_batch;
        model = setParam(model, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn' x{j}}, [1.3, 1.7, 0, 0, 0, 0, x_val(i)]);
        model = setParam(model, 'obj', {'GROWTH'}, 1);
        sol = solveLP(model,1);
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
    
    for k=1:2
        xvar = {glcin, growth};

        %% create figure (flux dependency on j=1 -> glucose influx or j=2 -> growthrate)

        figure
        plot(xvar{k}, o2, 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'O_2');%blue
        hold on
        plot(xvar{k}, co2, 'Color', [0.6150, 0.5150, 0.7450], 'LineWidth', 1.5, 'DisplayName', 'CO_2'); %purple
        plot(xvar{k}, eth, 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.5, 'DisplayName', 'ethanol'); %red
        plot(xvar{k}, ac, 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'acetate'); %green
    %     if j == 1
    %         plot(xvar{k}, growth, 'Color', [0.8500, 0.3250, 0.0980],
    %         'LineWidth', 1.5, 'DisplayName', 'Growthrate [h^{-1}]');%orange
        if k == 2
            plot(xvar{k}, glcin, 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5, 'DisplayName', 'glucose'); %yellow
        end

        %chose axis label
        if k == 1
            xlabel('glucose influx [mMol/h/g dw]');
        elseif k == 2
            xlabel('Growthrate [h^{-1}]');
        end
        ylabel('flux [mMol/h/g dw]');

        %place legend and turn box off
        legend('Location', 'west', 'FontSize', 14);
        legend('boxoff');

        % chose axis color, turn box off, chose font size
        f = gcf;

        if k == 1
            f.CurrentAxes.XLim = [0 23];
        elseif k == 2
            f.CurrentAxes.XLim = [0 0.7];
            f.CurrentAxes.YLim = [0 70];
        end
        
        f.CurrentAxes.Box = 'off';
        f.CurrentAxes.XColor = [0 0 0];
        f.CurrentAxes.YColor = [0 0 0];
        f.CurrentAxes.FontSize = 14;

        %save
        iter = {'glucIt', 'growthIt'};
        name = {'fluxVsGlcIn', 'fluxVsGrowth'};
        mkdir Verification
        print(['Verification/2', name{k}, iter{j}], '-dpng','-r0')
        close(f)
    end
    
    %% plot growthrate against glucose influx
    figure
    plot(glcin, growth, 'Color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5)

    %chose axis label
    xlabel('glucose influx [mMol/h/g dw]');
    ylabel('Growthrate [h^{-1}]');

    % chose axis color, turn box off, chose font size
    f = gcf;
    f.CurrentAxes.Box = 'off';
    f.CurrentAxes.XColor = [0 0 0];
    f.CurrentAxes.YColor = [0 0 0];
    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.XLim = [0 38];
    f.CurrentAxes.YLim = [0 0.75];

    %save
    iter = {'glucIt', 'growthIt'};
    print(['Verification/2GrowthVsGlc', iter{j}], '-dpng','-r0')
    close(f)

end

%% Plot flux distribution
growth = [0.01, 0.265, glcin(end)];
for i = 1:length(growth)
    model = ecModel_batch;
    model = setParam(model, 'ub', {'acOUT', 'glyOUT', 'ethIN', 'acIN', 'galIN', 'gluActiveIn' 'GROWTH'}, [1.3, 1.7, 0, 0, 0, 0, growth(i)]);
    model = setParam(model, 'obj', {'GROWTH'}, 1);
    sol = solveLP(model,1);
    
    figure
    b=bar(categorical(model.enzNames), sol.x(find(contains(model.rxns, model.enzymes))));
    b.FaceColor = [0.6,0.6,0.6];
    ylabel('Flux [mMol/h/g dw]')
    b = gcf;
    %b.CurrentAxes.DataAspectRatio = [0.1 1 1];
    b.CurrentAxes.PlotBoxAspectRatio = [30 4 1];
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96])
    b.CurrentAxes.FontSize = 9.5;
    b.CurrentAxes.TickLength = [0 0];
    b.CurrentAxes.Box = 'off';
    b.CurrentAxes.XColor = [0 0 0];
    b.CurrentAxes.YColor = [0 0 0];
    
    print(['Verification/FluxDistr', num2str(i), 'glcIN'], '-dpng','-r0')
    close(b)
end
    


