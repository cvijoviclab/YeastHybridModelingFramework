function[] = plotfluxes (o2, co2, eth, ac, glcin, growth, path)

%written by: Julia Münch
%date: 2019-12-11
%description: this function plots the fluxes for specific nutrients against
%the respective growthrate and the growthrate against the glucose influx
%arguments: 
%   1.-6. vectors containing fluxes
%   7. path to save figures
    
%%  plot nutrient fluxes against growthrate
    figure
    plot(growth(1:end-1), o2(1:end-1), 'Color', [0, 0.4470, 0.7410], 'LineWidth', 1.5, 'DisplayName', 'O_2');%blue
    hold on
    plot(growth(1:end-1), co2(1:end-1), 'Color', [0.6150, 0.5150, 0.7450], 'LineWidth', 1.5, 'DisplayName', 'CO_2'); %purple
    plot(growth(1:end-1), eth(1:end-1), 'Color', [0.6350, 0.0780, 0.1840], 'LineWidth', 1.5, 'DisplayName', 'ethanol'); %red
    plot(growth(1:end-1), ac(1:end-1), 'Color', [0.4660, 0.6740, 0.1880], 'LineWidth', 1.5, 'DisplayName', 'acetate'); %green
    plot(growth(1:end-1), glcin(1:end-1), 'Color', [0.9290, 0.6940, 0.1250], 'LineWidth', 1.5, 'DisplayName', 'glucose'); %yellow

    %chose axis label

    xlabel('Growthrate [h^{-1}]');
    ylabel('flux [mMol/h/g dw]');

    %place legend and turn box off
    legend('Location', 'west', 'FontSize', 14);
    legend('boxoff');

    % chose axis color, turn box off, chose font size
    f = gcf;
    f.CurrentAxes.Box = 'off';
    f.CurrentAxes.XColor = [0 0 0];
    f.CurrentAxes.YColor = [0 0 0];
    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.XLim = [0 0.7];
    
    print([path, 'fluxes'],'-dpng','-r0')
    close(f)
    
    %% plot growthrate against glucose influx    
    figure
    plot(glcin(1:end-1), growth(1:end-1), 'Color', [0.8500, 0.3250, 0.0980],'LineWidth', 1.5)

    %chose axis label
    xlabel('glucose influx [mMol/h/g dw]');
    ylabel('Growthrate [h^{-1}]');

    % chose axis color, turn box off, chose font size
    f = gcf;
    f.CurrentAxes.Box = 'off';
    f.CurrentAxes.XColor = [0 0 0];
    f.CurrentAxes.YColor = [0 0 0];
    f.CurrentAxes.FontSize = 14;
    f.CurrentAxes.XLim = [0 max(glcin)+2];
    f.CurrentAxes.YLim = [0 max(growth)+0.1];
    
    print([path, 'GlucVsGrowth'],'-dpng','-r0')
    close(f)
    
end
