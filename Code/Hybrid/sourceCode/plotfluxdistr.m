function[] = plotfluxdistr(glucflux, sol, model, path)

%written by: Julia Münch
%date: 2019-12-11
%description: this function plots the fluxes through the
%enzymes of interest for growthrates=0.01, 0.265 and 0.6 (before, at and after induction of Crabtree effect)
%arguments: 
%   1. glucose influx
%   2. flux distribution for specific gluc influx
%   3. FBA model
%   4. path to save figures
    
    %% plot flux distributions
    if glucflux == 0.01 || glucflux == 0.265 || glucflux == 0.6
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
        
        print([path, 'fluxdistr', strrep(num2str(glucflux), '.', ','), 'growth'],'-dpng','-r0')
        close(b)
    end
    
end