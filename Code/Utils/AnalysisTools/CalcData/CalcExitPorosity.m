function [ExPor] = CalcExitPorosity (td, m, plot_opt)

ExPorInd = m.Nv*(m.Nz-1) + m.blk.is.phi_g;
ExPor = td.y(ExPorInd,:);

if plot_opt == 1
    % extrusion flux time series
    
    PoroDat = LoadPorosityData;
    t = mean(PoroDat.Date,2);
    dt = PoroDat.Date(:,2) - t;
    
    figure;
    errorbar(t,PoroDat.Data,PoroDat.Err,PoroDat.Err,dt,dt,...
        '^', 'color', 0.7*ones(3,1),'linewidth',1,...
        'markersize',5,'markerfacecolor','k','markeredgecolor','k');
    hold on;
    plot(ConvertSecToYear(td.x), 100*ExPor, 'b-','linewidth',2);
    hold off;
    legend('Data', 'Model', 'location', 'southeast'); legend boxoff;
    ylabel('Exit Porosity (%)');
    xlabel('Time (years)');
    title('Exit porosity');
end


end