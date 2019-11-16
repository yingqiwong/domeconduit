function [VolExRate] = CalcExtrusionRate (td, m, plot_opt)

% calculates the rate of volume extruded (m3/s, dense rock equivalent)

tdvars = extract_y(td,m);

% calculate dense rock equivalent of rate of extruded volume by removing 
% gas vol fraction
VolExRate = pi*m.R^2*tdvars.v(end,:).*(1-tdvars.phi_g(end,:));

if plot_opt == 1
    % extrusion flux time series
    
    [VolRate] = LoadExtrusionRate;

    figure;
    errorbar(VolRate.Date,VolRate.ExtrRate,VolRate.ExtrRateErr,'^', ...
        'color', 0.7*ones(3,1),...
        'linewidth',1,'markersize',5,'markerfacecolor','k','markeredgecolor','k');
    hold on;
    plot(tdvars.tyr, VolExRate, 'b-','linewidth',2);
    hold off;
    legend('Data', 'Model', 'location', 'southeast'); legend boxoff;
    % datetick('x','mmmyy');
    ylabel('Extruded lava volume (x10^6 m^3)');
    xlabel('Time (years)');
    title('Extruded volume');
end


end