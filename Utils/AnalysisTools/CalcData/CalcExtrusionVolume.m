function [VolEx] = CalcExtrusionVolume (td, m, plot_opt)

% calculates the volume extruded using trapezidal rule for exit velocity
% in units of million m3 (10^6 m3, dense rock equivalent)

tdvars = extract_y(td,m);

% calculate dense rock equivalent of volume extruded by removing gas vol
% fraction
flux = 1e-6*pi*m.R^2*tdvars.v(end,:).*(1-tdvars.phi_g(end,:));
VolEx = cumtrapz(td.x, flux);

if plot_opt == 1
    % extrusion flux time series
    
    [FluxNoDate, CumVol, CumVolErr] = LoadExtrusionVol ();

    figure;
    errorbar(FluxNoDate,CumVol,CumVolErr,'^', 'color', 0.7*ones(3,1),...
        'linewidth',1,'markersize',5,'markerfacecolor','k','markeredgecolor','k');
    hold on;
    plot(ConvertSecToYear(td.x), VolEx, 'b-','linewidth',2);
    hold off;
    legend('Data', 'Model', 'location', 'southeast'); legend boxoff;
    % datetick('x','mmmyy');
    ylabel('Extruded lava volume (x10^6 m^3)');
    xlabel('Time (years)');
    title('Extruded volume');
end

end