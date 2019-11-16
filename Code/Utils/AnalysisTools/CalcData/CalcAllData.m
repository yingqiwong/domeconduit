function [VolExtd, JRO1td, CO2td] = CalcAllData (td, m, plot_opt, EruptDate, VolEx, JRO1, GasEm)

[VolExtd] = CalcExtrusionVolume(td, m, 0);
[ur, ut, uz] = CalcJRO1Def(td, m, 0); 
JRO1td = 1e3*ur;
[~, CO2td, ~] = CalcGasEmissions(td, m, 0);




if (plot_opt)
    figure;
    set(gcf,'Position',[398 284 1389 354]);
    
    subplot(131);
    errorbar(VolEx.Date,VolEx.Vol,VolEx.Err,'^', 'color', 0.7*ones(3,1),...
        'linewidth',1,'markersize',5,'markerfacecolor','k','markeredgecolor','k');
    hold on;
    plot(ConvertSecToYear(td.x), VolExtd, 'b-','linewidth',2);
    hold off;
    legend('Data', 'Model', 'location', 'southeast'); legend boxoff;
    % datetick('x','mmmyy');
    ylabel('Extruded lava volume (x10^6 m^3)');
    xlabel('Time (years)');
    title('Extruded volume');

    subplot(132);
    h1 = errorbar(JRO1.Date - EruptDate.Start, ...
         JRO1.rad_clean, 0.5*JRO1.rad_error, ...
         '.','color', 0.8*ones(3,1),'linewidth',1,...
         'markersize',5,'markerfacecolor', 'k', 'markeredgecolor', 'k', 'CapSize', 0);
    hold on; xlim([-1,4]); ylim auto; 
    h2 = plot(ConvertSecToYear(td.x), JRO1td, 'b-','linewidth',2);
    plot([0,0], ylim, 'k-');
    hold off;
    title('JRO1 radial deformation');
    xlabel('Time (years)'); ylabel('Radial displacements (mm)');
    legend([h1 h2], {'Data', 'Model'}); legend boxoff;
    
    subplot(133);
    hdat = errorbar(GasEm.Date,GasEm.Flux,GasEm.Err,'^', 'color', 0.7*ones(3,1),...
        'linewidth',1,'markersize',5,'markerfacecolor','k','markeredgecolor','k'); 
    hold on;
    hmod = plot(ConvertSecToYear(td.x), CO2td, 'b-', 'linewidth',2);
    hold off;
    ylabel('CO2 emission rate (ton/day)');
    xlabel('Time (years)');
    legend([hdat, hmod], {'Data', 'Model'}, 'location', 'northeast'); 
    legend boxoff;
    title('CO2 flux');
    
    suptitle({...
    ['V0 = ', num2str(1e-9*m.ch.V0), 'km3, ',...
    'Omega = ', num2str(m.ch.Omega*3600*1e6,4), 'm3/day/MPa, ',...
    'Aspect Ratio = ', num2str(m.ch.AR, 2), ', ', ...
    'pch0 = ', num2str(m.ch.pdeep*1e-6, 4), 'MPa, ', ...
    'radius = ', num2str(m.R), 'm, ', ...
    'f0 = ', num2str(m.fr.f0), ', ', ...
    'a = ', num2str(m.fr.a), ','], ...
    ['total h2o = ', num2str(m.chi_ch.total.h2o), 'wt%, ', ...
    'total co2 = ', num2str(m.chi_ch.total.co2), 'ppm, ', ...
    'kc = ', num2str(m.k_lat), 'm2, ', ...
    'phigc = ', num2str(m.phi_gc)]});
end





end

