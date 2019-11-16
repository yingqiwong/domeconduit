function [VolErr, JRO1DefErr, CO2Err] = CalcDataModelErr (td, m, EruptDate, VolEx, JRO1, GasEm, plot_opt)

tyr = ConvertSecToYear(td.x);
EruptDate.Span = EruptDate.End - EruptDate.Start;

% extrusion volume error
[VolExtd] = CalcExtrusionVolume(td, m, 0);
PredVol = interp1(tyr, VolExtd, VolEx.Date);
VolInds = find(VolEx.Date<=EruptDate.Span); % pick out data during eruption only to calculate error
VolErr = sum(((VolEx.Vol(VolInds) - PredVol(VolInds))./VolEx.Err(VolInds)).^2, 'omitnan');
% pick out data during eruption only to calculate error

% predict radial displacement
[ur, ut, uz] = CalcJRO1Def(td, m, 0); 
JRO1td = 1e3*ur;
PredDefm = interp1(tyr, JRO1td, JRO1.Date-EruptDate.Start);
DefmInds = find(JRO1.Date>=EruptDate.Start & JRO1.Date<=EruptDate.End);
JRO1DefErr = sum(((JRO1.rad_clean(DefmInds) - PredDefm(DefmInds))./JRO1.rad_error(DefmInds)).^2,'omitnan');

% predict CO2 flux
[~, CO2td, ~] = CalcGasEmissions(td, m, 0);
PredCO2 = interp1(tyr, CO2td, GasEm.Date);
CO2Inds = find(GasEm.Date<=EruptDate.Span);
CO2Err = (sum(1./GasEm.Err(CO2Inds).^2.*(GasEm.Flux(CO2Inds) - PredCO2(CO2Inds)).^2, 'omitnan'));




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
    title(sprintf('Extruded volume, Err = %.2f', VolErr));

    subplot(132);
    h1 = errorbar(JRO1.Date - EruptDate.Start, ...
         JRO1.rad_clean, 0.5*JRO1.rad_error, ...
         '.','color', 0.8*ones(3,1),'linewidth',1,...
         'markersize',5,'markerfacecolor', 'k', 'markeredgecolor', 'k', 'CapSize', 0);
    hold on; xlim([-1,4]); ylim auto; 
    h2 = plot(ConvertSecToYear(td.x), JRO1td, 'b-','linewidth',2);
    plot([0,0], ylim, 'k-');
    hold off;
    title(sprintf('JRO1 defm, Err = %.2f', JRO1DefErr));
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
    title(sprintf('CO2 flux, Err = %.2f', CO2Err));
    
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

