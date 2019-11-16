function [CO2out] = CalcCO2Flux (td, m, plot_opt)

tdvars = extract_y(td,m);

rho_g = tdvars.p.*(tdvars.mw./(m.Rw*m.T) + (1 - tdvars.mw)./(m.Rc*m.T)); 
rho_gtop = rho_g(end,:);

Gamma = (1 - tdvars.mw)./(m.B*tdvars.mw);
Gammatop = Gamma(end,:);

% get pressure gradient for vertical gas velocity
dz = abs(td.z(2) - td.z(1));
dpdz = 1/dz*(tdvars.p(2:end,:) - tdvars.p(1:end-1,:));

% vertical gas velocity
vgtop = zeros(1,length(td.x));
for ti = 1:length(td.x)
    vgtop(ti) = CalcVertGasVel(td.x(ti), td.z', tdvars.phi_g(:,ti), dpdz(:,ti), m);
end

% add components that moves with magma and that moves faster than magma
% first in kg/s/m2 (i.e. per unit area)
CO2kgs = Gammatop./(1+Gammatop).*rho_gtop.*tdvars.phi_g(end,:).*vgtop + ...
     Gammatop./(1+Gammatop).*rho_gtop.*tdvars.phi_g(end,:).*tdvars.v(end,:);
 
% convert to tons/day/m2
CO2tpd = CO2kgs*24*3600/1000;

% multiply output by conduit cross sectional area
CO2out = pi*m.R^2*CO2tpd;

if plot_opt ==1
    [GasNoDate, CO2, CO2Err, ~, ~] = LoadGasEmissions;
    
    figure;
    errorbar(GasNoDate,CO2,CO2Err,'k^',...
        'linewidth',1,'markersize',5,'markerfacecolor','k');
    hold on;
    plot(ConvertSecToYear(td.x), CO2out, 'b-');
    hold off;
    legend('Data', 'Model', 'location', 'southeast');
    ylabel('CO2 emission rate (ton/day)');
    xlabel('Time (years)');
end

end


function [vgtop] = CalcVertGasVel (t, z, phi_g, dpdz, m)

% calculate vertical permeability
zphigcInd = tdcFV('calc_zphigc', t, phi_g, m);
degas  = 1./(1+exp(-0.04*(z-z(zphigcInd))));
[kvert, ~] = tdcFV('calc_k_allz', m, phi_g, z, degas);

% get values at top of conduit (north face of final cell)
[kmagn, ~] = tdcFV('get_FaceVals', kvert/m.eta_g, m.dQUICKn);
dpdzn = [dpdz; 2*dpdz(end) - dpdz(end-1)];
vgdiff = -kmagn.*dpdzn;
vgtop = vgdiff(end);

end