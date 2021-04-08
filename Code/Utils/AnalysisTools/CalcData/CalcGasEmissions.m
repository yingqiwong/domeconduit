function [H2OFlux, CO2Flux, vgdiff, Gamma, rho_g] = CalcGasEmissions (td, m, plot_opt)
%
% Calculates gas flux out of the conduit assuming that all gas loss in the
% plug is vertical (no lateral gas escape through conduit walls).
%
% Compares with CO2 emissions data from Gerlach et al. (2008) and Dzurisin
% et al. (2015).  Output fluxes in ton/day

tdvars = extract_y(td,m);
[phigcVec,plugdepth] = PlotResults('calc_phigc_plugdepth',td, m, 0);

% just using upwinding to estimate p and mw depth gradients at continuity
% cell centers
dpdz = diff(tdvars.p)./diff(td.z');
dpdz = [dpdz(1,:); dpdz];
dmwdz = diff(tdvars.mw)./diff(td.z');
dmwdz = [dmwdz(1,:); dmwdz];

rho_g = tdvars.p/m.T.*(tdvars.mw./m.Rw + (1-tdvars.mw)./m.Rc);
drhogdz = rho_g./tdvars.p.*dpdz + ...
    tdvars.p*(1/(m.Rw*m.T) - 1/(m.Rc*m.T)).*dmwdz;

dz = td.z(2) - td.z(1);

vgdiff = zeros(1,length(td.x));

for ti = 1:length(td.x)
    
    pd = m.PlugDepth;
    
    if m.PlugDepth==0, pd = plugdepth(ti); end
    
    inds = pd:length(td.z);
    
    % set up system of equations to solve for (vg-v)(t)
    RHS = -tdvars.phi_g(inds,ti).*tdvars.v(inds,ti).*drhogdz(inds,ti);
    
    % first point of RHS should contain vg from just below plug depth
    kmag_pd = tdcFV('calc_k_allz', m, tdvars.phi_g(pd, ti), [], td.z(pd), 1);
    vg_pd = -kmag_pd/m.eta_g*dpdz(pd,ti);
    RHS(1) = RHS(1) + 1/dz*rho_g(pd,ti)*tdvars.phi_g(pd,ti)*vg_pd;
    
    A = 1/dz*(diag(rho_g(inds,ti).*tdvars.phi_g(inds,ti),0) + ...
        diag(-rho_g(inds(1:end-1),ti).*tdvars.phi_g(inds(1:end-1),ti),-1));
    
    vgvec = A\RHS;
    vgdiff(ti) = vgvec(end);
    
end

Gamma = (1-tdvars.mw)./m.B./tdvars.mw;

% gas flux in tons/day - need to check if we need phi_g term
H2OFlux = 1e-3*24*3600*pi*m.R^2.*1./(1+Gamma(end,:)).*...
    rho_g(end,:).*tdvars.phi_g(end,:).*vgdiff;
CO2Flux = Gamma(end,:).*H2OFlux;


if plot_opt == 1
    
    % load gas emissions data from Gerlach (2008), Dzurisin (2015)
    T = readtable('GasEmissionsTimeSeries.xlsx');
    ErrFrac = 0.2;
    
    date = table2array(T(:,1));
    GasDate = datenum(datestr(date,'dd-mmm-yy',1900));
    StartDate = datenum(datestr(EruptionStartDate,'dd-mmm-yy',1900));
    GasNoDate = 1/365*(GasDate - StartDate);
    
    CO2dat    = table2array(T(:,2));
    CO2datErr = ErrFrac*table2array(T(:,2));
    
    figure;
    hdat = errorbar(GasNoDate,CO2dat,CO2datErr,'^', 'color', 0.7*ones(3,1),...
        'linewidth',1,'markersize',5,'markerfacecolor','k','markeredgecolor','k'); 
    hold on;
    hmod = plot(ConvertSecToYear(td.x), CO2Flux, 'b-', 'linewidth',2);
    hold off;
    ylabel('CO2 emission rate (ton/day)');
    xlabel('Time (years)');
    legend([hdat, hmod], {'Data', 'Model'}, 'location', 'northeast'); 
    legend boxoff;
end


end