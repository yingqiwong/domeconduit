function varargout = PlotResults (varargin)
% Code to plot results from tdcFV, with several different plotting options

[varargout{1:nargout}] = feval(varargin{:});
end

function AllPlots (ss, td, m)
% PlotResults('AllPlots', ss, td, m);

z = -ss.z(1:2:end);

ss = smf_rad_dz('add_fields',ss);
rsmf_io('PlotResults',ss,500);

plot_td(td, m);

if td.x(end) == m.Tspan(2)
    plot_profiles(ss,td,m,0.8*m.Tspan(2));
else 
    plot_profiles(ss,td,m,td.x(end));
end

[phigcVec,plugdepth] = PlotResults('calc_phigc_plugdepth',td, m, 1);

end

function plot_td (td, m)
% PlotResults('plot_td', td, m);

is = m.blk.is;
z = td.z;
t = ConvertSecToYear(td.x);
p = td.y(is.p:m.Nv:end,:)*m.slv.sy(is.p);
v = td.y(is.v:m.Nv:end,:)*m.slv.sy(is.v);
phi_g = td.y(is.phi_g:m.Nv:end,:)*m.slv.sy(is.phi_g);

if m.Nv == 4
    mw = td.y(is.mw:m.Nv:end,:)*m.slv.sy(is.mw);
else
    mw = ones(size(p));
end

figure;
set(gcf, 'position', [465 48 1254 908]);
subplot(341); plot(t,1e-6*p(end,:)); title('Pressure at top (MPa)');
subplot(342); plot(t,v(end,:)); title('Velocity at top (m/s)'); 
    xlabel('Time (year)');
subplot(343); plot(t,100*phi_g(end,:)); title('Porosity at top (%)');
subplot(344); plot(t,mw(end,:)); title('Mole fraction water at top');
subplot(345); plot(t,1e-6*p(1,:)); title('Pressure at base (MPa)');
subplot(346); plot(t,v(1,:)); title('Velocity at base (m/s)'); 
    xlabel('Time (year)');
subplot(347); plot(t,100*phi_g(1,:)); title('Porosity at base (%)');
subplot(348); plot(t,mw(1,:)); title('Mole fraction water at base');
tmid = ConvertSecToYear(0.25*m.Tspan(2));
tmidi = find(t<tmid); tmidi = tmidi(end);
subplot(349); plot(1e-6*p(:,1),z,1e-6*p(:,tmidi),z,1e-6*p(:,end),z);
title('Pressure (MPa)'); ylabel('z (m)');
subplot(3,4,10); plot(v(:,1),z,v(:,tmidi),z,v(:,end),z);
title('Velocity (m/s)');
legend('t=0s','t~quarter','t=end','location','southeast');
subplot(3,4,11); plot(100*phi_g(:,1),z,100*phi_g(:,tmidi),z,100*phi_g(:,end),z);
title('Porosity');
subplot(3,4,12); plot(mw(:,1),z,mw(:,tmidi),z,mw(:,end),z);
title('Mole fraction water');

PlotTitle(m, 'Conduit properties in time')


end

function plot_timeseries (td, m)
% PlotResults('plot_timeseries', td, m);

t = ConvertSecToYear(td.x);
tdvars = extract_y(td,m);

figure;
set(gcf, 'position', [78 442 1580 469]);
subplot(241); plot(t,1e-6*tdvars.p(end,:)); 
title('Pressure at top (MPa)'); xlabel('Time (year)');
subplot(242); plot(t,tdvars.v(end,:)); 
title('Velocity at top (m/s)'); xlabel('Time (year)');
subplot(243); plot(t,100*tdvars.phi_g(end,:)); 
title('Porosity at top (%)'); xlabel('Time (year)');
subplot(244); plot(t,tdvars.mw(end,:)); 
title('Mole fraction water at top'); xlabel('Time (year)');
subplot(245); plot(t,1e-6*tdvars.p(1,:)); 
title('Pressure at base (MPa)'); xlabel('Time (year)');
subplot(246); plot(t,tdvars.v(1,:));
title('Velocity at base (m/s)'); xlabel('Time (year)');
subplot(247); plot(t,100*tdvars.phi_g(1,:)); 
title('Porosity at base (%)'); xlabel('Time (year)');
subplot(248); plot(t,tdvars.mw(1,:)); 
title('Mole fraction water at base'); xlabel('Time (year)');

PlotTitle(m, 'Conduit properties in time')


end


function plot_dydt (td, m)
% PlotResults('plot_dydt', td, m);

t = ConvertSecToYear(td.x);
tdvars = extract_y(td,m);

figure;
set(gcf, 'position', [78 406 1355 505]);
subplot(241); plot(t,tdvars.pp(end,:)); 
title('(a) dpdt at top (Pa/s)'); xlabel('Time (year)');
subplot(242); plot(t,tdvars.vp(end,:)); 
title('(b) dvdt at top (m/s2)'); xlabel('Time (year)');
subplot(243); plot(t,tdvars.phi_gp(end,:)); 
title('(c) dphig/dt at top (/s)'); xlabel('Time (year)');
subplot(244); plot(t,tdvars.mwp(end,:)); 
title('(d) dmwdt at top (/s)'); xlabel('Time (year)');
subplot(245); plot(t,tdvars.pp(1,:)); 
title('(e) dpdt at base (Pa/s)'); xlabel('Time (year)');
subplot(246); plot(t,tdvars.vp(1,:));
title('(f) dvdt at base (m/s2)'); xlabel('Time (year)');
subplot(247); plot(t,tdvars.phi_gp(1,:)); 
title('(g) dphigdt at base (/s)'); xlabel('Time (year)');
subplot(248); plot(t,tdvars.mwp(1,:)); 
title('(h) dmwdt at base (/s)'); xlabel('Time (year)');

PlotTitle(m, 'Rate of change of conduit properties in time')


end

function plot_profiles (ss, td, m, teval)
% PlotResults('plot_profiles',ss,td,m,0.8*m.Tspan(2));

z = -ss.z(1:2:end);
zMBE = z(1:end-1) + 0.5*(z(2)-z(1));
E = CalcFieldsAtOneTimeStep(td, teval, m);
if isempty(E), return; end % solution not evaluated at this time.

figure;
set(gcf, 'Position', [106 174 678 808]);
subplot(241); plot(1e-6*E.p, z); xlabel('Pressure (MPa)');
subplot(242); plot(1e-6*E.tauR, zMBE); xlabel('Shear Stress (MPa)');
subplot(243); plot(E.v.vv, zMBE, E.v.vf, zMBE, E.v.v(2:end), zMBE, '--');
xlabel('Velocity (m/s)');
legend({'Visc', 'Fric', 'Total'}, 'Location', 'South','fontsize',8)
subplot(244); plot(E.ug, z, E.vgn, z); xlabel('Gas Velocity (m/s)');
legend({'Lateral','Vertical'}, 'Location', 'South','fontsize',8)
subplot(245); plot(E.phi_l, z, E.phi_s, z, E.phi_g, z);
xlabel('Vol Fracs');
legend({'Melt', 'Solid', 'Gas'}, 'Location', 'South','fontsize',8)
subplot(246); plot(E.Chi_hd, z, 100*E.Chi_cd, z);
xlabel('Dissolved Volatiles');
subplot(247); semilogx(E.k_lat, z, 'LineWidth', 2); hold all
xlabel('Perm m^2');
semilogx(E.kvert, z, '--');
xlim([1e-20, 1e0]);
legend({'k_{lat}','k_{mag} model'},...
        'Location', 'South','fontsize',8);
if m.klw.use == 1
    semilogx(m.klw.top./(abs(z*1e-3).^m.klw.mi), z, ':');
    legend({'k_{lat}','k_{mag} model', 'k_{wall} model'},...
        'Location', 'South','fontsize',8);
end
subplot(248); semilogx(E.eta, zMBE); xlabel('Viscosity Pa-s');

PlotTitle(m, ['Profiles at ', num2str(ConvertSecToYear(teval), 2), ' years'])

end

function [phigcVec,plugdepth] = calc_phigc_plugdepth (td, m, plot_opt)
% [phigcVec,plugdepth] = PlotResults('calc_phigc_plugdepth',td, m, 0);

% calculate pressure gradient across faces.
% this has length Nz-1 since we only have (Nz-1) interior faces
Nvz = m.Nz*m.Nv;
p = td.y(m.blk.is.p:m.Nv:Nvz,:)*m.slv.sy(m.blk.is.p);
dz = td.z(2) - td.z(1);
dpdz = 1/dz*(p(2:end,:) - p(1:end-1,:));

phigcVec  = zeros(length(td.x),1);
plugdepth = zeros(length(td.x),1);

for ti = 1:length(td.x)
    yscl = td.y(1:Nvz,ti).*repmat(m.slv.sy(1:m.Nv),length(td.z),1);
    E = tdcFV('calc_exprs',td.x(ti),yscl,td.z',m,dpdz(:,ti),1);
    phigcVec(ti)  = E.zphigc;
    plugdepth(ti) = E.plugdepth;
end

if plot_opt == 1
figure; 
subplot(121); plot(ConvertSecToYear(td.x), td.z(phigcVec)); 
xlabel('Time (yr)'); ylabel('Percolation threshold depth (m)');
subplot(122); plot(ConvertSecToYear(td.x), td.z(plugdepth)); 
xlabel('Time (yr)'); ylabel('Plug depth (m)');

PlotTitle(m, 'Critical depth changes')
end
end

function [] = PlotProfiles (td, m, inds)
% PlotResults('PlotProfiles', td, m, []);

tdvars = extract_y(td, m);
if isempty(inds), inds = round(linspace(1,length(td.x),15)); end
% inds = 1:5:50;

figure; 
set(gcf,'position', [136 569 1384 416], ...
    'DefaultAxesColorOrder', parula(length(inds)+1));  
subplot(141); plot(1e-6*tdvars.p(:,inds), td.z); 
xlabel('Pressure (MPa)');

subplot(142); semilogx(tdvars.v(:,inds), td.z); 
xlabel('Velocity (m/s)');

subplot(143);
plot(100*tdvars.phi_g(:,inds), td.z); hold on;
plot(100*[m.phi_gc, m.phi_gc], ylim, 'k--'); hold off;
xlabel('Porosity (%)');
legend(num2str(ConvertSecToYear(td.x(inds))'),'location','southeast'); 
legend boxoff
    
subplot(144); semilogx(tdvars.mw(:,inds), td.z); 
xlabel('Mole fraction water');

PlotTitle(m, 'Solution profiles at different times')

end

function [] = PlotTitle (m, PlotType)
suptitle({PlotType, ...
    ['V0 = ', num2str(1e-9*m.ch.V0), 'km3, ',...
    'Omega = ', num2str(m.ch.Omega*3600*1e6,4), 'm3/day/MPa, ',...
    'Aspect Ratio = ', num2str(m.ch.AR, 2), ', '], ...
    ['pch0 = ', num2str(m.ch.pdeep*1e-6, 4), 'MPa, ', ...
    'radius = ', num2str(m.R), 'm, ', ...
    'f0 = ', num2str(m.fr.f0), ', ', ...
    'a = ', num2str(m.fr.a), ','], ...
    ['total h2o = ', num2str(m.chi_ch.total.h2o), 'wt%, ', ...
    'total co2 = ', num2str(m.chi_ch.total.co2), 'ppm, ', ...
    'kc = ', num2str(m.k_lat), 'm2, ', ...
    'phigc = ', num2str(m.phi_gc)]});
end












