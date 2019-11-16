function [COM, nv, h2o, co2] = CalcMassBalance (td, m, plot_opt, tquery)
% [COM, nv, h2o, co2] = CalcMassBalance(td,m,1,[1:100:401]);

tdvars = extract_y(td,m);
is = m.blk.is;
dz = td.z(2) - td.z(1);

vn = tdvars.v(2:end,:);
vs = tdvars.v(1:end-1,:);

% initialize variables
nv.mass.drho = zeros(size(vn));
nv.mass.c   = zeros(size(vn));  % cell center
nv.mass.in  = zeros(size(vn));  % mass flux coming in
nv.mass.out = zeros(size(vn));  % mass flux going out

h2o.mass     = nv.mass;
h2o.vert.in  = zeros(size(vn));
h2o.vert.out = zeros(size(vn));
h2o.lat.c = zeros(size(vn));

co2 = h2o;

for ti = 1:length(td.x)
    
    % transient terms
    massmat = tdcFV('calc_mass_matrix',td.x(ti),td.y(:,ti),td.z',m);
    LHS = massmat*td.yp(:,ti);
    nv.mass.drho(:,ti)  = LHS(m.Nv+is.p:m.Nv:end);
    h2o.mass.drho(:,ti) = LHS(m.Nv+is.phi_g:m.Nv:end);
    co2.mass.drho(:,ti) = LHS(m.Nv+is.mw:m.Nv:end);
    
    % calculate dependent variables
    E = CalcFieldsAtOneTimeStep(td, td.x(ti), m);
    
    % non-volatile components
    nv.mass.c(:,ti) = E.nv.mass(2:end);
    [nvn, nvs]  = tdcFV('get_FaceVals', E.nv.mass, m.dQUICKn);
    nv.mass.in(:,ti)  = nvs(2:end).*vs(:,ti);
    nv.mass.out(:,ti) = nvn(2:end).*vn(:,ti);
    
    % water
    watermass = E.h2o.liq + E.h2o.gas;
    h2o.mass.c(:,ti) = watermass(2:end);
    [h2on, h2os] = tdcFV('get_FaceVals', watermass, m.dQUICKn);
    h2o.mass.in(:,ti)  = h2os(2:end).*vs(:,ti);
    h2o.mass.out(:,ti) = h2on(2:end).*vn(:,ti);
    
    [h2ovn, h2ovs] = tdcFV('get_FaceVals', E.h2o.gas, m.dQUICKn);
    h2o.vert.in(:,ti)  = h2ovs(2:end).*E.vgs(2:end);
    h2o.vert.out(:,ti) = h2ovn(2:end).*E.vgn(2:end);
    
    h2o.lat.c(:,ti) = E.h2o.lat(2:end);
    
    % carbon dioxide
    co2mass = E.co2.liq + E.co2.gas;
    co2.mass.c(:,ti) = co2mass(2:end);
    [co2n, co2s] = tdcFV('get_FaceVals', co2mass, m.dQUICKn);
    co2.mass.in(:,ti)  = co2s(2:end).*vs(:,ti);
    co2.mass.out(:,ti) = co2n(2:end).*vn(:,ti);
    
    [co2vn, co2vs] = tdcFV('get_FaceVals', E.co2.gas, m.dQUICKn);
    co2.vert.in(:,ti)  = co2vs(2:end).*E.vgs(2:end);
    co2.vert.out(:,ti) = co2vn(2:end).*E.vgn(2:end);
    
    co2.lat.c(:,ti) = E.co2.lat(2:end);
end

COM.nv  = nv.mass.drho*dz + nv.mass.out - nv.mass.in;
COM.h2o = h2o.mass.drho*dz...
                + h2o.mass.out - h2o.mass.in ...
                + h2o.vert.out - h2o.vert.in ...
                + h2o.lat.c*dz;
COM.co2 = co2.mass.drho*dz...
                + co2.mass.out - co2.mass.in...
                + co2.vert.out - co2.vert.in...
                + co2.lat.c*dz;

if plot_opt == 1
    COMPlot(td, nv, h2o, co2, COM, tquery);
end
end

function [] = COMPlot(td, nv, h2o, co2, COM, tquery)

if isempty(tquery)
    tquery = round(linspace(1,length(td.x), 5));
else
    if td.x(tquery(end)) > td.x(end)
        tquery = round(linspace(1,length(td.x), length(tquery)));
    end
end

Nplots = length(tquery);
z = td.z(2:end);
dz = td.z(2) - td.z(1);

figure;
PlotInd = 0;
for tqi = tquery
    
    PlotInd = PlotInd + 1;

    subplot(3,Nplots,PlotInd);
    plot(nv.mass.drho(:,tqi)*dz, z); hold on;
    plot(nv.mass.in(:,tqi) - nv.mass.out(:,tqi), z,'--'); hold off;
    title(sprintf('t = %.1f year', ConvertSecToYear(td.x(tqi))));
    
    subplot(3,Nplots,Nplots+PlotInd);
    plot(h2o.mass.drho(:,tqi)*dz, z); hold on;
    plot(h2o.mass.in(:,tqi) - h2o.mass.out(:,tqi), z); 
    plot(h2o.vert.in(:,tqi) - h2o.vert.out(:,tqi), z); 
    plot(h2o.lat.c(:,tqi)*dz, z); hold off;
    
    subplot(3,Nplots,2*Nplots+PlotInd);
    plot(co2.mass.drho(:,tqi)*dz, z); hold on;
    plot(co2.mass.in(:,tqi) - co2.mass.out(:,tqi), z); 
    plot(co2.vert.in(:,tqi) - co2.vert.out(:,tqi), z); 
    plot(co2.lat.c(:,tqi)*dz, z); hold off;
    
end


subplot(3,Nplots,PlotInd); 
legend('Transient', 'Mass flux','location','southeast'); legend boxoff;
subplot(3,Nplots,Nplots+PlotInd);
legend('Transient', 'Mass flux', 'Vert gas', 'Lat gas','location','southeast'); 
legend boxoff;
subplot(3,Nplots,2*Nplots+PlotInd);
legend('Transient', 'Mass flux', 'Vert gas', 'Lat gas','location','southeast'); 
legend boxoff;
    
end














