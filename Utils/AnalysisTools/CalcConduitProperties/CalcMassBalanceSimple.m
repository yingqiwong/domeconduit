function [nv, h2o, co2] = CalcMassBalanceSimple (td, m, plot_opt, tquery)
% [nv, h2o, co2] = CalcMassBalanceSimple(td,m,1,[1:100:401]);

tdvars = extract_y(td,m);
is = m.blk.is;

% initialize variables
nv.drho = zeros(m.Nz, 1);
nv.dydt = zeros(m.Nz, 1);
h2o = nv; co2 = nv;

for ti = 1:length(td.x)
    
    % transient terms
    massmat = tdcFV('calc_mass_matrix',td.x(ti),td.y(:,ti),td.z',m);
    LHS = massmat*td.yp(:,ti);
    nv.drho(:,ti)  = LHS(is.p:m.Nv:end);
    h2o.drho(:,ti) = LHS(is.phi_g:m.Nv:end);
    co2.drho(:,ti) = LHS(is.mw:m.Nv:end);
    
    RHS = tdcFV('des',td.x(ti),td.y(:,ti),td.z',m);
    nv.dydt(:,ti)  = RHS(is.p:m.Nv:end);
    h2o.dydt(:,ti) = RHS(is.phi_g:m.Nv:end);
    co2.dydt(:,ti) = RHS(is.mw:m.Nv:end);
end

if plot_opt == 1
    COMPlot(td, nv, h2o, co2, tquery);
end
end

function [] = COMPlot(td, nv, h2o, co2, tquery)


if isempty(tquery)
    tquery = round(linspace(1,length(td.x), 9));
else
    if td.x(tquery(end)) > td.x(end)
        tquery = round(linspace(1,length(td.x), length(tquery)));
    end
end

Nplots = length(tquery);
z = td.z;
dz = td.z(2) - td.z(1);

figure;
PlotInd = 0;
for tqi = tquery
    
    PlotInd = PlotInd + 1;

    subplot(3,Nplots,PlotInd);
    plot(nv.drho(:,tqi), z); hold on;
    plot(nv.dydt(:,tqi), z,'--'); hold off;
    title(sprintf('t = %.1f year', ConvertSecToYear(td.x(tqi))));
    
    subplot(3,Nplots,Nplots+PlotInd);
    plot(h2o.drho(:,tqi), z); hold on;
    plot(h2o.dydt(:,tqi), z,'--'); hold off;
    
    subplot(3,Nplots,2*Nplots+PlotInd);
    plot(co2.drho(:,tqi), z); hold on;
    plot(co2.dydt(:,tqi), z,'--'); hold off;
    
end


subplot(3,Nplots,PlotInd); 
legend('LHS', 'RHS','location','southeast'); legend boxoff;
subplot(3,Nplots,Nplots+PlotInd);
legend('LHS', 'RHS','location','southeast'); legend boxoff;
subplot(3,Nplots,2*Nplots+PlotInd);
legend('LHS', 'RHS','location','southeast'); legend boxoff;

    
end














