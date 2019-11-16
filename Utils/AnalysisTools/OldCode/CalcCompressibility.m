
tdvars = extract_y(td, m);
for i = 1:length(td.x)
    beta(i) = tdcFV('CalcCompressibility', tdvars.p(1,i), m);
end

figure; plot(td.x, beta)

%%

E = PlotResults('AddFields', td, td.x(1), m);

pch = m.ch.pdeep;
p1 = pch + E.rho(1)*m.g*m.ch.a;
p2 = pch + 2300*m.g*m.ch.a;
mu = 20e9;

pcc = fzero(@(pcc) ChamberCenterPressure(pcc, pch, m), [pch, 300e6]);
rhocc = CalcRho(pcc, m);

pvec = pcc + 1e6*linspace(-5,5,51);
for i = 1:length(pvec)
    [rhovec(i)] = CalcRho(pvec(i), m);
end
% figure; plot(pvec, rhovec);

pdiff = pcc + [-1e6, 1e6];
rhodiff(1) = CalcRho(pdiff(1), m);
rhodiff(2) = CalcRho(pdiff(2), m);

% beta_mag = 1/rhocc*mean(diff(rhovec)./diff(pvec))
beta_mag = 1/rhocc*(diff(rhodiff)/diff(pdiff))
beta_chb = 3/4/mu

beta_sys = (beta_mag + beta_chb)

function [r] = ChamberCenterPressure (pcc, pch, m)

[rho] = CalcRho(pcc, m);
r = pcc - pch - rho*m.g*m.ch.a;

end

function [rho] = CalcRho (p, m)

[mw, phi_g] = tdcFV('calc_ch_gas', p, m);

% Volatile mass fractions. 
[Chi_hd, Chi_cd] = tdcFV('solubility_liu_explicit', p, m.T, mw);
Chi_hd = 1e-2*Chi_hd; 
Chi_cd = 1e-6*Chi_cd;
c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
c12 = c1.*c2;

% solid fraction
[chi_s] = tdcFV('calc_chi_s', m.pf, p*1e-6);
% chi_s = m.pf.chi(end); % hold chi_s fixed?
phi_s = chi_s.*m.rho_l.*c12.*(1 - phi_g)./...
    (m.rho_s + chi_s.*(c12*m.rho_l - m.rho_s));

% liquid fraction
phi_l = (1 - phi_s - phi_g).*c2; 

% gas density 
rho_g = p.*(mw./(m.Rw*m.T) + (1 - mw)./(m.Rc*m.T)); 

% bulk density
rho = m.rho_l*phi_l.*c1 + m.rho_s.*phi_s + rho_g.*phi_g; 

end