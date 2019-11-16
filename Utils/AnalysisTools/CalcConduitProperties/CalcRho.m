function [rho] = CalcRho (p, phi_g, mw, m)

% Volatile mass fractions. 
[Chi_hd, Chi_cd] = tdcFV('solubility_liu_explicit', p, m.T, mw);
Chi_hd = 1e-2*Chi_hd; 
Chi_cd = 1e-6*Chi_cd;
c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
c12 = c1.*c2;

% solid fraction
[chi_s] = tdcFV('calc_chi_s', m.pf, p*1e-6);

phi_s = chi_s.*m.rho_l.*c12.*(1 - phi_g)./...
    (m.rho_s + chi_s.*(c12*m.rho_l - m.rho_s));

% liquid fraction
phi_l = (1 - phi_s - phi_g).*c2; 
E.phi_l = phi_l;

% gas density 
rho_g = p.*(mw./(m.Rw*m.T) + (1 - mw)./(m.Rc*m.T)); 

% bulk density
rho = m.rho_l*phi_l.*c1 + m.rho_s.*phi_s + rho_g.*phi_g; 
end