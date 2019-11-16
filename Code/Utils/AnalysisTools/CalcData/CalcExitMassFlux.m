function [Q_kgs, v, rho] = CalcExitMassFlux (td, m)
% calculates the exit mass flux of the magma in kg/s

is = m.blk.is;

ExitCond = td.y(end-3:end,:).*m.slv.sy;
p = ExitCond(is.p, :);
v = ExitCond(is.v, :);
phi_g = ExitCond(is.phi_g, :);
mw = ExitCond(is.mw, :);

rho_g = p.*(mw./(m.Rw*m.T) + (1 - mw)./(m.Rc*m.T)); 
rho = rho_g.*phi_g + m.rho_s*(1-phi_g);

Q_kgs = rho.*pi*m.R^2.*v;


end