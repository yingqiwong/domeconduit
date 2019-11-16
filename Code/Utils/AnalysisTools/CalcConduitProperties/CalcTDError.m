function error = CalcTDError (td, m)

% calculates the error from the transient solution, taking LHS - RHS

ErrMat = zeros(size(td.y));

for ti = 1:length(td.x)
     massmat = tdcFV('calc_mass_matrix',td.x(ti),td.y(:,ti),td.z',m);
     dydt = tdcFV('des',td.x(ti),td.y(:,ti),td.z',m);
     ErrMat(:,ti) = massmat*td.yp(:,ti) - dydt;
end

is = m.blk.is;
error.p = ErrMat(is.p:m.Nv:end,:)*m.slv.sy(is.p);
error.v = ErrMat(is.v:m.Nv:end,:)*m.slv.sy(is.v);
error.phi_g = ErrMat(is.phi_g:m.Nv:end,:)*m.slv.sy(is.phi_g);

if m.Nv == 4
    error.mw = ErrMat(is.mw:m.Nv:end,:)*m.slv.sy(is.mw);
else 
    error.mw = zeros(size(error.p));
end

plot_error(error, td, m);
end

function plot_error (error, td, m)

t = ConvertSecToYear(td.x);

figure;
subplot(241); plot(t, 1e-6*error.p(end,:)); title('Pressure at top (MPa)');
subplot(242); plot(t, error.v(end,:)); title('Velocity at top (m/s)');
subplot(243); plot(t, 1e2*error.phi_g(end,:)); title('Porosity at top (%)');
subplot(244); plot(t, error.mw(end,:)); title('Mole fraction water at top');
subplot(245); plot(t, 1e-6*error.p(1,:)); title('Pressure at base (MPa)');
subplot(246); plot(t, error.v(1,:)); title('Velocity at base (m/s)');
subplot(247); plot(t, 1e2*error.phi_g(1,:)); title('Porosity at base (%)');
subplot(248); plot(t, error.mw(1,:)); title('Mole fraction water at base');
suptitle('Error through time');


end