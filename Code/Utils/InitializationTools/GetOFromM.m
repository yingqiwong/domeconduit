function [o] = GetOFromM (m)
% retrieves tdcFV options o from the structure m (definitions as in tdcFV)

o = tdcFV('setdef');

o.V0        = m.ch.V0;
o.Omega     = m.ch.Omega;
o.AR        = m.ch.AR;
o.L         = m.conduit_length;
o.op        = m.p_ch - m.p_top - 2200*9.81*m.conduit_length;
o.R         = m.R;
o.total_h2o = m.chi_ch.total.h2o;
o.total_co2 = m.chi_ch.total.co2;
o.phi_gc    = m.phi_gc;
o.k_lat     = m.k_lat;
o.f0        = m.fr.f0;
o.a         = m.fr.a;

end