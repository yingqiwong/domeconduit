function [o] = GetOFromOpts (opts)
% retrieves tdcFV options o from the structure opts (definitions as in tdcFV)

o = tdcFV('setdef');

o.V0        = opts.ch.V0;
o.Omega     = opts.ch.Omega;
o.AR        = opts.ch.AR;
o.L         = opts.conduit_length;
o.op        = opts.p_ch - opts.p_top - 2200*9.81*opts.conduit_length;
o.R         = opts.R;
o.total_h2o = opts.total_h2o;
o.total_co2 = opts.total_co2;
o.phi_gc    = opts.phi_gc;
o.k_lat     = opts.k_lat;
o.f0        = opts.fr.f0;
o.a         = opts.fr.a;

end