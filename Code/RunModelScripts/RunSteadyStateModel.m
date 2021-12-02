% this file contains various commands that you can use to run the
% steady-state models

AddPaths

set(0,'defaultlinelinewidth',2, 'defaultaxesfontsize', 14);

%% to only run steady-state model

o = tdcFV('setdef');
o.op = 15e6;
[ss, opts, ssflag] = tdcFV('run_ssc', o);
ss = smf_rad_dz('add_fields',ss);

%% if mex file doesn't work (depends on os)

o                = tdcFV('setdef');
opts             = tdcFV('ss_init',o);
opts.slv.use_mex = 0;

[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);
ss = smf_rad_dz('add_fields',ss);

%% turn off plug gas loss

o    = tdcFV('setdef');
% o.op = 18e6;
o.phi_gc = 0;
opts = tdcFV('ss_init',o);
% opts.plug_gas_loss = 0;

[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);
ss = smf_rad_dz('add_fields',ss);
