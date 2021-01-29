% this file contains various commands that you can use to run the
% steady-state models

AddPaths

set(0,'defaultlinelinewidth',2, 'defaultaxesfontsize', 14);

%% to only run steady-state model

o = tdcFV('setdef');
[ss, opts, ssflag] = tdcFV('run_ssc', o);

%% if mex file doesn't work (depends on os)

o                = tdcFV('setdef');
opts             = tdcFV('ss_init',o);
opts.slv.use_mex = 0;

[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);

