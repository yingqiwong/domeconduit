
clear all;
AddPaths

%%

o = tdcFV('setdef');
opts = tdcFV('ss_init',o);
[ss, ssflag] = smf_rad_dz('solve', opts);
[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
tic;
[td, m, flag] = tdcFV('run_tdc', y0, z, m);
tdRT = toc;
PlotResults('plot_timeseries',td,m);

%%
tic; 
[be, flag] = RunBE(y0, z, m, 1); 
beRT=toc;
CompareTDBE_Plotting('PlotCompareTDBE', td, m, be, m, 1e6);

%%

o = tdcFV('setdef');
dpVec = [0.1e6, 0.25e6, 0.5e6, 1e6, 2e6];
[td, m, tdRunTime, be, mbe, beRunTime] = CompareTDBE_dpVec('main', o, dpVec);
CompareTDBE_Plotting('PlotCompareTDBE', td, m, be, mbe, dpVec)

