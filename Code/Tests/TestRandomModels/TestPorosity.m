

o = tdcFV('setdef');
o.R = 30;
o.op = 29e6;
o.phi_gc = 0.2;
o.k_lat = 5e-12;
opts = tdcFV('ss_init',o);
opts.Nz = 801;
[ss, ssflag] = smf_rad_dz('solve', opts);
[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
[td, flag] = tdcFV('des_run', y0, z, m);

% PlotResults('plot_timeseries', td, m);
% PlotResults('PlotProfiles', td, m, []);
% [phigcVec,plugdepth] = PlotResults('calc_phigc_plugdepth', td, m, 1);
CalcExtrusionVolume(td,m,1);
CalcJRO1Def(td,m,1);
CalcExitPorosity(td,m,1);

[tdOut, mOut] = RecalcSolutionAtNewTimes(td, m, ConvertYearToSec(linspace(0,5,51)));
CalcExitPorosity(tdOut,mOut,1);
% PlotResults('PlotProfiles', tdOut, mOut, 1:11);
