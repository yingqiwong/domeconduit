
AddPaths

%% runs steady-state and time-dependent model

use_be   = 0;   % whether you want to use the backward Euler method
plot_opt = 1;   
[ss, td, m, flag] = tdcFV('main',use_be,0);

PlotResults('plot_timeseries', td, m);
PlotResults('PlotProfiles', td, m, []);

%% calculate Mount St. Helens 2004 datasets:
% extruded volume, deformation, gas emissions

VolEx = CalcExtrusionVolume(td,m,1);

%% to only run steady-state model
o = tdcFV('setdef');
[ss, opts, ssflag] = tdcFV('run_ssc', o);