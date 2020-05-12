
AddPaths

set(0,'defaultlinelinewidth',2, 'defaultaxesfontsize', 14);

%% to only run steady-state model

o = tdcFV('setdef');
[ss, opts, ssflag] = tdcFV('run_ssc', o);

%% runs steady-state and time-dependent model

use_be   = 0;   % whether you want to use the backward Euler method
plot_opt = 0;   
[ss, td, m, flag] = tdcFV('main', use_be, plot_opt);

%% calculate Mount St. Helens 2004 datasets:
% extruded volume, deformation, gas emissions

VolEx = CalcExtrusionVolume(td,m,1);
JRO1Def = CalcJRO1Def(td,m,1);
[~, CO2Flux, ~, ~, ~] = CalcGasEmissions(td,m,1);

%% some additional plotting functions

PlotResults('plot_timeseries', td, m);
PlotResults('PlotProfiles', td, m, []);
PlaytdSol(td, m, 26, 'tdNomModel')

