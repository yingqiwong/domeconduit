% this file contains various commands that you can use to run the
% time-dependent model

AddPaths

set(0,'defaultlinelinewidth',2, 'defaultaxesfontsize', 14);

%% runs steady-state and time-dependent model

[ss, td, m, flag] = tdcFV('main');

%% more verbose way of running the previous section 
% this allows you more control of solver options

% load options for steady-state model
o    = tdcFV('setdef');
opts = tdcFV('ss_init',o);

% run steady-state model
[ss, opts, ssflag] = tdcFV('run_ssc_opts', opts);

% load time-dependent options
[m, y0, z]    = tdcFV('td_init', ss.m, ss, 0, 0);
[td, m, flag] = tdcFV('run_tdc', y0, z, m);

%% calculate Mount St. Helens 2004 datasets:
% extruded volume, deformation, gas emissions

VolEx        = CalcExtrusionVolume(td,m,1);
JRO1Def      = CalcJRO1Def(td,m,1);
[~, CO2Flux] = CalcGasEmissions(td,m,1);

%% some plotting functions

PlotResults('plot_timeseries', td, m);
PlotResults('PlotProfiles', td, m, []);
PlaytdSol(td, m, 26, 'tdNomModel')

