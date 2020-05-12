
addpath('../../../../ModelTestFiles/');
load JRO1TrialError.mat
set(0,'defaultlinelinewidth',2);

%%

isol = 1;
o = GetOFromM(sols(isol).m);
[ss, opts, ssflag] = tdcFV('run_ssc', o);
[m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
m.tdep.p_bot_tvary = 0;
[td, m, flag] = tdcFV('run_tdc', y0, z, m);
PlotResults('plot_timeseries', td, m);
PlotResults('PlotProfiles', td, m, []);

%%

PlotResults('plot_timeseries',sols(isol).td,sols(isol).m);

tdv = extract_y(td,m);
figure;
plot(ss.y(3,:), -ss.z)
hold on;
plot(tdv.phi_g(:,end), tdv.z);
