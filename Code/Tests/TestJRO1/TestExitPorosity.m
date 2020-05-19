
addpath('../../../../ModelTestFiles/');
load JRO1TrialError.mat
set(0,'defaultlinelinewidth',2);

%%
isol = 4;
ss = sols(isol).ss;
td = sols(isol).td;
m = sols(isol).m;
Porosity = CalcExitPorosity(td,m,1);