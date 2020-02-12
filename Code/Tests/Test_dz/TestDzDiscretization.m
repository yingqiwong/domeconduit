
AddPaths;

%%

NzSS = [301, 401, 501, 601, 701, 801, 901, 1001];
NzTD = (NzSS-1)*0.5 + 1;
Nz = length(NzSS);

o = tdcFV('setdef');
opts = tdcFV('ss_init',o);

RunTime = zeros(size(NzSS));

for zi = 1:Nz
    fprintf('Nz = %d.\n', NzSS(zi)); 
    
    opts.Nz = NzSS(zi);
    
    [ss, ssflag] = smf_rad_dz('solve', opts);
    [m, y0, z] = tdcFV('td_init', ss.m, ss, 0, 0);
    
    tic;
    [td, m, flag] = tdcFV('run_tdc', y0, z, m);
    RunTime(zi) = toc;
    
    sol(zi).ss = ss;
    sol(zi).m  = m;
    sol(zi).td = td;
end

%% add data

for zi = 1:Nz
    sol(zi).tdv = extract_y(sol(zi).td, sol(zi).m);
    [Vol, Def, CO2] = CalcAllData(sol(zi).td, sol(zi).m, 0);
    sol(zi).tdv.Vol = Vol;
    sol(zi).tdv.Def = Def';
    sol(zi).tdv.CO2 = CO2;
end

%% plot results

figure;
set(gcf,'Position',[550 284 1170 516],'defaultlinelinewidth',2,...
    'defaultaxescolororder',parula(Nz+1));

subplot(231); PlotField(sol, 'p', ones(1,Nz), 'Chamber pressure');
hl = legend(num2str(NzTD')); legend boxoff;
title(hl, 'Nz');
subplot(232); PlotField(sol, 'v',       NzTD,       'Exit velocity');
subplot(233); PlotField(sol, 'phi_g',   NzTD,       'Exit porosity');
subplot(234); PlotField(sol, 'Def',     ones(1,Nz), 'JRO1 deformation');
subplot(235); PlotField(sol, 'Vol',     ones(1,Nz), 'Extruded volume')
subplot(236); PlotField(sol, 'CO2',     ones(1,Nz), 'CO2 emissions')

%%

VolErr = zeros(Nz,1);
DefErr = zeros(Nz,1);
phigErr= zeros(Nz,1);
for zi = 1:Nz
    VolErr(zi) = abs((sol(zi).tdv.Vol(end) - sol(end).tdv.Vol(end))/sol(end).tdv.Vol(end));
    DefErr(zi) = abs((sol(zi).tdv.Def(end) - sol(end).tdv.Def(end))/sol(end).tdv.Def(end));
    phigErr(zi)= abs((sol(zi).tdv.phi_g(end) - sol(end).tdv.phi_g(end))/sol(end).tdv.phi_g(end));
end

figure; 
subplot(221); plot(NzTD, RunTime, '+-'); title('Run Time (s)');
subplot(222); plot(NzTD, 100*VolErr, '+-'); title('total extruded vol err%');
subplot(223); plot(NzTD, 100*DefErr, '+-'); title('total JRO1 def err%');
subplot(224); plot(NzTD, 100*phigErr, '+-'); title('final porosity err%');

%%

function [] = PlotField (sol, field, ind, fieldname)

colors = parula(length(sol)+1);
for si = 1:length(sol)
    plot(sol(si).tdv.tyr, sol(si).tdv.(field)(ind(si),:), '-', 'color', colors(si,:));
    hold on;
end
hold off;
axis tight
title(fieldname);

end


