function [] = PlotDzSols (td, m)

Nz = length(td);

for zi = 1:Nz
    sol(zi).tdv = extract_y(td(zi), m(zi));
    [Vol, Def, CO2] = CalcAllData(td(zi), m(zi), 0);
    sol(zi).tdv.Vol = Vol;
    sol(zi).tdv.Def = Def';
    sol(zi).tdv.CO2 = CO2;
    NzTD(zi) = length(td(zi).z);
end


figure;
set(gcf,'Position',[550 284 1170 516],'defaultlinelinewidth',2,...
    'defaultaxescolororder',parula(Nz+1));

subplot(231); PlotField(sol, 'p', ones(1,Nz), 'Chamber pressure');
hl = legend(num2str(NzTD')); legend boxoff;
title(hl, 'Nz');
subplot(232); PlotField(sol, 'v',       NzTD,       'Exit velocity');
set(gca,'YScale','log');
subplot(233); PlotField(sol, 'phi_g',   NzTD,       'Exit porosity');
subplot(234); PlotField(sol, 'Def',     ones(1,Nz), 'JRO1 deformation');
subplot(235); PlotField(sol, 'Vol',     ones(1,Nz), 'Extruded volume')
subplot(236); PlotField(sol, 'CO2',     ones(1,Nz), 'CO2 emissions')
set(gca,'YScale','log');

end


function [] = PlotField (sol, field, ind, fieldname)

colors = parula(length(sol)+1);
for si = 1:length(sol)
    plot(sol(si).tdv.tyr, sol(si).tdv.(field)(ind(si),:), '-', 'color', colors(si,:));
    hold on;
end
hold off;
% axis tight
title(fieldname);

end


