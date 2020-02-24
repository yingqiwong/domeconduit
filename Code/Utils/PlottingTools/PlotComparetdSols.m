function [] = PlotComparetdSols (td, m, lgdtxt, Nsols)

if nargin<4
    Nsols = min([length(td), 4]);
    fprintf('Plotting first %d solutions only.\n', Nsols);
end

if nargin<3
    lgdtxt = repmat({''}, Nsols, 1);
end

lineColors = parula(Nsols+1);

figure;
set(gcf, 'position', [78 406 1355 505], ...
    'DefaultLineLineWidth',2, 'DefaultLineMarkerSize',3)

for isol = 1:Nsols
    tdv = extract_y(td(isol), m(isol));
    
    subplot(241); plot(tdv.tyr, tdv.p(end,:),   'color', lineColors(isol,:));   hold on;
    subplot(242); semilogy(tdv.tyr, tdv.v(end,:),   'color', lineColors(isol,:));   hold on;
    subplot(243); plot(tdv.tyr, tdv.phi_g(end,:), 'color', lineColors(isol,:)); hold on;
    subplot(244); plot(tdv.tyr, tdv.mw(end,:),  'color', lineColors(isol,:));   hold on;
    subplot(245); plot(tdv.tyr, tdv.p(1,:),     'color', lineColors(isol,:));   hold on;
    subplot(246); semilogy(tdv.tyr, tdv.v(1,:),     'color', lineColors(isol,:));   hold on;
    subplot(247); plot(tdv.tyr, tdv.phi_g(1,:), 'color', lineColors(isol,:));   hold on;
    subplot(248); plot(tdv.tyr, tdv.mw(1,:),    'color', lineColors(isol,:));   hold on;

end


subplot(241); title('(a) Pressure at top (MPa)'); xlabel('Time (year)');
subplot(242); title('(b) Velocity at top (m/s)'); xlabel('Time (year)');
subplot(243); title('(c) Porosity at top (%)'); xlabel('Time (year)');
subplot(244); title('(d) Mole fraction water at top'); xlabel('Time (year)');
subplot(245); title('(e) Pressure at base (MPa)'); xlabel('Time (year)');
subplot(246); title('(f) Velocity at base (m/s)'); xlabel('Time (year)');
subplot(247); title('(g) Porosity at base (%)'); xlabel('Time (year)');
subplot(248); title('(h) Mole fraction water at base'); xlabel('Time (year)');
legend(lgdtxt, 'location', 'best', 'box', 'off');



figure;
set(gcf, 'position', [78 406 1355 505], ...
    'DefaultLineLineWidth',2, 'DefaultLineMarkerSize',3)
for isol = 1:Nsols
    tdv = extract_y(td(isol), m(isol));
    
    subplot(141); plot(tdv.p(:,end), tdv.z, 'color', lineColors(isol,:));   hold on;
    subplot(142); plot(tdv.v(:,end), tdv.z, 'color', lineColors(isol,:));   hold on;
    subplot(143); plot(tdv.phi_g(:,end), tdv.z, 'color', lineColors(isol,:));   hold on;
    subplot(144); plot(tdv.mw(:,end), tdv.z, 'color', lineColors(isol,:));   hold on;
end

subplot(141); title('(a) Pressure (MPa)'); ylabel('Depth (m)');
subplot(142); title('(b) Velocity (m/s)'); 
subplot(143); title('(c) Porosity (%)'); 
subplot(144); title('(d) Mole fraction water'); 
% suptitle('Results at t = 3.3 years');


% calculate predicted data
CalcExtrusionVolume(td(1), m(1), 1); hold on;
for isol = 2:Nsols
    tyr = ConvertSecToYear(td(isol).x);
%     Vol = CalcExtrusionVolume(td(isol), m(isol), 0);
    Vol = td(isol).y(end,:);
    plot(tyr, Vol, 'linewidth', 2, 'color', lineColors(isol,:)); 
end

CalcJRO1Def(td(1), m(1), 1); hold on
for isol = 2:Nsols
    tyr = ConvertSecToYear(td(isol).x);
    JRO1 = CalcJRO1Def(td(isol), m(isol), 0);
    plot(tyr, 1e3*JRO1, 'linewidth', 2, 'color', lineColors(isol,:));
end

CalcGasEmissions(td(1), m(1), 1); hold on;
for isol = 2:Nsols
    tyr = ConvertSecToYear(td(isol).x);
    [~,CO2] = CalcGasEmissions(td(isol), m(isol), 0);
    plot(tyr, CO2, 'linewidth', 2, 'color', lineColors(isol,:));
end










end