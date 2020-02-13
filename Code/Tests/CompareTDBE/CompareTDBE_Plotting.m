function varargout = CompareTDBE_Plotting (varargin)
[varargout{1:nargout}] = feval(varargin{:});
end

function [td, be] = CalcBEData (td, m, be, mbedef)

[td.Vol, td.Def, td.CO2] = CalcData(td, m);
td.tdv = extract_y(td,m);

for mi = 1:numel(be)
    
    if isempty(be(mi).x), continue; end
        
    [be(mi).Vol, be(mi).Def, be(mi).CO2] = CalcData(be(mi), mbedef);
    be(mi).bev = extract_y(be(mi),mbedef);
end
end

function [vol, def, co2] = CalcData (td, m)
vol = td.y(end,:);

[ur, ut, uz] = CalcJRO1Def(td, m, 0); 
def = 1e3*ur;

[~,co2] = CalcGasEmissions(td,m,0);
end

function [] = PlotCompareTDBE (td, m, be, mbedef, dpVec)

[td, be] = CalcBEData (td, m, be, mbedef);

beplt = [];
for ibe = 1:length(be)
    if ~isempty(be(ibe).x), beplt = [beplt, ibe]; end
end

figure;
set(gcf,'Position',[550 284 1170 516],'defaultlinelinewidth',2,...
    'defaultaxescolororder',parula(length(beplt)+2));

subplot(231);
plot(td.tdv.tyr, td.tdv.p(1,:), '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).bev.p(1,:), '+-');
end
hold off;
legend([{'td'}; strcat({'dp = '}, num2str(dpVec(beplt)'*1e-6), {' MPa'})]); 
legend boxoff;
title('Chamber pressure');

subplot(232);
plot(td.tdv.tyr, td.tdv.v(end,:), '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).bev.v(end,:), '+-');
end
hold off;
title('Exit velocity');

subplot(233);
plot(td.tdv.tyr, td.tdv.phi_g(end,:), '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).bev.phi_g(end,:), '+-');
end
hold off;
title('Exit porosity');

subplot(234);
plot(td.tdv.tyr, td.Def, '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).Def, '+-');
end
hold off;
title('JRO1 deformation');

subplot(235);
plot(td.tdv.tyr, td.Vol, '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).Vol, '+-');
end
hold off;
ylim([0,1.2*td.Vol(end)])
title('Extruded volume');

subplot(236);
plot(td.tdv.tyr, td.CO2, '+-'); hold on;
for ibe = beplt
    plot(be(ibe).bev.tyr, be(ibe).CO2, '+-');
end
% set(gca,'yscale','log');
hold off;
title('CO2 emissions');



end



