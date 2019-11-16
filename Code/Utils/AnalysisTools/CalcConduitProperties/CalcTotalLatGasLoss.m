function [LatFlux] = CalcTotalLatGasLoss (td, m, plot_opt)

dz = abs(td.z(2) - td.z(1));
LatFlux = zeros(1,length(td.x));

for ti = 1:length(td.x)
    E           = PlotResults('AddFields', td, td.x(ti), m);
    LatFluxdz   = 2*pi*m.R*E.ug;
    LatFlux(ti) = dz*trapz(LatFluxdz);
end

if plot_opt == 1
    figure;
    plot(ConvertSecToYear(td.x), LatFlux);
    xlabel('Time (year)'); ylabel('Lateral Gas Flux (m3/s)');
end

end