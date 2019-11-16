function [MBE] = CalcMBEParts (td, m, plot_opt, tquery)
%[MBE] = CalcMBEParts (td, m, 1, [1:100:801]);

Nt = length(td.x);

MBE.dpdz = zeros(m.Nz-1, Nt);
MBE.rhon = zeros(m.Nz-1, Nt);
MBE.tauR = zeros(m.Nz-1, Nt);
MBE.etan = zeros(m.Nz-1, Nt);
MBE.sigc = zeros(m.Nz-1, Nt);

for ti = 1:length(td.x)
    E = PlotResults('AddFields', td, td.x(ti), m);
    
    [rhon, ~] = tdcFV('get_FaceVals', E.rho, m.dQUICKn);
    
    MBE.dpdz(:,ti) = E.dpdz;
    MBE.rhon(:,ti) = rhon(1:end-1);
    MBE.tauR(:,ti) = E.tauR;
    MBE.etan(:,ti) = E.eta;
    MBE.sigc(:,ti) = m.sig.slope*abs(E.zMBE) + m.sig.offset; 
end

MBE.zMBE = E.zMBE;

if plot_opt == 1
    MBEPlot(MBE, td, m, tquery)
end

end

function [] = MBEPlot (MBE, td, m, tquery)

if isempty(tquery)
    tquery = round(linspace(1,length(td.x), 5));
else
    if td.x(tquery(end)) > td.x(end)
        tquery = round(linspace(1,length(td.x), length(tquery)));
    end
end


Nplots = length(tquery);
z = MBE.zMBE;

figure;
PlotInd = 0;
for tqi = tquery
    PlotInd = PlotInd + 1;
    
    subplot(1,Nplots,PlotInd);
    plot(-MBE.dpdz(:,tqi), z, MBE.rhon(:,tqi)*m.g, z);
end

subplot(1,Nplots,PlotInd);
legend('dpdz', 'rho*g', 'location', 'southeast'); legend boxoff;

end











