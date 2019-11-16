function [SignChangeTime, Trat] = CalcPchParts (td, m, plot_opt)

is = m.blk.is;
pch = td.y(is.p,:)*m.slv.sy(is.p);
vch = td.y(is.v,:)*m.slv.sy(is.v);

tdvars = extract_y(td,m);

for ti = 1:length(td.x)
    Trat(ti,:) = CalcTimescaleRatio(tdvars.p(:,ti), tdvars.v(:,ti), m);
end

% find where dpchdt changes sign (== peaks and troughs of pch)
dpchdt = td.yp(is.p,:)*m.slv.sy(is.p);
temp = diff(sign(dpchdt')); % detect sign change
SignChangeInd = find(temp~=0);
SignChangeTime = td.x(SignChangeInd);

% parts of the pch evolution equation
Term1 = m.ch.Omega*(m.ch.pdeep-pch)/m.ch.V0/m.ch.beta;
Term2 = pi*m.R^2*vch/m.ch.V0/m.ch.beta;

if plot_opt == 1 
% plotting
figure; 
subplot(411); plot(td.x, 1e-6*pch); 
hold on; plot_vline(td.x(SignChangeInd)); hold off;
ylabel('Chamber Pressure (MPa)'); 
subplot(412); plot(td.x, vch); 
hold on; plot_vline(td.x(SignChangeInd)); hold off;
ylabel('Chamber velocity (m/s)'); 
subplot(413); plot(td.x, Term1, td.x, Term2); 
hold on; plot_vline(td.x(SignChangeInd)); hold off;
ylabel('dpch/dt (Pa/s)'); 
legend('Influx', 'Outflux'); legend boxoff
subplot(414); plot(td.x, 1e-6*dpchdt);
hold on; plot_vline(td.x(SignChangeInd)); hold off;
ylabel('dpchdt (MPa/s)'); xlabel('Time (s)'); grid on; 
% ylim([-0.2, 0.2]*1e-5);

suptitle(sprintf(...
    'Chamber pressure evolution components\n pch = %d MPa, kc = %.1e m2, Omega = %.1e m3/day/MPa, beta = %.1e /Pa',...
    m.ch.pdeep*1e-6, m.k_lat, m.ch.Omega*(24*3600)*1e6, m.ch.beta));


figure;
subplot(211); plot(td.x, Trat(:,1)); title('pch timescale 1');
hold on; plot_vline(td.x(SignChangeInd)); hold off;  xlim([0,3e7]);
subplot(212); plot(td.x, Trat(:,2)); title('pch timescale 2');
hold on; plot_vline(td.x(SignChangeInd)); hold off; xlim([0,3e7]);


suptitle(sprintf(...
    'Chamber pressure evolution components\n pch = %d MPa, kc = %.1e m2, Omega = %.1e m3/day/MPa, beta = %.1e /Pa',...
    m.ch.pdeep*1e-6, m.k_lat, m.ch.Omega*(24*3600)*1e6, m.ch.beta));

end

end

function [] = TheoretTc ()

E1 = PlotResults('AddFields', td, 5e7, m);
E2 = PlotResults('AddFields', td, 5.2e7, m);
detadt = (E2.eta(1) - E1.eta(1))/(td.x(2) - td.x(1));

omeg2 = -pi*m.R^4/8/m.ch.V0/m.ch.beta/E1.eta(1)^2/m.conduit_length*detadt;
f = sqrt(omeg2)/2/pi;

zeta = m.ch.Omega/m.ch.V0/m.ch.beta + ...
    pi*m.R^4/8/E1.eta(1)/m.conduit_length/m.ch.V0/m.ch.beta;


end


function [] = plot_vline (xvals)

for i = 1:length(xvals)
    plot([xvals(i), xvals(i)], ylim, 'k:');
end

end