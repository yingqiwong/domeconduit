function [Trat, Tch, Tasc] = CalcTimescaleRatio (t, p, v, m)

% outputs the timescale ratios from transient solution to determine if we
% are approximating the steady-state solutions.  
% if Trat >> 1, we are close to the steady-state solution.

if any(size(p) == 1)
    p = reshape(p, numel(p), 1);
    v = reshape(v, numel(v), 1);
end

vmean = mean(v,1);
p1 = p(1,:);
v1 = v(1,:);

beta = zeros(size(p1));

for i = 1:length(p1)
    beta(i) = tdcFV('CalcCompressibility', t(i), p1(i), m);
end

Tmag1 = m.ch.V0*beta.*p1./pi/m.R^2./v1;
Tmag2 = m.ch.V0*beta/m.ch.Omega;

Tasc = (m.conduit_length./vmean);

Toutflux = Tmag1'./Tasc';
Tinflux = Tmag2'./Tasc';

Trat = [Toutflux, Tinflux];
Tch = [Tmag1', Tmag2'];
Tasc = Tasc';
end
