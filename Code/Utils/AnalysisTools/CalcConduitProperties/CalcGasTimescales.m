function [Tg, zphigc] = CalcGasTimescales (td, m)

tdvars = extract_y(td,m);

vars = {'ug', 'dpdz', 'kvert','zphigc'};
GasVel = CalcConduitProperties (td, m, vars, 0, []);
zphigc = GasVel.zphigc;

dpdz = [GasVel.dpdz(1,:); GasVel.dpdz];

vgdiff = -GasVel.kvert/m.eta_g.*dpdz;

Magma   = m.conduit_length./mean(tdvars.v);
LatGas  = m.R./GasVel.ug;
VertGas = abs(td.z(GasVel.zphigc))./vgdiff;

Tg = zeros(3,length(td.x));
% row 1: magma timescale
% row 2: lateral gas loss timescale
% row 3: vertical gas loss timescale

for ti = 1:length(td.x)
    Tg(1,ti) = Magma(ti);
    Tg(2,ti) = mean(LatGas(GasVel.zphigc(ti):m.PlugDepth,ti));
    Tg(3,ti) = mean(VertGas(GasVel.zphigc(ti):m.PlugDepth,ti));
end

end