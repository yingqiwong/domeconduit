function [p_mag] = CalcMagmastaticPressure (td, m)

Rho = CalcConduitProperties(td, m, {'rho'}, 0, []);

dz = td.z(2) - td.z(1);
p_mag = m.p_top + dz*trapz(Rho.rho*m.g);

end