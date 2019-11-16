function [op] = CalcOverpressure (td, m)

[p_mag] = CalcMagmastaticPressure(td, m);

op = m.p_ch - p_mag;

end