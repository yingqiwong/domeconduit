function [op] = CalcOPfromPch (m)
% calculates the excess pressure op in Pa from the solution options
% structure m

op = m.p_ch - m.rho_l*m.g*m.conduit_length - m.p_top;

end