function [h2o_d, co2_d, h2o_e, co2_e] = CalcSolubilityLiu (p, T, total_h2o, total_co2)
% evaluates the amount of dissolved co2 and h2o according to empirical
% solubility relation in Liu et al. (2005)
%
% INPUTS
% p             Pressure in MPa
% T             temperature in Kelvins
% total_h2o     total h2o in wt%
% total_co2     total co2 in ppm
%
% OUTPUTS
% h2o_d         dissolved h2o in wt%
% co2_d         dissolved co2 in ppm
% h2o_e         exsolved h2o in wt%
% co2_e         exsolved co2 in ppm

% get number of moles of each species and mole fraction
Nmol_h2o = 1e-2*total_h2o/18.02;
Nmol_co2 = 1e-6*total_co2/44.01;
mw = Nmol_h2o./(Nmol_h2o + Nmol_co2);

% partial pressures
Pw = mw.*p;
Pc = (1-mw).*p;

% precompute some stuff for speed
sqrtPw = sqrt(Pw);
Pw15 = Pw.*sqrtPw; % Pw^1.5

% dissolved volatiles
h2o_d = (354.94*sqrtPw + 9.623*Pw - 1.5223*Pw15)/T + ...
    0.0012439*Pw15 + Pc.*(-1.084e-4*sqrtPw - 1.362e-5*Pw);
co2_d = Pc.*(5668 - 55.99*Pw)/T + Pc.*(0.4133*sqrtPw + 2.041e-3*Pw15);

% exsolved volatiles 
h2o_e = total_h2o - h2o_d;
co2_e = total_co2 - co2_d;


end