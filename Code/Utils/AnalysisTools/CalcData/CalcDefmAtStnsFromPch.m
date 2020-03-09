function [uOut] = CalcDefmAtStnsFromPch (tyr, dpchIn, m, time, PosRTZ, StnInd)

% material parameters
mu      = m.ch.mu;              % shear modulus
nu      = 0.25;                 % Poisson's ratio
lambda  = 2*mu*nu/(1-2*nu);     % first lame constant
matrl   = [lambda;mu;nu];

% geometry of ellipsoid
plunge  = deg2rad(89.99);       % plunge angle of ellipsoid, deg (90 = prolate)
strike  = deg2rad(0);           % strike of ellipsoid, deg (0 = aligned north)

params = [m.ch.x0; m.ch.y0; m.ch.depth; 0; m.ch.a; m.ch.a*m.ch.AR; plunge; strike];

if ~isempty(tyr)
    dpch = interp1(tyr, dpchIn, time, 'pchip');
else
    dpch = dpchIn;
end

u = zeros(length(time), 3);
for ti = 1:length(time)
    params(4) = dpch(ti);
    [u(ti,1), u(ti,2), u(ti,3)] = ...
        fcn_yangM(params, PosRTZ(StnInd(ti),1), 0, matrl, 0);
end

uOut = 1e3*u; % convert to mm

end
