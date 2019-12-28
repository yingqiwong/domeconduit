function [be, flag] = RunBE (y0, z, m, dtvary)
% runs the conduit model using the backward Euler method. 
% Inputs
% y0        Initial condition, derived from steady-state code (smf)
% z         Vector of depths
% m         structure of options
% dtvary    Choice of time-stepping scheme:
%                   0 fixed time step
%                   1 varying time step by velocity scale
%                   2 quadratically increasing time step
%                   3 varying time step by dpdt

% initialize variables 
be = []; flag = -14;

% reinitialize persistent variables in tdcFV
[zphigcInd] = tdcFV('calc_zphigc', 0, y0(m.blk.is.phi_g:4:end), m);

% include the pattern of the jacobian for speed.
b = ones(length(y0),40);
dFdy = spdiags(b, -20:19, length(y0), length(y0));

% options for fsolve
optsBE = optimoptions('fsolve','JacobPattern',dFdy,'Display','off',...
    'OptimalityTolerance', 1e-10, 'FunctionTolerance',1e-10,...
    'StepTolerance',1e-10,...
    'Algorithm','trust-region','SubproblemAlgorithm','factorization');

if (dtvary==1)
    [be, flag] = BE_dtVary(y0, z, m, optsBE);
else
    [be, flag] = BE_dtFixed(y0, z, m, optsBE);
end

end

function [be, flag] = BE_dtFixed (y0, z, m, optsBE)
% backward Euler with fixed time stepping

t  = linspace(m.Tspan(1),m.Tspan(2),m.Nt);
Nt = length(t);
fprintf('Number of time steps = %d.\n', Nt);

dt = t(2) - t(1);
tiOut = Nt;

yMat = zeros(length(y0), Nt);
yMat(:,1) = y0;

for ti = 2:Nt
    [yMat(:,ti), fval, exitflag, output] = ...
        fsolve(@(yN) des(t(ti), dt, yN, yMat(:,ti-1), z, m), yMat(:,ti-1), optsBE);
    
    if exitflag<=0, tiOut = ti - 1; break; end
end

be.x = t(1:tiOut);
be.y = yMat(:,1:tiOut);
be.z = z';
be.yp = [zeros(size(y0)),1./dt.*diff(yMat(:,1:tiOut),1,2)];
be.xe = [];
be.ye = [];
be.ie = [];

if ti == m.Nt,  flag = 1;
else,           flag = -1;
end

end

function [be, flag] = BE_dtVary (y0, z, m, optsBE)
% backward Euler with dt determined by dp
% rationale - variations are usually fastest in the initial part of the
% eruption, so we need smaller time steps for better approximations then.
% Important for good volume estimation in the early part of eruption

be = []; flag = -114;

dp = m.be.dp;
[dt0, dpdt] = calc_dt(m.Tspan(1), y0, m.be.dp, m);
NtMax = round(diff(m.Tspan)/dt0);
if m.slv.verbose
    fprintf('Initial time step = %.4f yr for %.2f MPa change.\n Max time steps = %d.\n', ...
        ConvertSecToYear(dt0), 1e-6*m.be.dp, NtMax);
end

% adjust dt so that it is not too small or too big.
dtRange = m.be.dtFrac*diff(m.Tspan);
if dt0<dtRange(1)
    dt0 = dtRange(1);
    dp = -dpdt*dt0;
    
    if m.slv.verbose
        fprintf('Too many time steps. Adjusted to dt = %.4f yr for dp = %.2f MPa change.\n', ...
            ConvertSecToYear(dt0), 1e-6*dp);
    end
    
elseif dt0>dtRange(2)
    dt0 = dtRange(2);
    dp = -dpdt*dt0;
    
    if m.slv.verbose
        fprintf('Too few time steps. Adjusted to dt = %.4f yr for dp = %.2f MPa change.\n', ...
            ConvertSecToYear(dt0), 1e-6*dp);
    end
end

% max length of solution
NtMax = round(diff(m.Tspan)/dt0);
if m.slv.verbose, fprintf('Adjusted max time steps = %d.\n', NtMax); end

% initialize
t     = zeros(1,NtMax);
dt    = dt0*ones(1,NtMax);
yMat  = zeros(length(y0), NtMax);

vScale = m.slv.sy(m.blk.is.v)*ones(1,NtMax);

% t = 0
ti    = 1;
t(ti) = m.Tspan(1);
yMat(:,ti) = y0;

while t(ti) < m.Tspan(2)
    
    ti = ti + 1;
    t(ti) = t(ti-1) + dt(ti-1);
    
    if t(ti) > m.Tspan(2)
        t(ti) = m.Tspan(2);
        dt(ti-1) = t(ti) - t(ti-1);
    end
    
%     fprintf('t = %.4f yr.\n', ConvertSecToYear(t(ti)));
    [yMat(:,ti), fval, exitflag, output] = fsolve(@(yN) ...
        des(t(ti), dt(ti-1), yN, yMat(:,ti-1), z, m), yMat(:,ti-1), optsBE);
    
    if any(any(imag(yMat(:,ti)))) || any(any(yMat(:,ti)<0)), exitflag = -10; end
    if exitflag<=0, break; end
    
    % find next time step using dpdt at t(ti)
    dt(ti) = calc_dt(t(ti), yMat(:,ti), dp, m);

    % if the dp is too large, tend to get wrong results so reduce dt as
    % dp/dt decreases
    if (dp>2*m.be.dp) 
        if dt(ti)>dtRange(1)
            dp = max([0.5*dp, m.be.dp]);
            dt(ti) = calc_dt(t(ti), yMat(:,ti), dp, m); 
        end
    end
%     fprintf('dp = %.2f MPa.\n', 1e-6*dp);
    
    [yMat(:,ti), m, vScale(ti)] = CheckVelScales(yMat(:,ti), yMat(:,1), m);    
end

yMat(m.blk.is.v:m.Nv:end,:) = yMat(m.blk.is.v:m.Nv:end,:).*vScale/vScale(1);

if exitflag<=0, ti = max([ti-1,1]); end

be.x = t(1:ti);
be.y = yMat(:,1:ti);
be.z = z';
be.yp = [zeros(size(y0)),1./dt(1:ti-1).*diff(yMat(:,1:ti),1,2)];
be.xe = [];
be.ye = [];
be.ie = [];

if t(ti)>=m.Tspan(2),  flag = 1;
else,                  flag = -1;
end

end

function [dt, dpdt] = calc_dt (t, y, dp, m)
pch  = y(m.blk.is.p)*m.slv.sy(m.blk.is.p);
vch  = y(m.blk.is.v)*m.slv.sy(m.blk.is.v);
beta = tdcFV('CalcCompressibility', t, pch, m);
dpdt = -pi*m.R^2*vch/m.ch.V0/beta;
dt   = -dp/dpdt;
end

function [r] = des (t, dt, yi, yiprev, z, m)
% non-linear system of equations to be solved after time discretization
% using backward Euler.

dydt = 1/dt*(yi - yiprev);
r = tdcFV('calc_mass_matrix',t,yi,z,m)*dydt - tdcFV('des',t,yi,z,m);

end

function [yNew, mNew, vScale] = CheckVelScales (y, y0, m)

v1 = y(m.blk.is.v)*m.slv.sy(m.blk.is.v);

mNew = m;
mNew.slv.sy(mNew.blk.is.v) = 10^floor(log10(v1));

vScale = mNew.slv.sy(mNew.blk.is.v);
yNew = RescaleVars(y, m, mNew.slv.sy);

end

