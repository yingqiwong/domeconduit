function varargout = tdcFV (varargin)
% Time dependent problem for conduit flow using finite volume with
% staggered grid
[varargout{1:nargout}] = feval(varargin{:});
end

function [ss, td, m, flag] = main (use_be, plot_opt)
% [ss, td, m, flag] = tdcFV('main',0,1);

if nargin == 0, use_be = 0; plot_opt = 0; end

o = tdcFV('setdef');
[ss, td, m, flag] = tdcFV('run_tdcFV', o, use_be, plot_opt);

end

function [ss, td, m, flag] = run_tdcFV (o, use_be, plot_opt)
% [ss, td, m, flag] = tdcFV('run_tdcFV',o, 0, 1);
dbstop if error

ss = []; td = []; m = []; flag = -14;

[ss, opts, ssflag] = tdcFV('run_ssc', o);
if ssflag ~= 1
    flag = ssflag; return;
end

% set options for time stepping
[m, y0, z] = tdcFV('td_init', ss.m, ss, use_be, 0);

if (opts.slv.verbose), disp('Starting time dependent solution...'); end
if (use_be)
    tic; [td, flag] = RunBE(y0, z, m, 1); toc;
else
    [td, m, flag] = tdcFV('run_tdc', y0, z, m);
end
if (opts.slv.verbose)
    fprintf('Time dependent solution solved with flag %d. \n', flag);
end

if (flag == 1) && (plot_opt == 1)
    PlotResults('plot_timeseries', td, m);
    PlotResults('PlotProfiles', td, m, []);
end

end

function [ss, opts, ssflag] = run_ssc (o)

% initialize
ss = [];

opts = tdcFV('ss_init',o);
if opts.conduit_length < 1e3, return; end

[ss, opts, ssflag] = run_ssc_opts(opts);
end

function [ss, opts, ssflag] = run_ssc_opts (opts)

if (opts.slv.verbose), disp('Starting steady state solution...'); end
tic;
[ss, ssflag] = smf_rad_dz('solve', opts);
sstime = toc;
if (opts.slv.verbose), fprintf('Elapsed time = %.4f seconds.\n', sstime); end

if (ssflag~=1)
    if (opts.slv.verbose), disp('No steady state solution'); end
    ssflag = -10; return;
elseif ss.y(2,end)>1
    if (opts.slv.verbose)
        fprintf('ss exit velocity too fast = %.3e m/s. \n', ss.y(2,end));
    end
    ssflag = -11; return;
elseif ss.y(2,end)<1e-6
    if (opts.slv.verbose)
        fprintf('ss exit velocity too low = %.3e m/s. \n', ss.y(2,end));
    end
    ssflag = -12; return;
end

if opts.slv.verbose
    fprintf('Steady state exit velocity = %.3e m/s. \n', ss.y(2,end));
    disp('Steady-state model initialized');
end
end

function [td, m, flag] = run_tdc (y0, z, m)

% initialize
td = [];

m.tStart = tic;
[td, flag] = tdcFV('des_run', y0, z, m);
tdtime = toc(m.tStart);
if (m.slv.verbose), fprintf('Elapsed time = %.4f seconds.\n', tdtime); end

end

function o = setdef ()

% default set
o.V0        = 100e9;    % m3
o.Omega     = 0;        % m3/s/Pa, factor /1e6/24/3600
o.AR        = 0.66;
o.L         = 4e3;
o.op        = 19e6;     % Pa
o.R         = 50;       % m
o.total_h2o = 5;        % wt%
o.total_co2 = 2000;     % ppm
o.phi_gc    = 0.3;
o.k_lat     = 1e-13;    % m2. actually kc (magma perm scale), but here using legacy name from older codes
o.f0        = 0.1;
o.a         = 0.01;

end

function opts = ss_init (o)
% some options
opts = smf_rad_dz('opt_setdef');

% solver options
opts.slv.use_mex    = 1;
opts.slv.max_time   = 60;
opts.slv.zero.tolx  = 1e-10;
opts.slv.de.reltol  = 1e-10;
opts.slv.de.abstol  = 1e-20;
opts.slv.verbose    = 0;
opts.plug_gas_loss  = 1;
opts.p_top          = 1e5;
opts.Nz             = 601;
opts.use_phi_s_ratio= 1;
opts.eta.use_strain = 1;
opts.newtonian.vvfrac_thr = 0.5;

opts.conduit_length = o.L;
if opts.conduit_length < 1e3, fprintf('Chamber too big.\n'); return; end
opts.R        = o.R;

% model parameters
opts.op        = o.op;    % overpressure (not truly op since rho unknown)
opts.p_ch      = 2200*9.81*opts.conduit_length + opts.op + opts.p_top;
opts.fr.v_ref  = 1e-5;    % units m/s
opts.fr.f0     = o.f0;     % dimensionless
opts.fr.a      = o.a;    % dimensionless

% chamber parameters
opts.ch.V0    = o.V0;               % m3
opts.ch.Omega = o.Omega;            % m3/s/Pa
opts.ch.pdeep = opts.p_ch;
opts.ch.AR    = o.AR;
opts.ch.a     = (3*opts.ch.V0/4/pi/(opts.ch.AR^2))^(1/3);
opts.ch.depth = o.L + opts.ch.a;    % o.chdepth
opts.ch.chi_s = 0.3116;
opts.ch.mu    = 20e9;
opts.ch.beta_fixed = 0;
opts.ch.beta  = 1e-8;
% opts.ch.beta_ch = 3/4/opts.ch.mu*(1 + 1/3*log10(1/opts.ch.AR)); % straight line interpolation
opts.ch.beta_ch = Get_beta_ch(opts.ch.AR, opts.ch.depth, opts.ch.a, opts.ch.mu);
opts.ch.x0    = 0;  %-2.160634633581430e3;
opts.ch.y0    = 0;  %-2.907450243039151e3;

opts.total_h2o = o.total_h2o;
opts.total_co2 = o.total_co2;

% degassing parameters
opts.phi_gc     = o.phi_gc;
opts.use_k_vert = 1;
opts.k_lat      = o.k_lat;
opts.k_lat_wall.L   = opts.R;
opts.k_lat_wall.top = 1e-14;

end

function [tm, y0, z] = td_init (sm, ss, use_be, mms)
% set default values for model parameters
% tm - time-dependent model options
% sm - steady-state model options

tm = sm;

% check number of variables
if sm.chi_ch.total.co2 == 0, tm.Nv = 3; tm.mw_ch = 1;
else, tm.Nv = 4; end

% options for time-dependent boundary conditions
tm.tdep.p_top_tvary = 0;
tm.tdep.p_bot_tvary = 1;
tm.tdep.h2och_tvary = 0;


tm.Tspan = ConvertYearToSec([0,5]); % years converted to seconds
tm.Nt    = 201;

% some things for finding events
% if velocity is small relative to starting velocity, break and start
% integration again
tm.vfrac = 0;
tm.vmin  = 1e-9; % stop integration if velocity drops below threshold
tm.slv.max_time  = 1000;
tm.slv.max_order = 5; % default for ode15s is 5
tm.slv.dtFixed   = 0;

% additional options for backward euler method
tm.be.use    = use_be;
tm.be.dp     = 1e6;         % for time-step determination
tm.be.dtFrac = [1e-4,0.05]; % range of acceptable dt's (here x Tspan)

%scaling for the solution components
if (mms) % for method of manufactured solutions
    tm.slv.sy = ones(4,1);
    z     = linspace(-tm.conduit_length, 0, tm.Nz)';
    tm.mw = sm.mw_ch;
    y0    = [];
    tm.v0 = 1e-3;
    
else  % for actual solutions
    tm.slv.sy = zeros(4,1);
    if (use_be) % using backward Euler method
        tm.slv.sy(tm.blk.is.p)     = 0.1*10^floor(log10(ss.y(1)));
        tm.slv.sy(tm.blk.is.v)     = 10^floor(log10(ss.y(2)));
        tm.slv.sy(tm.blk.is.phi_g) = 1;
        tm.slv.sy(tm.blk.is.mw)    = 1;
    else % using ode15s here
        tm.slv.sy = zeros(4,1);
        tm.slv.sy(tm.blk.is.p)     = 10^ceil(log10(ss.y(1)));
        tm.slv.sy(tm.blk.is.v)     = 10^ceil(log10(ss.y(2)));
        tm.slv.sy(tm.blk.is.phi_g) = 1;
        tm.slv.sy(tm.blk.is.mw)    = 1;
    end
    
    z = -ss.z(1:2:end)';
    dz = abs(ss.z(2) - ss.z(1));
    tm.Nz = length(z);
    
    % some time integration options
    tm.slv.reltol = 1e-5;
    tm.slv.abstol = [repmat([1e-20;1e-10;1e-10;1e-10],tm.Nz,1); 1e-10];      % different abstol for p,v,phi_g
    
    % fill in initial conditions
    ssy = reshape(ss.y(1:tm.Nv,1:2:end), tm.Nz*tm.Nv, 1);
    v   = [ss.y(2,2:2:end-1),interp1(ss.z, ss.y(2,:), -dz, 'pchip', 'extrap')];
    ssy(tm.blk.is.v:tm.Nv:end) = v;
    
    y0 = [ssy./repmat(tm.slv.sy(1:tm.Nv),length(z),1); 0];
    
    tm.v0 = y0(2);
    tm.zphigc = find(ss.state(1,1:2:end)==1,1);
    if isempty(tm.zphigc), tm.zphigc = tm.Nz; end
    
    tm.PlugDepth = find(ss.state(2,1:2:end) == 1,1);
    if isempty(tm.PlugDepth), tm.PlugDepth = tm.Nz; end
    %     if (use_be), tm.PlugDepth = 0; end
end

tm.dQUICKn = QUICK(z);

end

function [td, flag] = des_run (y0, z, m, MassOn)
% runs the time-dependent equations
% td = tdcFV('des_run', y0, z, m)
%
% -- Inputs --
% Tspan:time span for integration
% y0:   initial condition
% m:    model parameters
%
% -- Outputs --
% td:   solution

if (nargin<4), MassOn = 1; end

% specify Jacobian pattern for speed.
% Could be sparser if desired, here overpredicting number of non-zero elements
Nvz = m.Nv*m.Nz;
b = ones(length(y0),20);
dFdy = spdiags(b, -10:9, length(y0), length(y0));

if MassOn == 1
    % new one with symbolic inputs
    options = odeset('RelTol',m.slv.reltol, 'AbsTol',m.slv.abstol,...
        'Mass', @(t,y) calc_mass_matrix(t,y,z,m),'MassSingular', 'yes',...
        'MStateDependence', 'strong', 'MaxOrder', m.slv.max_order, ...
        'Events', @(t,y) des_event (t, y, z, m), 'Jpattern', dFdy);
else
    % with no mass matrix, i.e. ss time stepping
    is = m.blk.is;
    Mtemp = zeros(length(y0));
    Mtemp(is.p,is.p) = 1*m.slv.sy(is.p);
    Mtemp(end-m.Nv+is.p,Nvz-m.Nv+is.p) = 1*m.slv.sy(is.p);
    options = odeset('RelTol',m.slv.reltol, 'AbsTol',m.slv.abstol,...
        'Mass', Mtemp,'MassSingular', 'yes',...
        'Jpattern', dFdy);
end


t = [];  dt = (m.Tspan(2) - m.Tspan(1))/(m.Nt-1); t0 = m.Tspan(1);
y = [];  yp = [];
xe = []; ye = []; ie = [];
IntAttempt = 1;

while (1)
    
    %fprintf('Integration attempt %d\n', IntAttempt);
    
    Tspan = [t0, m.Tspan(2)];
    
    tdtmp = ode15s(@(t,y) des(t,y,z,m),Tspan,y0,options);
    
    
    if (m.slv.dtFixed)
        % Record data for this segment
        ti = ceil(t0/dt)*dt:dt:tdtmp.x(end);
        if ~isempty(ti)
            [yi, ypi] = deval(tdtmp, ti);
            t  = [t, ti];
            y  = [y, yi];
            yp = [yp, ypi];
        end
    else
        t  = [t, tdtmp.x];
        y  = [y, tdtmp.y];
%         yp = [yp, tdtmp.yp];
    end
    
    % Made it to end of integration time? If so, we're done
    % else if velocity gets too low (event 2), the solver will also break
    % else if velocity goes negative (event 3), the solver will break
    if (tdtmp.x(end) == m.Tspan(2))
        flag = 1;
        break;
    elseif (tdtmp.ie == 2)
        if (m.slv.verbose)
            fprintf('Event 2 at t=%.1e, velocity below %.1e m/s.\n',...
                tdtmp.xe, m.vmin);
        end
        flag = 2;
        break;
    elseif (any(tdtmp.ie == 3))
        if (m.slv.verbose)
            fprintf('Event 3 at t=%.1e, velocity goes negative.\n',...
                tdtmp.xe);
        end
        flag = 3;
        break;
    elseif (any(tdtmp.ie == 4))
        if (m.slv.verbose)
            fprintf('Event 4 at t=%.1e, solver timeout after %.2f seconds.\n',...
                tdtmp.xe, toc(m.tStart));
        end
        flag = 4;
        break;
    end
    
    % If not, prepare next segment
    t0 = tdtmp.x(end);
    y0 = tdtmp.y(:,end);
    m.v0 = y0(2);
    
    % rescale to real units
    yscl    = y0(1:Nvz).*repmat(m.slv.sy(1:m.Nv),length(z),1);
    p       = yscl(m.blk.is.p:m.Nv:end);
    
    % calculate pressure gradient across faces.
    % this has length Nz-1 since we only have (Nz-1) interior faces
    dz = abs(z(2) - z(1));
    dpdz = 1/dz*(p(2:end) - p(1:end-1));
    
    E = tdcFV('calc_exprs',t0,yscl,z,m,dpdz,[]);
    m.PlugDepth = E.plugdepth;
    
    xe = [xe, tdtmp.xe];
    ye = [ye, tdtmp.ye];
    ie = [ie, tdtmp.ie];
    
    IntAttempt = IntAttempt+1;
    if tdtmp.ie == 1
        if (m.slv.verbose)
            fprintf('Event 1 at t=%.1e, velocity below %d%% of IC.\n',...
                tdtmp.xe, m.vfrac*100);
        end
    end
end

td.x = t;
td.y = y;
td.z = z';
td.yp = yp;
td.xe = xe;
td.ye = ye;
td.ie = ie;

end

function F = des_implicit (t, y, yp, z, m)

MassMat = calc_mass_matrix (t,y,z,m);
RHS  = des (t, y, z, m);
F = MassMat*yp - RHS;

end

function [MassMat] = calc_mass_matrix (t,y,z,m)
% tic; massmat = tdcFV('calc_mass_matrix',0,y0,z,m); toc;
%
% uses symbolic result for mass matrix
% input y in normalized units
% update 2017 Dec 18: editted to calculate elements at each depth in
% vectorized approach

if (0)
    [MassMat] = AnlyMassMatODE('main',t,y,z,m);
    return;
end

is = m.blk.is;
Nvz  = m.Nv*m.Nz;
yscl = y(1:Nvz).*repmat(m.slv.sy(1:m.Nv),length(z),1);
p = yscl(is.p:m.Nv:Nvz);

% make vector of chibreaks and coefbreaks for calculating chi_s in the
% mass matrix function
breaks = m.pf.pp.breaks;
coefs  = m.pf.pp.coefs;

chibreaks = zeros(size(p));
coefbreak  = zeros(length(p), 4);

for BreakInd = length(breaks):-1:2
    % which chibreak for each pressure
    IndReplace = p<1e6*breaks(BreakInd);
    NReplace = length(IndReplace(IndReplace==1));
    
    % replace elements in chibreaks vector
    chibreaks(IndReplace)   = breaks(BreakInd-1);
    coefbreak(IndReplace,:) = repmat(coefs(BreakInd-1,:),NReplace,1);
end

% interior points
[Diags, k] = MassMatFunc_getDiags(yscl,m,chibreaks,coefbreak);
MassMat    = sparse(length(y), length(y));
MassMat(1:Nvz, 1:Nvz) = spdiags(Diags, k, Nvz, Nvz);

% equations at boundary/boundary conditions
MassMat(is.p, :)    = 0;
MassMat(is.p, is.p) = 1*m.slv.sy(is.p);
MassMat(is.phi_g,:) = 0;

zi = length(z);
Row = m.Nv*(zi-1);
MassMat(Row+is.v, Row+is.p) = 1*m.slv.sy(is.p);
MassMat(end,end) = 1e6;

if m.Nv == 4
    MassMat(is.mw,:) = 0;
end
end


function [val, isterm, dir] = des_event (t, y, z, m)

is = m.blk.is;

% if the velocity drops below some threshold, it's faster to break and
% integrate over a new time interval with different initial conditions
val(1) = y(is.v)/m.v0 - m.vfrac;

% if the base velocity drops too low
val(2) = (y(is.v)*m.slv.sy(is.v) < m.vmin) - 1;

% if any solution component goes negative
% val(3) = any(y(is.v:m.Nv:end)<0) - 1;
val(3) = 1;

% time out if the solver takes too long
if isfield(m,'tStart')
    val(4) = double((toc(m.tStart) - m.slv.max_time)<0);
end

isterm = ones(size(val));
dir = [];

end

function dydt = des (t, y, z, m)
% Function that contains the rhs of the time derivatives
% des_dydt = tdcFV('des',0,y0,z,m);
%
% -- Inputs --
%  t:   time
%  y:   vector of 4 unknowns: [p,v,phi_g,mw]x Nz + Flux
%  z:   depth of evaluation
%  m:   model parameters
%
% -- Outputs --
% dydt: time derivatives

is = m.blk.is;
Nv = m.Nv;
Nvz= m.Nv*m.Nz;
dz = abs(z(2) - z(1));
dydt = zeros(size(y));

% rescale to real units
yscl    = y(1:Nvz).*repmat(m.slv.sy(1:m.Nv),length(z),1);
p       = yscl(is.p:Nv:end);
v       = yscl(is.v:Nv:end);
phi_g   = yscl(is.phi_g:Nv:end);
if m.Nv == 3, mw = m.mw_ch;
else,         mw = yscl(is.mw:Nv:end); end

% calculate pressure gradient across faces.
% this has length Nz-1 since we only have (Nz-1) interior faces
dpdz = 1/dz*(p(2:end) - p(1:end-1));
dpdzn = [dpdz; 2*dpdz(end) - dpdz(end-1)];
dpdzs = [2*dpdz(1)-dpdz(2); dpdz];

% face velocities for continuity cells
vn = v(2:end);
vs = v(1:end-1);

E = tdcFV('calc_exprs',t,yscl,z,m,dpdz,[]);

[nvn,   nvs]       = tdcFV('get_FaceVals', E.nv.mass,  m.dQUICKn);
[h2oliqn, h2oliqs] = tdcFV('get_FaceVals', E.h2o.liq,  m.dQUICKn);
[h2ogasn, h2ogass] = tdcFV('get_FaceVals', E.h2o.gas,  m.dQUICKn);
[kmagn, kmags]     = tdcFV('get_FaceVals', E.kvert/m.eta_g, m.dQUICKn);

vgn = -kmagn.*dpdzn;
vgs = -kmags.*dpdzs;

h2o.liq = 1/dz*(h2oliqn(2:end).*vn - h2oliqs(2:end).*vs);

h2o.gas = 1/dz*(h2ogasn(2:end).*(vn+vgn(2:end)) - h2ogass(2:end).*(vs+vgs(2:end)))...
    + E.h2o.lat(2:end);

% interior points
dydt(Nv+is.p:Nv:Nvz)     = -1/dz*(nvn(2:end).*vn - nvs(2:end).*vs);
dydt(Nv+is.phi_g:Nv:Nvz) = -h2o.liq - h2o.gas;

dydt(is.v:Nv:Nvz-Nv) = v(1:end-1) - E.v.vv - E.v.vf;

if m.Nv == 4
    [co2liqn, co2liqs] = tdcFV('get_FaceVals', E.co2.liq, m.dQUICKn);
    [co2gasn, co2gass] = tdcFV('get_FaceVals', E.co2.gas, m.dQUICKn);
    
    co2.liq = 1/dz*(co2liqn(2:end).*vn - co2liqs(2:end).*vs);
    co2.gas = 1/dz*(co2gasn(2:end).*(vn+vgn(2:end)) - co2gass(2:end).*(vs+vgs(2:end)))...
        + E.co2.lat(2:end);
    
    
    dydt(Nv+is.mw:Nv:Nvz) = -co2.liq - co2.gas;
end

if m.plug_gas_loss
    plugdepth = m.PlugDepth;
    
    if m.PlugDepth == 0
        vvfrac = E.v.vv./(E.v.vv + E.v.vf);
        plugdepth = find(vvfrac < m.newtonian.vvfrac_thr, 1);
    end
    
    pluglog  = 1./(1+exp(-0.05*(z-z(plugdepth))));
%     pluglog = zeros(size(p)); pluglog(z>z(plugdepth)) = 1;
    
    [h2ogasn2, h2ogass2] = tdcFV('get_FaceVals', 1./(1+E.Gamma).*phi_g,  m.dQUICKn);
    h2o.plug = E.rho_g(2:end)/dz.*(h2ogasn2(2:end).*vn - h2ogass2(2:end).*vs);
    
    waterRHS = -h2o.plug;
    indstart = 1;
    h2oind   = (plugdepth-indstart)*Nv + is.phi_g;
    %     dydt(h2oind:m.Nv:end) = waterRHS(plugdepth-indstart:end);
    dydt(h2oind:m.Nv:Nvz) = (1-pluglog((plugdepth-indstart+1):end)).*dydt(h2oind:m.Nv:Nvz) + ...
        pluglog((plugdepth-indstart+1):end).*waterRHS((plugdepth-indstart):end);
    
    if m.Nv == 4
        [co2gasn2, co2gass2] = tdcFV('get_FaceVals', ...
            E.Gamma./(1+E.Gamma).*phi_g,  m.dQUICKn);
        co2.plug = E.rho_g(2:end)/dz.*(co2gasn2(2:end).*vn - co2gass2(2:end).*vs);
        
        co2RHS   = - co2.plug;
        co2ind = (plugdepth-indstart)*Nv + is.mw;
        %         dydt(co2ind:m.Nv:end) = co2RHS(plugdepth-indstart:end);
        dydt(co2ind:m.Nv:Nvz) = (1-pluglog((plugdepth-indstart+1):end)).*dydt(co2ind:m.Nv:Nvz) + ...
            pluglog((plugdepth-indstart+1):end).*co2RHS((plugdepth-indstart):end);
        
    end
end


% boundary conditions
% chamber boundary
[mw_ch, phi_g_ch] = calc_ch_gas(t, p(1), m);
[beta]            = CalcCompressibility(t, p(1), m);
dydt(is.phi_g) = phi_g(1) - phi_g_ch;
if m.Nv == 4,     dydt(is.mw) = mw(1) - mw_ch;  end

vch            = 1.5*v(1) - 0.5*v(2);
if m.tdep.p_bot_tvary == 1
    dydt(is.p) = ((m.ch.Omega*(m.ch.pdeep-p(1)) - pi*m.R^2*vch)/m.ch.V0/beta);
    %     tc = ConvertYearToSec(1); dydt(is.p) = -1e6/tc*sin(t/tc);
else
    dydt(is.p) = 0;
end

% top boundary condition
if m.tdep.p_top_tvary == 1
    tcdome = ConvertYearToSec(0.09);
    dydt(Nvz-Nv+is.v) = 5050188/tcdome*exp(-t/tcdome);
else
    dydt(Nvz-Nv+is.v) = 0;
end
vExit = 0.5*(v(end-1)+v(end));
dydt(end) = pi*m.R^2*vExit*(1-phi_g(end));

if isfield(m, 'Source')
    % for method of manufactured solutions
    % input the correct bc for dpdt(1)
    switch m.Source.Type
        case 'zexp_tsine'
            dydt(is.p) = m.pdiff/m.p_tc*cos(t/m.p_tc);
        case 'zexp_texp'
            dydt(is.p) = -m.pdiff/m.p_tc.*exp(-t/m.p_tc);
    end
    
    zSrc = linspace(-m.conduit_length, 0.5*dz, 2*(m.Nz-1)+2)';
    [y_ex,dydt_ex,dydz_ex,dpdz2,m] = MMSTestSols(zSrc, t, m, m.Source.Type, 0);
    
    SymSource = MMSFuncs('CalcSourceFromSyms', zSrc, t, y_ex, dydt_ex, dydz_ex, dpdz2, m);
    SymSourceStag = MMSFuncs('ConvertCellCenterToStagger', SymSource, m);
    dydt = dydt + SymSourceStag;
end

end

function E = calc_exprs (t, y, z, m, dpdz, AllFields)
% similar to first part of des_calc_exprs in smf_rad.m.
% This function calculate dependent variables other than those listed in y
% E = tdcFV('calc_exprs',0,yscl,ztd,m,dpdz,1);
%
% -- Inputs --
%  y:     vector of 4 unknowns in real SI units: [p,v,phi_g,mw]
%  z:     depth of evaluation
%  m:     model parameters
%  degas: switch for degassing. =1 when percolation threshold is reached.
%
% -- Outputs --
%  E:     structure containing various properties

if (isempty(AllFields)), AllFields = 0; end

is = m.blk.is;
Nv = m.Nv;

p     = y(is.p:Nv:end);
v     = y(is.v:Nv:end);
phi_g = y(is.phi_g:Nv:end);
if m.Nv == 3, mw = m.mw_ch;
else,         mw = y(is.mw:m.Nv:end); end

zphigc = calc_zphigc(t, phi_g, m);
degas  = 1./(1+exp(-0.04*(z-z(zphigc))));

% if m.PlugDepth == 0, degas = zeros(size(p)); degas(z>z(zphigc)) = 1; end

if isfield(m, 'Source')
    degas = ones(size(p));
end

% Volatile mass fractions.
[Chi_hd, Chi_cd] = solubility_liu_explicit(p, m.T, mw);
Chi_hd = 1e-2*Chi_hd;
Chi_cd = 1e-6*Chi_cd;
c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
c12 = c1.*c2;

% solid fraction
[chi_s] = tdcFV('calc_chi_s', m.pf, p*1e-6);

% simulate cooling by increasing chi_s with time
if isfield(m, 'cooling')
    if m.cooling == 1
        Fac = 1 + 0.2/ConvertYearToSec(4)*t;
        chi_s = Fac*chi_s;
        chi_s(chi_s>1) = 1;
    end
end

E.chi_s = chi_s;
phi_s = chi_s.*m.rho_l.*c12.*(1 - phi_g)./...
    (m.rho_s + chi_s.*(c12*m.rho_l - m.rho_s));

% liquid fraction
phi_l = (1 - phi_s - phi_g).*c2;
E.phi_l = phi_l;

% gas density
rho_g = p.*(mw./(m.Rw*m.T) + (1 - mw)./(m.Rc*m.T));
E.rho_g = rho_g;

% bulk density
rho = m.rho_l*phi_l.*c1 + m.rho_s.*phi_s + rho_g.*phi_g;
phi_s_eta = phi_s./(1-phi_g);

Gamma = (1 - mw)./(m.B*mw);
E.Gamma = Gamma;

% gas loss parameters
phyd = abs(m.phydro.slope*z);
[kvert, klat] = calc_k_allz(m, phi_g, phi_l, z, degas, m.plug_gas_loss);
E.kvert = kvert;

ug = klat/m.eta_g.*(p-phyd)./(m.R+m.klw.L);

% components in mass conservation equations
E.nv.mass  = m.rho_l*phi_l + m.rho_s*phi_s;

E.h2o.liq = Chi_hd*m.rho_l.*phi_l.*c1;
E.co2.liq = Chi_cd*m.rho_l.*phi_l.*c1;

E.h2o.gas =     1./(1+Gamma).*rho_g.*phi_g;
E.co2.gas = Gamma./(1+Gamma).*rho_g.*phi_g;

E.h2o.lat  = 2*       rho_g.*phi_g.*ug./m.R./(1+Gamma);
E.co2.lat  = 2*Gamma.*rho_g.*phi_g.*ug./m.R./(1+Gamma);

% prep varables for MBE
[rhon, ~]       = tdcFV('get_FaceVals', rho,        m.dQUICKn);
[phi_s_etan, ~] = tdcFV('get_FaceVals', phi_s_eta,  m.dQUICKn);
[Chihdn, ~]     = tdcFV('get_FaceVals', Chi_hd,     m.dQUICKn);

if (m.eta.use_strain==1)
    gdot = 2*abs(v)/m.R;
elseif (m.eta.use_strain==0)
    gdot = zeros(size(v));
else
    gdot = zeros(size(v));
end

[eta,eta_m,eta_phi] = calc_eta(Chihdn(1:end-1), phi_s_etan(1:end-1), m.T, m.eta.melt, gdot(1:end-1));

% momentum balance
tauR    =-0.5*m.R*(dpdz + rhon(1:end-1)*m.g);
E.v.vv  = 0.25*tauR*m.R./eta;
zMBE    = z(1:end-1) + 0.5*(z(2)-z(1));
sig_c   = m.sig.slope*abs(zMBE) + m.sig.offset;
E.v.vf  = 2*m.fr.v_ref*exp(-m.fr.f0/m.fr.a)*sinh(tauR./m.fr.a./sig_c);

vvfrac = E.v.vv./(E.v.vv + E.v.vf);
E.plugdepth = find(vvfrac < m.newtonian.vvfrac_thr, 1);

% output fields if requested
if AllFields==1
    E.dpdz  = dpdz;
    E.eta   = eta;
    E.eta_m = eta_m;
    E.eta_phi = eta_phi;
    E.phi_s = phi_s;
    E.phi_s_eta = phi_s_eta;
    E.phi_l = phi_l;
    E.phi_g = phi_g;
    E.rho_g = rho_g;
    E.rho   = rho;
    E.c1    = c1;
    E.c2    = c2;
    E.Gamma = Gamma;
    E.gdot  = gdot;
    E.tauR  = tauR;
    E.Chi_hd= Chi_hd;
    E.Chi_cd= Chi_cd;
    E.k_lat = klat;
    E.ug    = ug;
    E.zphigc = zphigc;
    E.degas = degas;
    E.zMBE  = zMBE;
    E.sig_c = sig_c;
end

end

function [zphigcInd] = calc_zphigc (t, phi_g, m)
% function that returns the percolation threshold depth.
% It remembers the percolation threshold depth from a previous time step,
% checks the new depth at later time, and adjusts if deeper

% ensures variable values are remembered in this fn for later time
persistent tprev zphigc_prev

% initialize
if isempty(tprev) || t == 0
    tprev = t;
    zphigc_prev = m.zphigc;
    %     zphigc_prev = find(phi_g > m.phi_gc, 1);  % may not equal m.zphigc
    %     if isempty(zphigc_prev), zphigc_prev = m.Nz; end
end

if t>tprev  % ensures that we only adjust for later time
    
    zphigcInd = find(phi_g > m.phi_gc, 1);
    
    % only adjust if the new threshold depth is deeper
    if zphigcInd > zphigc_prev
        zphigcInd = zphigc_prev;
    elseif isempty(zphigcInd)
        zphigcInd = zphigc_prev;
    end
    
    % reassign new time and percolation depth
    tprev = t;
    zphigc_prev = zphigcInd;
    
else        % maintain some depth if time is = or < tprev
    zphigcInd = zphigc_prev;
end

end


function [Fn, Fs] = get_FaceVals (F, dQUICKn)
% get values at north and south faces of each cell

Nz = length(F);
% dQUICKs = [1.5, -0.5, zeros(1,Nz-2); dQUICKn(1:end-1,:)];
% dQUICKs = [15/8, -5/4, 3/8, zeros(1,Nz-3); dQUICKn(1:end-1,:)];

Fn = dQUICKn*F;
% Fn(1) = interp1(2:20, Fn(2:20), 1, 'spline', 'extrap');
% Fs = [1.5*F(1) - 0.5*F(2); Fn(1:end-1)];
% Fs(1) = interp1(2:20, Fs(2:20), 1, 'spline', 'extrap');
Fs = [15/8*F(1) - 5/4*F(2) + 3/8*F(3); Fn(1:end-1)];

end

function [dQUICKn] = QUICK (z)
% [dQUICKn] = tdcFV('QUICK',z);

% using QUICK to estimate variable values at the north cell face,
% using cell-based integrals

Nz = length(z);

% DiagVecs = [-1/6*ones(Nz,1), 5/6*ones(Nz,1), 2/6*ones(Nz,1)];
DiagVecs = [-1/8*ones(Nz,1), 3/4*ones(Nz,1), 3/8*ones(Nz,1)];

dQUICK = spdiags(DiagVecs, -1:1, Nz, Nz);
dQUICKn = dQUICK;
dQUICKn(1,:) = [0.5, 0.5, zeros(1,Nz-2)];
% dQUICKn(1,:) = [3/8, 6/8, -1/8, zeros(1,Nz-3)];
% dQUICKn(end,:) = [zeros(1,Nz-2), -0.5, 1.5];
dQUICKn(end,:) = [zeros(1,Nz-3), 3/8, -5/4, 15/8];

end

function [mw_ch, phi_g_ch] = calc_ch_gas (t, pch, m)
% evaluates the two equations that relate total water and CO2 content in
% the chamber to mw_ch and phi_g_ch.
% if total volatile content changes with time, expand input variables
% [mw_ch, phi_g_ch] = tdcFV('calc_ch_gas', t, pch, m);

mw_ch = fzero(@(mw) calc_mw_ch(t,mw,m,pch),[1e-2,1]);
% if m.tdep.h2och_tvary == 1      % h2o varying with time
%     mw_ch = m.mw_ch;
% end

[Chi_hd, Chi_cd] = solubility_liu_explicit(pch, m.T, mw_ch);
Chi_hd = 1e-2*Chi_hd;
Chi_cd = 1e-6*Chi_cd;
chi_ch = CalcVolatileContent(t, m);
total_exsolved = 1e-2*chi_ch.h2o + ...
    1e-6*chi_ch.co2 - Chi_hd - Chi_cd;
c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
c12 = c1.*c2;
[chi_s] = calc_chi_s (m.pf,pch*1e-6);
rho_g = mw_ch.*pch/(m.Rw*m.T) + (1 - mw_ch).*pch/(m.Rc*m.T);
d = (total_exsolved.*m.rho_l.*c12).*...
    (1 - chi_s.*m.rho_l.*c12./(m.rho_s + chi_s.*(c12.*m.rho_l - m.rho_s)));

phi_g_ch = d/(rho_g + d);

end

function r = calc_mw_ch (t, mw, m, pch)
Gamma = (1 - mw)./(m.B*mw);
[Chi_hd Chi_cd] = solubility_liu_explicit(pch, m.T, mw);
chi_ch = CalcVolatileContent(t, m);
r = Gamma*1e-2*(chi_ch.h2o - Chi_hd) -...
    1e-6*(chi_ch.co2 - Chi_cd);
end

function [chi_ch] = CalcVolatileContent (t, m)
% determine the total water content at conduit base

chi_ch = m.chi_ch.total;

if m.tdep.h2och_tvary == 1      % h2o varying with time
    chi_ch.h2omin = 4;
    chi_ch.tc = ConvertYearToSec(2);
    chi_ch.h2o = chi_ch.h2omin + (m.chi_ch.total.h2o-chi_ch.h2omin)*exp(-t/chi_ch.tc);
end

end

function [beta] = CalcCompressibility (t, pch, m)
% calculates the system compressibility (Pa^-1).

if m.ch.beta_fixed == 1
    beta = m.ch.beta;
else
    
    % find the pressure at the chamber center, which determines density
    pMax  = 2*(m.ch.depth*m.sig.slope + m.sig.offset);
    pcc   = fzero(@(pcc) ChamberCenterPressure(t, pcc, pch, m), [pch, pMax]);
    rhocc = CalcRho(t, pcc, m);
    
    % magma compressibility calculated using finite difference
    pdiff      = pcc + [-1e6, 1e6];
    rhodiff(1) = CalcRho(t, pdiff(1), m);
    rhodiff(2) = CalcRho(t, pdiff(2), m);
    
    beta_mag = 1/rhocc*(diff(rhodiff)/diff(pdiff));
    beta     = (beta_mag + m.ch.beta_ch);   % system compressibility
end

end

function [r] = ChamberCenterPressure (t, pcc, pch, m)
% returns the residual for fzero to calculate the pressure at the center of
% the chamber

[rho] = CalcRho(t, pcc, m);
r = pcc - pch - rho*m.g*m.ch.a;
end

function [rho] = CalcRho (t, p, m)
% calculates the density of the magma given pressure at center of chamber

[mw, phi_g] = calc_ch_gas(t, p, m);

% Volatile mass fractions.
[Chi_hd, Chi_cd] = solubility_liu_explicit(p, m.T, mw);
Chi_hd = 1e-2*Chi_hd;
Chi_cd = 1e-6*Chi_cd;
c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
c12 = c1.*c2;

% solid fraction
% [chi_s] = tdcFV('calc_chi_s', m.pf, p*1e-6);
chi_s = m.ch.chi_s; % fixed chi_s.  Correct? Otherwise beta_mag < 0
phi_s = chi_s.*m.rho_l.*c12.*(1 - phi_g)./...
    (m.rho_s + chi_s.*(c12*m.rho_l - m.rho_s));

% liquid fraction
phi_l = (1 - phi_s - phi_g).*c2;

% gas density
rho_g = p.*(mw./(m.Rw*m.T) + (1 - mw)./(m.Rc*m.T));

% bulk density
rho = m.rho_l*phi_l.*c1 + m.rho_s.*phi_s + rho_g.*phi_g;

end

function [chi] = calc_chi_s (pf,p)
% similar to pfl_val in smf_rad
% [chi] = tdcFV('calc_chi_s',m.pf,1e-6*s.p);
%
% -- Inputs --
%  pf:      Structure containing pf.p[MPA],pf.chi[fraction],pf.pp(polynomial)
%  p:       Pressure [MPA]

% -- Outputs --
% chi:      Non-vapor mass fraction of solids

chi = zeros(size(p));
chi(p <= pf.p(1)) = pf.chi(1);
chi(p >= pf.p(end)) = pf.chi(end);
mask = p > pf.p(1) & p < pf.p(end);
chi(mask) = ppval(pf.pp, p(mask));
assert(all(chi <= 1));

end

function [kvert, klat] = calc_k_allz(m, phi_g, phi_l, z, degas, plug_gas_loss)

kmag = m.k_lat.*(phi_g.^3).*degas;

if (m.klw.use && any(kmag > 0))
    k_lat_wall = m.klw.top./(1e-3*abs(z)).^m.klw.mi;
    klat = (m.klw.L+m.R)./(m.R./kmag + m.klw.L./k_lat_wall);
else
    klat = zeros(size(phi_g));
end

if m.use_k_vert == 0
    kvert = zeros(size(phi_g));
else
    kvert = kmag;
end

end


function k_lat = calc_k_lat (m, z, phi_g, magma_on)
kmag = calc_k(m, phi_g, magma_on);
% What was computed so far is k_lat_magma. If ~isinf(k_lat_wall), then the
% wall has to be accounted for.
if (m.klw.use && any(kmag > 0))
    k_lat_magma = kmag;
    k_lat_wall = m.klw.top./(1e-3*abs(z)).^m.klw.mi;
    k_lat = (m.klw.L+m.R)./(m.R./k_lat_magma + m.klw.L./k_lat_wall);
end
end

function k = calc_k (m, phi_g, magma_on)
k = m.k_lat.*(phi_g.^3).*magma_on;
end

function klw = klw_init (oklw)
klw.use = ~isinf(oklw.top);
if (~klw.use) return; end
klw.top = oklw.top;
klw.mi = oklw.mi;
klw.L = oklw.L;
end

function b = is_thr_perm (k_lat_model)
b = any(k_lat_model == [2 3]);
end

function [eta,eta_m,eta_phi] = calc_eta (c, phi, T, meltmodel, gdot)
% This function calculates the viscosity of a silicate melt as a function of
% dissolved water content (Hess & Dingwell) and crystal content (Costa).
%
%  INPUTS
%       c:    Dissolved water content (ie., .05, *not* 5%), in MASS
%             CONCENTRATION (relative to melt phase only). This is normal.
%       phi:  VOLUME fraction of crystals in the melt. Do not confuse with
%             mass fraction!
%       T:    Melt temperature, in KELVIN
%  OUTPUTS
%       eta:  Viscosity in Pa-s
%
%  From K. Anderson's routine.

switch meltmodel
    case {1}
        % Hess and Dingwell 1996:
        lc = log(100*c);
        eta_m = 10.^((-3.545 + 0.833*lc) + (9601 - 2368*lc)./...
            (T - (195.7 + 32.25*lc)));
    case {2}
        % Whittington et al 2009
        c = 100*c;
        eta_m = 10.^(-4.43 + (7618.3-17.25*log10(c+0.26))./...
            (T-406.1+292.6*log10(c+0.26)));
    otherwise
        error('not a eta_model');
end

if gdot==0
    % From Costa 2005b and 2007b
    kappa = 0.999916;
    phistar = 0.673;
    gamma = 3.98937;
    delta = 16.9386;
else
    % From Costa 2005 and Caricchi 2007
    phistar = -0.066499*tanh(0.913424*log10(gdot)+3.850623) + 0.591806;
    delta   = -6.301095*tanh(0.818496*log10(gdot)+2.86)     + 7.462405;
    kappa   = -0.000378*tanh(1.148101*log10(gdot)+3.92)     + 0.999572;
    gamma   =  3.987815*tanh(0.8908*log10(gdot)+3.24)       + 5.099645;
end


pih = 1.772453850905515882;
B = 2.5;

eta_phi = (1 + (phi./phistar).^delta) ./ ...
    (1 - kappa.*erf((pih.*phi)./(2*kappa.*phistar).* ...
    (1 + (phi./phistar).^gamma)) ...
    ).^(B*phistar);


eta = eta_m.*eta_phi;

end

function [h2o_d co2_d h2o_e co2 co2_e] = solubility_liu_explicit ( ...
    p, T, Mh, h2o, B)
%  Given pressure, temperature, mole fraction H20, and total H20, this function
%  will calculate H2O and CO2 solubilities (this does not require the total
%  water). Then, using total water and the relationship between mole fractions
%  and exsolved volatile content, the function will calculate remaining unknowns
%  such that we have total, dissolved, and exsolved concentrations for both H2O
%  and CO2. Expressions come from Liu et al [2005], and see Anderson and Segall
%  [2011].
%
%  This function is vectorized in p.
%
%  WARNING: Careful with units!
%
%  -- Inputs --
%  p:               Pressure [PA]
%  T:               Temperature [KELVIN]
%  Mh:              Mole fraction water [FRACTION FROM 0 to 1]
%  h2o:             Total water [WT %]
%  B:               Ratio of molecular mass of H2O to CO2
%
%  -- Outputs --
%  h2o_d            Dissolved water [WT%]
%  h2o_e            Exsolved water [WT%]
%  co2              Total CO2 [PPM]
%  co2_d            Dissolved CO2 [PPM]
%  co2_e            Exsolved CO2 [PPM]
%
% April 18, 2013: Initial coding (based on much earlier code) (KA)

% Convert to MPa (as required by Liu expressions).
p = p*1e-6;

% Precompute for speed.
Pw = p.*Mh;
Pc = p.*(1-Mh);
sqrtPw = sqrt(Pw);
Pw15 = Pw.*sqrtPw; % Pw^1.5

% These equations require MPa, Kelvin, and return wt% and ppm in weight.
% Assuming saturation!
h2o_d = (354.94*sqrtPw + 9.623*Pw - 1.5223*Pw15)/T + ...
    0.0012439*Pw15 + Pc.*(-1.084e-4*sqrtPw - 1.362e-5*Pw);
co2_d = Pc.*(5668 - 55.99*Pw)/T + Pc.*(0.4133*sqrtPw + 2.041e-3*Pw15);

if (nargin == 5)
    h2o_e = h2o - h2o_d; % exsolved water
    co2_e = (h2o_e/B.*(1 - Mh)./Mh)*1e-4; % 1e-4 converts from wt% to ppm
    co2 = co2_d + co2_e;
    
    Ih = h2o < h2o_d;
    Ic = co2 < co2_d;
    if (sum(Ih + Ic) > 0), keyboard; warning('Melt is undersaturated!'); end
end
end






