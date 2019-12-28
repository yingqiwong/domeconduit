function varargout = smf_rad_dz (varargin)
% (F)ull steady-state conduit model solved using the (s)hooting (m)ethod.
%
% o = opt_setdef(): Default options.
% o = opt_convert(o): Convert from slvof options.
% de = solve(o, v_ch_guess, verbose, check_options):
%   Solve the problem. v_ch_guess is optional. verbose is 0 (off) by
%   default. Pass v_ch_guess = [] if you want to specify verbose = 1 but not
%   provide v_ch_guess. check_options is also optional; it checks that the
%   options in o are OK.
% de = add_fields(de): Add a bunch of fields to the de struct.

% Model and code by K. Anderson, A.M. Bradley (ambrad@cs.stanford.edu),
%                   P. Segall (segall@stanford.edu).
% This file:
%   Oct 2013. AMB. Initial.
%   Jun 2014. AMB. Added a lot of comments and msg output. Made
% add_measurement_fields use add_fields if use_mex is 0. Set tauY to 0 rather
% than NaN when visc_only is 1.
%   1 Jul 2014. AMB. visc_only -> model = 'newtonian'. Plug gas loss for
% 'newtonian'.
%   3 Apr 2017 - YQW.  Editted to have equal z steps in solution. Different
%   from file in smf_2016 folder.
  
%   Coordinates:
%     z_norm = z / m.conduit_length
%     x = 1 - z_norm.
%
%   A key quantity is F:
%     zeta(x) = c/(a + b x);
%     F(x) = 0 if zeta(x) >= 1  (Pure plug.)
%            1 if zeta(x) <= 0  (Pure viscous flow.)
%            1 + (zeta^4 - 4 zeta)/3 else;  (Mixed.)
%
%   State variables track switches in quantities that lead to tricky C^0
% solutions. (Not every cause of a C^0 solution is associated with a state
% and ode event.)
%     k (permeability) state is 0 (k = 0),
%                               1 (k can be > 0).
%     plug state is 1 (plug is total and plug_gas_loss is being used and the
%                      effect is active),
%                   0 (otherwise).

% I. Math outline.
%   The system of continuity equations f in fields y(z) has the form
%     d f(z, y)/dz = g(z, y).
% The change in f(z, y) within a parcel dz is balanced by the source g(y). To
% turn the continuity equations into a form suitable for odme15i, ie, an ODE in
% dy/dz, take the partial derivative w.r.t. y and z:
%     f_y(z, y) dy/dz + f_z(z, y) = g(z, y).
%   The other equation is the momentum balance equation (MBE):
%     vv F + vf (1 - F)^p = v.
% vv, vf, and F all depend on dp/dz. Together, these equations create an index-1
% DAE. We turn the DAE into an implicit ODE by solving the MBE ourselves at each
% call from ode15i. (Testing of various approaches shows this one is the best.)
%   In this implementation, y = [p, v, phi_g, mw]. The implicit ODE system
%     F(y, dy/dz, z) = 0
% is solved to integrate y from the chamber to the surface. F has four
% equations. The first is
%     dp/dz - x = 0,
% where x is obtained by solving the MBE for dp/dz. The other three are
%     f_y(z, y) dy/dz + f_z(z, y) - g(z, y) = 0.
%   We are solving a 2-point boundary value problem. v_ch is unspecified at
% the chamber, and p must match p_top at the surface. We solve this problem
% using the shooting method. Define the scalar function
%     fsm(v_ch) = p(z = surface) - p_top.
% The shooting method solves the equation
%     fsm(v_ch) = 0
% for v_ch.
%
% II. Code outline.
%   solve is the top-level solver routine. Its purpose is to prevent any
% otherwise uncaught errors from disrupting MCMC or something like it.
%   solve calls solve1. It sets up and solves fsm(v_ch) = 0 for v_ch. It
% initializes the model from the user options struct. It find an interval
% [v_ch_low, v_ch_high] such that sign(f(v_ch_low)) ~= sign(f(v_ch_high)). Then
% it calls fzero with the interval and a handle to the function
% sm_calc_fsm. fzero tries to find v_ch so that f(v_ch) = 0. Finally, it
% deals with various exit conditions.
%   solve_calc_fsm integrates the DAE from the chamber to the surface if it
% can. If pressure goes negative, the DAE stops being integrated, and
% solve_calc_fsm extrapolates a value to the surface. In any case, some value
%     fsm(v_ch) = p(surface) - p_top
% is returned.
%   sm_calc_fsm calls de_run to solve the DAE. sm_calc_fsm uses ode15i. The
% state variables for permeability k and the plug segment the conduit. These
% segments are not known in advance. sm_calc_fsm use a while loop to march up
% the conduit. In each iteration of the while loop, a new segment is
% integrated until a state variable changes based on an event. Three function
% handles are passed to ode15i: des_calc_F, des_calc_F_jac,
% des_calc_event. The 2nd and 3rd are passed through odeset. Before ode15i
% can be called, consistent initial conditions y and yp must be
% found, where yp is dy/dz. des_init solves a linear equation to obtain yp in
% terms of y.
%   des_calc_F calculates F(y, yp, z). It calls des_calc_exprs to evaluate all
% the science formulas. It calls des_calc_f_{y|z} to approximate f_y(z, y) and
% f_z(z, y) by finite differences.
%   des_calc_F_jac calculates the Jacobian matrices F_y and F_yp.
%   des_calc_event determines when state variables switch values.
%   des_calc_exprs calls
%     solubility_liu_explicit to calculate Chi_hd and Chi_cd,
%     pfl_val to calculate phi_s(p),
%     calc_eta to evaulate viscocity,
%     calc_yield to calculate the yield strength,
%     calc_k to calculate permeability, and
%     mb_solve to solve the momentum balance equation for dp/dz.
%
% III. Data structures.
%   o is the struct of user-controlled options.
%   m is the model initialized by o and v_ch.
%   de is a DAE solution for the whole conduit.
%   des (internal) is a DAE solution for one segment.
%   pf (internal) contains the data for the chi-pressure curves.

  [varargout{1:nargout}] = feval(varargin{:});
end

function o = opt_setdef ()
  o.Tindex = 3;
  o.p_top = 1e5;
  
  o.v_init_exp = -4;
  o.p_ch = []; % If provided, this has precedence over v_init_exp.
  o.total_h2o = 5;    % percent
  o.total_co2 = 2000; % ppm
  o.R = 100;
  o.conduit_length = 3500;
  
  % Solver.
  o.slv.fd = 1e-6;
  o.slv.fdo = 2;
  o.slv.zero.tolx = 1e-5;
  o.slv.de.reltol = 1e-5;
  o.slv.de.abstol = 1e-10;
  o.slv.verbose = 1;
  % Max time [s] per solve attempt.
  o.slv.max_time = 60;
  o.Nz = 1001;  % number of points in depth
  o.slv.use_mex = 1;
  [~,cpname] = system('hostname');
  if strcmp(cpname(1:4), 'cees')
      o.slv.mexFunc = 'smf_rad_mex_cees';
  else
      o.slv.mexFunc = 'smf_rad_mex';
  end

  % Permeability.
  o.k_lat_model = 2;
  o.L_diff = 1; % [m]
  % Magma only.
  o.k_lat = 1e-15; % [m^2]
  o.phi_gc = 0.3;
  % Wall only.
  o.k_lat_wall.L = 100; % [m] length scale for pw-phyd
  o.k_lat_wall.top = inf;
  o.k_lat_wall.mi = 3.2; % Manning-Ingebritsen exponent
  % If used, overrides k_lat when F < 1.
  o.plug_gas_loss = 1;
  o.use_k_vert = 1;
  o.eta_g = 1e-4;
  o.use_phi_s_ratio = 0;

  % Yield strength.
  o.ys.phi0 = 0.5;
  o.ys.phim = 0.75;
  o.ys.tauc = 1e6;     % constant that multiplies expression in phi_s
  o.ys.tauys = 1.5e6;  % yield strength of solid, bound on YS
  o.ys.gamma = 1;      % exponent
  
  % The nominal model:
%   o.model = 'bingham';
  % Alternative model. Yield strength is 0. However, friction remains
  % active. Mathematically, the original momentum balance equation
  %     vv F + vf (1 - F) = v
  % is modified to
  %     vv + vf = v.
  % The motivation for this alternative model is to get rid of all parameters
  % having to do with yield strength, and also to remove the imposed no-slip
  % boundary condition.
  o.model = 'newtonian';
  % Control plug_gas_loss behavior in the model = 'newtonian' case. Turn on gas
  % loss when vv/(vv + vf) <= vvfrac_thr.
%   o.newtonian.vvfrac_thr = 1e-3;
  o.newtonian.vvfrac_thr = 0.5;
  
  % Magma rheology
  o.eta.melt = 2;  %1=Hess&Dingwell; 2=Whittington
  o.eta.use_strain = 1;  %1 = using strain dependent rheology

  % Parameters for friction
  o.fr.v_ref = 0.05; % units m/s
  o.fr.f0 = 0.4;  % dimensionless
  o.fr.a  = 0.1;  % dimensionless
end

function opt_check (o, ot, str)
% Check o to make sure it doesn't have any extraneous or missing fields. Issues
% error on failure.
  if (nargin == 1)
    ot = opt_setdef();
    str = 'Extraneous field o';
  end
  fns = fieldnames(o);
  for (i = 1:numel(fns))
    if (~isfield(ot, fns{i}))
      error(sprintf('%s.%s\n', str, fns{i}));
    end
    if (isstruct(o.(fns{i})))
      opt_check(o.(fns{i}), ot.(fns{i}), sprintf('%s.%s', str, fns{i}));
    end
  end
  if (nargin == 1)
    opt_check(ot, o, 'Missing field o');
  end
end

function o = opt_convert (o1)
% slvof -> smf options struct. No longer needed, as we've moved away from
% solvof entirely.
  if (o1.testing == 1 || (o1.testing == 0 && o1.etac < 1))
    error('etac < 1 or testing not supported.');
  end
  if (o1.eta_factor ~= 1) error('eta_factor ~= 1 not supported'); end
  o = opt_setdef();
  flds = {'Tindex' 'v_init_exp' 'total_co2' 'total_h2o' 'k_lat_model', ...
          'plug_gas_loss' 'k_lat' 'L_diff' 'phi_gc' 'p_top'};
  for (i = 1:numel(flds)) o.(flds{i}) = o1.(flds{i}); end
  o.ys = o1.ys;
  o.k_lat_wall = o1.k_lat_wall;
  o.fr = struct('v_ref', o1.v_ref, 'f0', o1.f0, 'a', o1.a);
end

function m = m_init (o, v_ch)
% m (for model) is a struct that contains everything about the
% model. Initialize it based on the options struct and the BC v_ch.
  if (nargin < 2) v_ch = 1e-1; end
    
  m.ys = o.ys;
  m.fr = o.fr;
  
  m.p_top = o.p_top; % Pressure at surface (atmospheric) [Pa].
  
%   m.R = 50;          % [m]
%   m.R = 100;          % [m]
  m.R = o.R;          % [m]
  rho = 2500;        % [kg/m^3]
  rho_litho = 2700;  % [kg/m^3]
  rho_wat = 1000;    % [kg/m^3]
  m.g = 9.81;        % [m/s^2]
  m.conduit_length = o.conduit_length; % [m]
  m.ch = o.ch;
  
  % Normal stress. Add some surface parallel stress so > 0.
  m.sig.slope = (rho_litho - rho_wat)*m.g;
  m.sig.offset = 1e6;

  eta = 1e8;   % [Pa-s]
  m.expn = 1;  % exponent used in Bingham function
  
  Ttable = [950; 900; 850; 800; 750]; % [deg C]
  m.Ti = o.Tindex;
  m.T = Ttable(m.Ti) + 273.15;
  m.pf = pfl_init(m.Ti);
  
  m.chi_ch.total.h2o = o.total_h2o; % [percent]
  m.chi_ch.total.co2 = o.total_co2; % [ppm]
  
  m.Rw = 461.5000;
  m.Rc = 188.9000;
  m.B = 18.02/44.01;
  m.rho_l = 2200;
  m.rho_hd = m.rho_l;
  m.rho_cd = m.rho_l;
  m.rho_s = 2600;
  m.rho_c = 741;

  % Magma rheology
  m.eta.melt = o.eta.melt; 
  m.eta.use_strain = o.eta.use_strain; 
  
  m.k_lat_model = o.k_lat_model;
  m.plug_gas_loss = o.plug_gas_loss;
  m.k_lat = o.k_lat;
  m.eta_g = o.eta_g;
  m.use_phi_s_ratio = o.use_phi_s_ratio;
  m.phydro.slope = rho_wat*m.g;
  m.L_diff = m.R*o.L_diff;
  m.phi_gc = o.phi_gc;
  m.klw = klw_init(o.k_lat_wall);
  m.use_k_vert = o.use_k_vert;
  
  m.model = lower(o.model);
  if (~any(strcmp(m.model, {'bingham', 'newtonian'})))
    error('model must be one of ''bingham'', ''newtonian''');
  end
  m.newtonian = o.newtonian;
  
  % BCs
  if (isempty(o.p_ch))
    m.v_ch_guess = 10^o.v_init_exp; % [m/s]
    cc = 8*eta*m.v_ch_guess/m.R^2;
    m.p_ch = (rho*m.g + cc)*m.conduit_length; % [Pa]
  else
    m.p_ch = o.p_ch;
    cc = m.p_ch/m.conduit_length - rho*m.g;
    m.v_ch_guess = cc*m.R^2/(8*eta);
    if (m.v_ch_guess < 0)
      m.v_ch_guess = 1e-4;
    end
  end
  [m.mw_ch m.phi_g_ch] = solve_ch_bcs(m, o);
  
  % Info about the y vector.
  m.blk.is = struct('p', 1, 'v', 2, 'phi_g', 3, 'mw', 4);
  m.blk.lb = [0.01*m.p_top; 0; 0; 1e-3];
  m.blk.ub = [1.1*m.p_ch; inf; 1 - 1e-6; 1];
  % And the residual.
  m.blk.fs = struct('mb', 1, 'nv', 2, 'h2o', 3, 'co2', 4);
  
  % Scaling for numerical reasons.
  m.Nz = o.Nz;
  m.sy = [0.1*m.p_ch v_ch 0.1 1].';
  m.syp = -m.sy/m.conduit_length;
  m.sF = [m.p_ch;
          v_ch*[0.5*(m.rho_l + m.rho_s); m.rho_l; m.rho_l]]/m.conduit_length;
  
  % Solver options.
  m.slv = o.slv;
  
  if (m.slv.use_mex) feval(m.slv.mexFunc,'m_init', m); end
end

function m = m_rescale (m, v_ch)
% Scale y and yp.
  m.v_ch = v_ch;
  m.sy(m.blk.is.v) = v_ch;
  m.syp = -m.sy/m.conduit_length;
  
  if (m.slv.use_mex) feval(m.slv.mexFunc,'m_init', m); end
end

% sm: shooting method
%   The differential equations are solved from the chamber to the surface. At
% the chamber, there is one unspecified degree of freedom: v_ch. At the
% surface, there is the boundary condition p_top. Solve the scalar equation
%     fsm(v_ch) = p(z == 0) - p_top = 0
% for v_ch, thereby matching the surface p BC. Each evaluation of fsm(v_ch)
% requires solving the DAE, ie, integrating the equations from the chamber to
% the surface.

function [de flag msg] = solve (o, va, verbose, check_options)
% Caller-level function that solves the problem for a set of
% parameters. Final line of defense against bugs. The idea here is that MCMC,
% say, will never get interrupted by a coding bug in smf.m.
  if (nargin < 2) va = []; end
  if (nargin < 3) verbose = 0; end
  if (nargin < 4) check_options = 0; end
  
  % If requested, check the user's options struct.
  if (check_options)
    try
      opt_check(o);
    catch me
      disp_MException(me);
      msg = 'opt_check failed.';
      flag = -21;
      return;
    end
  end
  
  msg = '';
  try
    % Attempt to solve the problem.
    [de flag msg] = solve1(o, va, verbose);
    % If the first attempt didn't go well, try a higher tolerance on the ODE
    % integration, which sometimes helps.
    tol = 1e-6;
    if (any(flag == [-13 -10]) && o.slv.de.reltol > 5*tol)
      o.slv.de.reltol = tol;
      [de flag msg] = solve1(o, va, verbose);
    end
  catch me
    disp_MException(me);
    de.o = o;
    de.m = m_init(o);
    msg = 'Unplanned failure: probably a bug in smf.m.';
    flag = -20;
  end
end

function [de flag msg] = solve1 (o, va, verbose)
  if (nargin < 2) va = []; end
  if (nargin < 3) verbose = 0; end
  
  msg = '';
  % Print output to the command window if requested. In any case, record the
  % message.
  function pr (varargin)
    msg = sprintf(varargin{:});
    if (verbose) fprintf(msg); end
  end
  
  % Initialize the model.
  m = m_init(o);
  if (isempty(va)) va = m.v_ch_guess; end
  
  % Keep a global struct.
  clear global gsm;
  global gsm;
  gsm.f = inf; %
  gsm.de = []; % Reuse the two final solutions from the bracketing step.
  gsm.vi = []; %
  gsm.cnt = 0;       % Number of calls to the fsm(v_ch) function.
  gsm.tmr.t = tic(); % Timer to enforce a timeout (m.slv.max_time).
  
  de = [];
  
  % Bracket the solution if possible. This establishes a valid interval to
  % initialize fzero.
  [vi fi flag msg] = solve1_get_initial_interval(@pr, m, va);
  
  if (flag)
    pr('When bracketing: exiting with flag %d (%s).\n', flag, msg);
    return;
  end
  if (numel(vi) == 2)
    pr('Interval: [%1.2e %1.2e] with f [%1.2e %.1e2] %d evals\n', ...
       vi(1), vi(2), fi(1), fi(2), gsm.cnt);
  else
    pr('v_ch %1.2e %d evals\n', vi, gsm.cnt);
  end

  % Remember these f values so we don't recalc for fzero.
  gsm.vi = vi;
  gsm.fi = fi;

  % Now solve fsm(v_ch) = 0 using fzero.
  disp = 'off';
  if (verbose) disp = 'iter'; end
  fo = optimset('display', disp, 'tolx', m.slv.zero.tolx*min(vi), ...
                'funvalcheck', 'on');
  try
    % sm_calc_fsm implements fsm(v_ch) = p(z = 0) - p_top. Find v_ch so
    % fsm(v_ch) = 0 (to a tolerance).
    [v_ch f flag] = fzero(@(v_ch) sm_calc_fsm(m, v_ch), vi, fo);
  catch me
    disp_MException(me);
    de.o = o; de.m = m;
    if (toc(gsm.tmr.t) > m.slv.max_time)
      flag = -14;
      msg = 'time to solve > m.slv.max_time';
    else
      flag = -15;
      msg = '';
    end
    pr('flag %d; exiting early: %s\n', flag, msg);
    return;
  end
  
  % Gather output.
  if (v_ch == gsm.de.y(2, 1))
    de = gsm.de;
  else
    m.slv.max_time = inf;
    de = de_run(m, v_ch);
  end
  de.o = o;
  de.m = m;
  
  % Check some additional exit conditions.
  msg = 'See ''help fzero'' for interpretation of this flag value.';
  if (flag == 1)
    if (de.z(end) ~= 0)
      flag = -13;
      msg = 'Did not and probably can not integrate to surface.';
    elseif (abs(de.y(1,end) - de.m.p_top)/de.m.p_top > 1.5)
      msg = 'Surface p boundary condition not matched very well.';
      flag = -10;
    end
  end
  
  % Done.
  pr('flag %d (%s);\nfound v_ch %1.5e with p_top %1.2e after %d evals\n', ...
     flag, msg, v_ch, de.y(m.blk.is.p, end), gsm.cnt);
  clear global gsm;
end

function [vi fi flag msg] = solve1_get_initial_interval (pr, m, va)
% Get an initial interval. I think I can do it faster than fzero b/c I know
% more. (Experiments show that indeed the following procedure typically uses 2-4
% evaluations vs fzero's many more.)
  vi = []; % v_ch at each end of the interval.
  fi = []; % fsm(v_ch) at the values in vi.
  flag = 0; msg = '';
  fa = sm_calc_fsm(m, va);
  if (isnan(fa))
    % Bad user guess. We take matters into our own hands.
    pr('Input guess gives nan.\n');
    va = 0.1;
    while (va >= 1e-6)
      fa = sm_calc_fsm(m, va);
      if (~isnan(fa)) break; end
      va = va/10;
    end
    if (isnan(fa))
      pr('0.1 to 1e-6 gives nan.\n');
      va = 0.1;
      while (va <= 1e2)
        fa = sm_calc_fsm(m, va);
        if (~isnan(fa)) break; end
        va = 2*va;
      end
    end
    if (isnan(fa))
      msg = 'Failed to find va.';
      pr([msg '\n']);
      de = [];
      flag = -11;
      return;
    end
  end
  % Now vb.
  if (fa < 0)
    vb = va/10;
    while (vb < va)
      fb = sm_calc_fsm(m, vb);
      if (isnan(fb)) vb = 2*vb; else break; end
    end
    if (isnan(fb))
      vb = [];
    elseif (sign(fb) == sign(fa))
      fac = 10;
      for (i = 1:10)
        vb = vb/fac;
        fb = sm_calc_fsm(m, vb);
        if (isnan(fb)) break; end
        if (sign(fb) ~= sign(fa)) break; end
        fac = fac*2;
      end
      if (isnan(fb)) vb = []; end
      if (sign(fb) == sign(fa))
        msg = sprintf(['Failed to find vb because fb is decreasing; ', ...
                       'smallest vb is %1.2e.'], vb);
        pr([msg '\n']);
        flag = -12;
%         
        return;
      end
    end
    vi = [vb va];
    if (isnan(fb)) fb = []; end
    fi = [fb fa];
  else
    vb = va;
    fac = 10;
    while (1)
      vb = fac*vb;
      fb = sm_calc_fsm(m, vb);
      if (isnan(fb)) break; end
      if (sign(fa) ~= sign(fb)) break; end
    end
    if (isnan(fb)) vb = []; fb = []; end
    vi = [va vb];
    fi = [fa fb];
  end
  
end

function de = add_fields (de)
% Add some fields useful for analysis, but not actually solving the equation,
% to de.
  if (~isfield(de, 'z') || isempty(de.z)) return; end
        
  
  flds = {};
  for (i = 1:numel(de.z))
    state = struct('k', de.state(1,i), 'plug', de.state(2,i));
%     if (de.m.plug_gas_loss), state.Kplug = de.Kplug; end

    y = put_in_bounds(de.m, de.y(:,i));
    [~, c] = des_calc_exprs(de.m, state, de.z(i), y, [], 1:3);
    if (isempty(flds))
      flds = setdiff( ...
        fieldnames(c), ...
        {'y' 'yp' 'yp_final' 'z' 'eqs' 'state' 'm' 'c' 'db' 'args' 'pp0'});
      for (j = 1:numel(flds))
        de.(flds{j}) = zeros(numel(de.z), 1);
      end
    end
    for (j = 1:numel(flds))
      if (isfield(c, flds{j}))
        de.(flds{j})(i) = c.(flds{j});
      else
        de.(flds{j})(i) = nan;
      end
    end
  end
  de.z = de.z.';
  
  for (i = 1:size(de.y,1))
    is = find(de.y(i,:) < de.m.blk.lb(i));
    if (any(is ~= numel(de.z)))
      warning(sprintf('Outside of bound in field %d at index %d\n', ...
                      i, is(1)));
    end
  end
end

function de = add_measurement_fields (de)
% Like add_fields, but for MCMC, add just what is needed using a mex
% implementation. If the mex file is not available, call add_fields.
  if (de.m.slv.use_mex)
    % Speedup.
    feval(m.slv.mexFunc,'m_init', de.m);
    de = feval(m.slv.mexFunc,'add_measurement_fields', de);
  else
    % Go with the slow approach.
    de = add_fields(de);
  end
end

function [f flag de] = sm_calc_fsm (m, v_ch)
% Calculate fsm(v_ch).
  flag = 0;
  
  % See if we can reuse the values from the bracketing step.
  global gsm;
  if (~isempty(gsm.vi))
    i = find(v_ch == gsm.vi);
    if (isempty(i))
      gsm.vi = [];
    else
      f = gsm.fi(i);
      return;
    end
  end
  gsm.cnt = gsm.cnt + 1;
  
  % Integrate the DAE.
  try
    [de flag] = de_run(m, v_ch);
  catch me
    disp_MException(me);
    f = nan;
    flag = -3;
    return;
  end
  if (flag ~= 1)
    if (flag == -2)
      % Didn't get the initial condition.
      p_top = nan;
    else
      % Extrapolate p to the top.
      p_top = de.y(1,end) - de.yp_final(1)*de.z(end);
    end
  else
    p_top = de.y(1,end);
  end
  
  % We succeeded in the integration, so now calculate fsm(v_ch):
  f = p_top - m.p_top;
  
%   semilogx(v_ch, f, '+'); hold on; drawnow;
  
  % Record this solution in the global struct.
  global gsm;
  if (abs(f) < gsm.f)
    gsm.f = abs(f);
    gsm.de = de;
  end
end

% de: differential equation. Solve F(y, yp, z) = 0. See des_calc_F for more.
function [de flag] = de_run (m, v_ch)
% Top level routine to integrate the DAE from the chamber to the surface. The
% DAE is implemented in des_calc_F.
  
  % Set initial (in the context of the DAE; boundary in the context of the PDE)
  % conditions at the chamber.
  m = m_rescale(m, v_ch);
  m.v_ch = v_ch;
  y0 = [m.p_ch; v_ch; m.phi_g_ch; m.mw_ch];
  yp0 = zeros(4, 1);
  yp0(m.blk.is.p) = y0(m.blk.is.p)/m.conduit_length;
  
  % Set state variables.
  if (is_thr_perm(m.k_lat_model) && m.phi_g_ch < m.phi_gc)
    state.k = 0;
  else
    state.k = 1;
  end
  state.plug = 0;
  
  dz = 1/(m.Nz-1);
  de.z = []; de.y = []; de.yp = []; de.state = [];
  x0 = 0;
  flag = 1;
  first = 1;
  % Integrate over segments, where each segment is determined by state
  % variables. The segments are not known in advance. Each segment involves a
  % separate call to ode15i.
  while (1)
    % Initialize the ode solver and data for this segment.
    [des flag] = des_init(m, state, x0, y0, yp0);
    state = des.state;
    if (flag ~= 1)
      if (first) flag = -2; else flag = -1; end
      break;
    end
    % Integrate this segment.
    try
      od = des_run(m, des);
    catch
      % Timed out.
      if (first) flag = -2; else flag = -1; end
      break;
    end
    first = 0;
    
    % Record data for this segment.
    if (numel(od.x) == 1)
      % Difficulties.
      assert(isempty(od.ie));
      flag = -1;
      break;
    end
    
    if od.x(1) == x0
        zi = ceil(x0/dz)*dz:dz:od.x(end);
    else
        zi = (de.z(end)+dz):dz:od.x(end);
    end
    
    if ~isempty(zi)
        [yi,ypi] = deval(od,zi);
        de.z = [de.z zi]; de.y = [de.y yi]; de.yp = [de.yp ypi];
        de.state = [de.state repmat([state.k; state.plug], 1, numel(zi))];
        de.yp_final = od.extdata.ypfinal;
    %db_plot(m, od, de, des);
    end
    
    % Made it to the surface? If so, we're done.
    if (od.x(end) == 1) break; end
    
    % No, so prepare next segment.
    ye = od.y(:,end);
    x0 = od.x(end);
    y0 = od.y(:,end).*m.sy;
    yp0 = de.yp_final.*m.syp;
    
    % Determine why the previous segment ended and change state accordingly.
    if (isempty(od.ie))
      % Integration failure.
      flag = -1;
      break;
    end
    for (i = 1:numel(od.ie))
      switch (od.ie(i))
        case 1
          state.k = 1;
        case 2
          % Pressure reached 0. Can't proceed, so quit.
          flag = -1;
          break;
        case 3
%           zplug  = des_transform (m, x0);
%           c = des_calc_exprs(m, state, zplug, y0, yp0, 1:3);
%           state.Kplug = c.rho_g*(c.mbe.vv + c.mbe.vf + c.vgdiff);
%           m.rho_g_plug = db.rho_g;
%           m.Kplug.K = (c.h2o.f + c.co2.f)./y0(m.blk.is.phi_g);
%           rho_g_top = m.p_top/m.Rw/m.T;
%           m.Kplug.ktop = 1e-9;
%           m.Kplug.zc = -50;
%           m.Kplug.ktop = m.eta_g/db.pp*(m.Kplug.K/rho_g_top - db.v);
%           m.Kplug.zc = zplug/log(db.k_vert/m.Kplug.ktop);
          state.plug = mod(state.plug + 1, 2);
        otherwise
          error('Not an event ie.');
      end
    end
    if (flag < 0) break; end
  end

  % Scale to science units.
  if (~isempty(de.z))
    de.z = m.conduit_length*(1 - de.z);
    de.y = repmat(m.sy, 1, size(de.y, 2)).*de.y;
    de.yp_final = de.yp_final.*m.syp; 
    de.yp = repmat(m.syp,1,size(de.yp,2)).*de.yp;
%     if m.plug_gas_loss, de.Kplug = state.Kplug; end
  end
end

function db_plot (m, od, de, des)
% Debug plot.
  plot(de.z, de.y, '.-'); drawnow;
  yl = ylim();
  if (yl(2) > 10) ylim([0 10]); end
  if (0)
    global gdb;
    hold on;
    x = linspace(0, 1, gdb.s.m.nsteps);
    plot(x, mop([gdb.s.db.p gdb.s.db.velocity gdb.s.db.phi_g gdb.s.db.mw],...
                '/', m.sy), '--');
    hold off; drawnow;
  end
  hold on; plot(od.x(1), od.y(:,1), 'ko'); drawnow;
end

% des: de segment.
function [des flag] = des_init (m, state, x0, y0, yp0)
% Initialize ode15i, and get consistent initial conditions, for this segment.
  des.x0 = x0;
  des.state = state;
  des.o = odeset('reltol', m.slv.de.reltol, 'abstol', m.slv.de.abstol);
  %des.o.OutputFcn = @odeplot;
  
  des.o.Events = @(x, y, yp) des_calc_event(m, des.state, x, y, yp);
  des.o.Jacobian = @(x, y, yp) des_calc_F_jac(m, des.state, x, y, yp);

  z = des_transform(m, x0);
  c = des_calc_exprs(m, state, z, y0, [], 2:3);
  if (m.plug_gas_loss && des.state.plug == 0 && c.F == 0)
    des.state.plug = 1;
    c = des_calc_exprs(m, state, z, y0, [], 2:3);
  end
  f = [c.nv.f; c.h2o.f; c.co2.f];
  g = [c.nv.g; c.h2o.g; c.co2.g];
  f_y = des_calc_f_y(m, state, z, y0, [], f);
  M = [1 0 0 0; f_y];
  rhs = [c.pp; g];

  wid = 'MATLAB:nearlySingularMatrix';
  ws = warning('query', wid);
  warning('off', wid);
  yp = M \ rhs;
  warning(ws.state, wid);
  
  des.y0 = y0./m.sy;
  des.yp0 = yp./m.syp;
  flag = 1;
end

function od = des_run (m, des)
% Wrapper around ode15i. Handle warning spam.
  wid = 'MATLAB:ode15i:IntegrationTolNotMet';
  % Turn off the warning so stdout doesn't get spammed. We'll detect failure
  % in the outer loop.
  ws = warning('query', wid);
  warning('off', wid);  
  od = ode15i(@(x, y, yp) des_calc_F(m, des.state, x, y, yp), ...
              [des.x0 1], des.y0, des.yp0, des.o);
  warning(ws.state, wid);
end

function [z y yp] = des_transform (m, x, y, yp)
% Transform from math (scaling for numerical reasons) to science units.
  z = m.conduit_length*(1 - x);
  if (nargout == 1) return; end
  y = y.*m.sy;
  if (nargout > 2 && ~isempty(yp)) yp = yp.*m.syp; end
end

function a = out_of_bounds (m, y)
  a = any(y < m.blk.lb | y > m.blk.ub);
end

function a = not_good (y)
  a = any(isinf(y) | isnan(y));
end

function y = put_in_bounds (m, y)
  mask = y < m.blk.lb;
  y(mask) = m.blk.lb(mask);
  mask = y > m.blk.ub;
  y(mask) = m.blk.ub(mask);
end

function F = des_calc_F (m, state, x, y, yp)
% Function passed to ode15i.
  if (m.slv.use_mex)
    F = feval(m.slv.mexFunc,'calc_F', state, x, y, yp);
    return;
  end
  % Scale and keep variables numerically in bounds.
  [z y yp] = des_transform(m, x, y, yp);
  y = put_in_bounds(m, y);
  % The system F(y, yp, z) = 0 is constructed as follows.
  %   The momentum balance equation is solved for dp/dz, and then the first
  % equation in F = 0 is
  %     dp/dz - yp(m.blk.is.p) = 0.
  % It gives the difference between dp/dz that is consistent with the m.b.e. and
  % the value that the ODE integrator currently has.
  %   The continuity equations have the form
  %     d f(z, y)/dz - g(y) = 0,
  % where f, g \in Real^3 since there are three continuity equations. To turn
  % this equation into an implicit ODE in dy/dz, take the partial derivative
  % of f to get
  %     f_y(z, y) dy/dz + f_z(z, y) - g(y) = 0.
  c = des_calc_exprs(m, state, z, y, yp, 2:3);
  % f(z, y) in the above.
  f = [c.nv.f; c.h2o.f; c.co2.f];
  % g(y) in the above.
  g = [c.nv.g; c.h2o.g; c.co2.g];
  % f_y(z, y) in the above.
  f_y = des_calc_f_y(m, state, z, y, yp, f);
  % f_z(z, y) in the above.
  f_z = des_calc_f_z(m, state, z, y, yp, f);
  % Assemble F(y, yp, z).
  F = [yp(m.blk.is.p) - c.pp;
       f_y*yp + f_z - g];
  % Scale.
  F = F./m.sF;
end

function [F_y F_yp] = des_calc_F_jac (m, state, x, y, yp)
% Jacobian function (of F, not f) passed to ode15i through odeset. Use a
% combination of finite differences and problem structure. Read 'help odeset'
% for usage.
  global gsm;
  if (~isempty(gsm) && toc(gsm.tmr.t) > m.slv.max_time) error('timeout'); end
  [z y yp] = des_transform(m, x, y, yp);
  % Put back in bounds. I suppose that's ok for the Jacobian.
  y = put_in_bounds(m, y);
  
  c = des_calc_exprs(m, state, z, y, yp, 2:3);
  h = -c.pp; % Drop the yp(m.blk.is.p) part b/c it's constant.
  f = [c.nv.f; c.h2o.f; c.co2.f];
  g = [c.nv.g; c.h2o.g; c.co2.g];
  f_y = des_calc_f_y(m, state, z, y, yp, f);
  
  h_y = zeros(1, 4);
  g_y = zeros(3, 4);
  % A(:,j) = partial( f_y(z, y) dy/dz, y(j) ) = partial( f_y(z, y), y(j) ) dy/dz.
  A = zeros(3, 4);
  % B(:,j) = partial( f_z(z, y), y )
  B = zeros(3, 4);
  for (i = 1:4)
    delta = m.slv.fd*m.sy(i);
    yd = y;
    yd(i) = y(i) + delta;
    on_bound = yd(i) > m.blk.ub(i);
    if (on_bound)
      delta = -delta;
      yd(i) = y(i) + delta;
    end
    c = des_calc_exprs(m, state, z, yd, yp, 2:3);
    hd = -c.pp;
    fd = [c.nv.f; c.h2o.f; c.co2.f];
    gd = [c.nv.g; c.h2o.g; c.co2.g];
    done = 0;
    if (m.slv.fdo == 2 && ~on_bound)
      ydm = y;
      ydm(i) = y(i) - delta;
      if (ydm(i) >= m.blk.lb(i))
        cm = des_calc_exprs(m, state, z, ydm, yp, 2:3);
        hdm = -cm.pp;
        fdm = [cm.nv.f; cm.h2o.f; cm.co2.f];
        gdm = [cm.nv.g; cm.h2o.g; cm.co2.g];
        delta = 2*delta;
        h_y(i) = (hd - hdm)/delta;
        A(:,i) = (des_calc_f_y(m, state, z, yd, yp, fd) - ...
                  des_calc_f_y(m, state, z, ydm, yp, fdm))*yp/delta;
        B(:,i) = (des_calc_f_z(m, state, z, yd, yp, fd) - ...
                  des_calc_f_z(m, state, z, ydm, yp, fdm))/delta;
        g_y(:,i) = (gd - gdm)/delta;
        done = 1;
      end
    end
    if (~done)
      h_y(i) = (hd - h)/delta;
      A(:,i) = (des_calc_f_y(m, state, z, yd, yp, fd) - f_y)*yp/delta;
      g_y(:,i) = (gd - g)/delta;
    end
  end  
  F_y = [h_y; A + B - g_y];
  
  h_yp = zeros(1, 4);
  h_yp(m.blk.is.p) = 1;
  F_yp = [h_yp; f_y];
  
  S = diag(1./m.sF);
  F_y = S*F_y*diag(m.sy);
  F_yp = S*F_yp*diag(m.syp);
  
  if (0)
    % Do some code checks.
    [t.F_y t.F_yp] = des_calc_F_jac_fd(m, state, x, y./m.sy, yp./m.syp);
    if (1 || relerr(t.F_y, F_y) > 1e-3 || relerr(t.F_yp, F_yp) > 1e-3)
      pr('%e %e\n', cond(F_y), cond(F_yp));
      pr('  %e %e\n', relerr(t.F_y, F_y), relerr(t.F_yp, F_yp));
    end
    if (isinf(cond(F_y)) || isinf(cond(F_yp))) kb; end
    %F_y = t.F_y; F_yp = t.F_yp;
  end
end

function [F_y F_yp] = des_calc_F_jac_fd (m, state, x, y, yp)
% Check the more clever routine des_calc_F_jac above. Not used in practice.
  F = des_calc_F(m, state, x, y, yp);
  F_y = zeros(4); F_yp = zeros(4);
  for (i = 1:4)
    delta = m.slv.fd;
    yd = y;
    yd(i) = y(i) + delta;
    if (yd(i)*m.sy(i) > m.blk.ub(i))
      delta = -delta;
      yd(i) = y(i) + delta;
    end
    F_y(:,i) = (des_calc_F(m, state, x, yd, yp) - F)/delta;
    delta = m.slv.fd;
    ypd = yp;
    ypd(i) = yp(i) + delta;
    F_yp(:,i) = (des_calc_F(m, state, x, y, ypd) - F)/delta;
  end
end

function [val isterm dir] = des_calc_event (m, state, x, y, yp)
% Event function passed to ode15i through odeset. Read 'help odeset' and
% 'help od15i' for usage.
  if (nargin < 5) yp = []; end
  val = ones(3, 1); isterm = val; dir = val;
  [z y yp] = des_transform(m, x, y, yp);

  if (is_thr_perm(m.k_lat_model) && state.k == 0)
    val(1) = y(m.blk.is.phi_g) - m.phi_gc;
    dir(1) = 1;
  end

  % Make sure p >= 0. This is treated as an event so the integrator quits if
  % pressure goes negative.
  val(2) = y(m.blk.is.p) - m.blk.lb(m.blk.is.p);
  dir(2) = -1;

  if (m.plug_gas_loss)
    y = put_in_bounds(m, y);
    switch (m.model)
      case 'bingham'
        c = des_calc_exprs(m, state, z, y, yp, []);
        if (state.plug == 0)
          % When going from F > 0 to F == 0, switch based on F because F is set
          % to 0, while zeta is calculated.
          val(3) = c.F;
          dir(3) = -1;
        else
          % But when going in the opposite direction, we want to detect the
          % switch based on zeta going from > 1 to 1. F == 0 is uninformative
          % because it is set.
          val(3) = c.zeta - 1;
          dir(3) = -1;
        end
      case 'newtonian'
        [c] = des_calc_exprs(m, state, z, y, yp, 1);
        %nb If the definition of vvfrac changes here, change also the code in calc_exprs
        % and smfmex.cpp.
        val(3) = c.mbe.vv/(c.mbe.vv + c.mbe.vf) - m.newtonian.vvfrac_thr;
        if (state.plug == 0)
          dir(3) = -1;
        else
          dir(3) = 1;
        end
    end
  end
end

function [c, db] = des_calc_exprs (m, state, z, y, yp, eqs)
% The science is implemented here. To some extent, you can modify the
% expressions here and in routines called from here without worrying about the
% higher-level operations being done in the solver. One obvious exception to
% this is if you introduce new physics that causes additional nonsmoothness in
% the solution. There may be other exceptions.
%   NB that yp is used only in the following way: yp(m.blk.is.p) is used as
% an initial guess when call mb_solve. In particular, none of these
% expressions is a function of yp. Consequently, we have f(y) and not f(y,
% dy/dz).
  if (m.slv.use_mex && nargout == 1)
    c = feval(m.slv.mexFunc,'calc_exprs', state, z, y, yp, eqs);
    return;
  end
  p = y(m.blk.is.p);
  v = y(m.blk.is.v);
  phi_g = y(m.blk.is.phi_g);
  mw = y(m.blk.is.mw);
  
  % Solubility ratio.
  Gamma = (1 - mw)/(mw*m.B);
  
  % Solid volume fraction.
  % assuming volatile saturation
  [Chi_hd Chi_cd] = solubility_liu_explicit(p, m.T, mw);
  Chi_hd = 1e-2*Chi_hd; Chi_cd = 1e-6*Chi_cd;
  c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
  c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
  c12 = c1*c2;
  phi_s_of_p = pfl_val(m.pf, p*1e-6);
  phi_s = phi_s_of_p*m.rho_l*c12*(1 - phi_g)/...
            (m.rho_s + phi_s_of_p*(c12*m.rho_l - m.rho_s));

  phi_l = (1 - phi_s - phi_g)*c2; % liquid volume fraction
  rho_g = p*(mw/(m.Rw*m.T) + (1 - mw)/(m.Rc*m.T)); % gas density
  rho = m.rho_l*phi_l*c1 + m.rho_s*phi_s + rho_g*phi_g; % bulk density
  
  % Solve the momentum balance equation for dp/dz.
  if (m.use_phi_s_ratio)
    phi_s_eta = phi_s/(1 - phi_g); %editted Oct 2, 2017
  else
    phi_s_eta = phi_s;
  end

  if (m.eta.use_strain==1)
      gdot = 2*v/m.R;
  elseif (m.eta.use_strain==0)
      gdot = 0;
  else 
      gdot = m.eta.use_strain;
  end
  
  eta = calc_eta(Chi_hd, phi_s_eta, m.T, m.eta.melt, gdot);
  Visc = 0.25*m.R/eta;
  Psi = 2*m.fr.v_ref/exp(m.fr.f0/m.fr.a);
  Den = m.fr.a*(m.sig.slope*z + m.sig.offset);
  
  pp0 = [];
  if (~isempty(yp)) pp0 = yp(m.blk.is.p); end
  tolx = 1e2*eps;
  tolf = 1e2*eps*v;
  switch (m.model)
    case 'bingham'
      tauY = calc_yield(m.ys, phi_s, phi_l);
      [c.pp f c.F c.zeta] = mb_solve( ...
        -0.5*m.R*rho*m.g, 0.5*m.R, tauY, Visc, Psi, 1/Den, -v, m.expn, tolx, ...
        tolf, pp0);
     case 'newtonian'
       [c.pp f c.F c.zeta] = mbv_solve( ...
         -0.5*m.R*rho*m.g, 0.5*m.R, Visc, Psi, 1/Den, -v, tolx, tolf, pp0);
       tauY = 0;
  end
  
  for (i = 1:numel(eqs))
    switch (eqs(i))
      case 1 % Momentum balance.
        % Since the m.b. equation is solved for dp/dz, it is satisfied within
        % the current context by construction. However, we may use the key
        % quantities from the m.b.e. in the event function, so I provide
        % these here.
        tauR = -0.5*m.R*(-c.pp + rho*m.g);
        switch (m.model)
          case 'bingham'
            c.mbe.vv = tauR*Visc*c.F;
            c.mbe.vf = Psi*sinh(tauR/Den)*(1 - c.F)^m.expn;
          case 'newtonian'
            c.mbe.vv = tauR*Visc;
            c.mbe.vf = Psi*sinh(tauR/Den);
        end

      case 2 % Non-volatile continuity; RHS only.
        c.nv.f = (m.rho_l*phi_l + m.rho_s*phi_s)*v;
        c.nv.g = 0;
        
      case 3 % Water and CO2 continuity; RHS only.
        C_fac = Gamma/(1 + Gamma);
        
        if (state.plug == 1)
          c.h2o.f = phi_g/(1 + Gamma)*v;
          c.h2o.g = 0;

          c.co2.f = phi_g*C_fac*v;
          c.co2.g = 0;
 
%           if isempty(yp)
%               drhogdz = rho_g/p*(c.pp);
%           else
%               drhogdz = rho_g/p*(c.pp) + p*(1/(m.Rw*m.T) - 1/(m.Rc*m.T))*(-yp(m.blk.is.mw));
%           end
%                      
%           c.h2o.f = (Chi_hd*m.rho_l*phi_l*c1 + phi_g*rho_g/(1 + Gamma))*v;
%           c.h2o.g = 1/(1+Gamma)*phi_g*v*drhogdz;
%           
%           c.co2.f = (Chi_cd*m.rho_l*phi_l*c1 + phi_g*rho_g*C_fac)*v;
%           c.co2.g = C_fac*phi_g*v*drhogdz;

        else
          
          % Differential vertical gas velocity.
          k_vert = 0;
          if (m.use_k_vert && state.k == 1)
            k_vert = calc_k(m, phi_g, state);
            % c.pp should be numerically equal to, but is not formally, dp/dz.
            % This distinction is important. We still have d f(y)/dz = g(z, y)
            % rather than d f(y, dy/dz)/dz = g(z, y).
            vgdiff = (k_vert/m.eta_g)*c.pp;
%             if state.plug == 1
%                 vgdiff = state.Kplug/rho_g - v;
%             end
          else
            vgdiff = 0;
          end
          
          k_lat = calc_k_lat(m, z, phi_g, state);
%           k_lat = calc_k_lat(m, z, phi_l, phi_g, c.pp, rho_g, v, state);
          u_lat = k_lat/m.eta_g*(p - m.phydro.slope*z)/(m.R+m.klw.L);

          % added c1 Oct 16, 2017
          water_mass = Chi_hd*m.rho_l*phi_l*c1 + phi_g*rho_g/(1 + Gamma); 
          water_vert = phi_g*rho_g*vgdiff/(1 + Gamma);
          water_loss = 2/m.R*rho_g*phi_g*u_lat/(1 + Gamma);
          
          c.h2o.f = water_mass*v + water_vert;
          c.h2o.g = water_loss;
          
          C_mass = Chi_cd*m.rho_l*phi_l*c1 + phi_g*rho_g*C_fac;
          C_vert = phi_g*rho_g*vgdiff*C_fac;
          C_loss = 2/m.R*rho_g*phi_g*u_lat*C_fac;
          
          c.co2.f = C_mass*v + C_vert;
          c.co2.g = C_loss;
          
          c.vgdiff = vgdiff;
          c.rho_g = rho_g;
          
%           if state.plug == 1 
%               c.h2o.g = 0; c.co2.g = 0;
%           end
        end

      otherwise
        error('eq must be in {1, 2, 3}');
    end
  end

  if (nargout > 1)
    pp = c.pp;
    F = c.F; zeta = c.zeta;
    vv = c.mbe.vv; vf = c.mbe.vf;
    % Compute vvfrac regardless of model. No reason not to in this slow Matlab
    % version.
    vvfrac = vv/(vv + vf);
    db = [];
    vars = who();
    for (i = 1:numel(vars))
      eval(sprintf('db.%s = %s;', vars{i}, vars{i}));
    end
  end
end

function f_y = des_calc_f_y (m, state, z, y, yp, f)
% Calculate f(y)_y using finite differences.
  if (m.slv.use_mex)
    f_y = feval(m.slv.mexFunc,'calc_f_y', state, z, y, yp, f);
    return;
  end
  f_y = zeros(3, 4);
  yp = [];
  for (i = 1:4)
    delta = m.slv.fd*m.sy(i);
    yd = y;
    yd(i) = y(i) + delta;
    on_bound = yd(i) > m.blk.ub(i);
    if (on_bound)
      delta = -delta;
      yd(i) = y(i) + delta;
    end
    c = des_calc_exprs(m, state, z, yd, yp, 2:3);
    if (i == 1) yp = [c.pp; zeros(3, 1)]; end
    fd = [c.nv.f; c.h2o.f; c.co2.f];
    done = 0;
    if (m.slv.fdo == 2 && ~on_bound)
      yd(i) = y(i) - delta;
      if (yd(i) >= m.blk.lb(i))
        c = des_calc_exprs(m, state, z, yd, yp, 2:3);
        fdm = [c.nv.f; c.h2o.f; c.co2.f];
        f_y(:,i) = (fd - fdm)/(2*delta);
        done = 1;
      end
    end
    if (~done)
      f_y(:,i) = (fd - f)/delta;
    end
  end
end

function f_z = des_calc_f_z (m, state, z, y, yp, f)
% Calculate f(z, y)_z using finite difference - central difference.
  f_z = zeros(3, 1);
  delta = m.slv.fd*m.conduit_length;
  zd = z + delta;
  on_bound = zd > m.conduit_length;
  if (on_bound)
    delta = -delta;
    zd = z + delta;
  end
  yp = [];
  c = des_calc_exprs(m, state, zd, y, yp, 2:3);
  fd = [c.nv.f; c.h2o.f; c.co2.f];
  done = 0;
  if (m.slv.fdo == 2 && ~on_bound)
    zd = z - delta;
    if (zd >= 0)
      yp = [c.pp; zeros(3, 1)];
      c = des_calc_exprs(m, state, zd, y, yp, 2:3);
      fdm = [c.nv.f; c.h2o.f; c.co2.f];
      f_z = (fd - fdm)/(2*delta);
      done = 1;
    end
  end
  if (~done)
    f_z = (fd - f)/delta;
  end
end

function disp_MException (me)
  fprintf('%s\n', me.message);
  for (i = 1:numel(me.stack))
    fprintf('%s (%d): %s\n', me.stack(i).file, me.stack(i).line, ...
            me.stack(i).name);
  end
end

function pf = pfl_init (Ti)
% Interp/fit the chi-pressure curves with a local representation. This uses
% the data in xtl_tables.mat, which were obtained by tracing the curves in
% Fig A4 of the Schneider et al SI. A somewhat hand-massaged, ad hoc
% procedure is used to come up with satisfactory curves for use in the
% solving the DAE. The objective is to come up with curves that (a) agree
% with Fig A4 to essentially a visual level and (b) are as smooth as I can
% make them to make the DAE integration faster and easier.
  global gpfl;
  if (~isfield(gpfl, 'pfs'))
    lo = load('xtl_tables.mat');
    gpfl.pfs = cell(1, numel(lo.ts));
    chi_thrs = linspace(0.6, 0.7, 5);
    for (i = 1:numel(lo.ts))
      p = flipud(lo.ts{i}(:,2));
      chi = flipud(lo.ts{i}(:,1));
      gpfl.pfs{i} = pfl_fit(p, chi, i == 1, chi_thrs(i));
    end
  end
  pf = gpfl.pfs{Ti};
end

function pf = pfl_fit (p, chi, fit_line, chi_thr)
  % I see that the ends are not well done, and in the paper, the ends are
  % basically lines.
  %   For the T = 950 curve, we need to fit a line on one end.
  %   For all other ends, extrap.
  %   For all curves, fit a line for chi >= chi_thr and chi < ~0.95.
  mask = chi >= 0.95;
  if (fit_line)
    n = 6;
    chi(end-n+1:end) = pfl_fit_line(p(end-n+1:end), chi(end-n+1:end));
  else
    mask = mask | p >= 190;
  end
  p(mask) = [];
  chi(mask) = [];
  mask = chi >= chi_thr;
  chi(mask) = pfl_fit_line(p(mask), chi(mask));

  is = pfl_select_is(p, chi, 7e-4);
  pf.p = p(is);
  pf.chi = chi(is);
  
  % Linearly extrapolate in each direction.
  chi0 = 1;
  pe = 200;
  a = (pf.p(1) - pf.p(2))/(pf.chi(1) - pf.chi(2));
  p0 = pf.p(1) + a*(chi0 - pf.chi(1));
  a = (pf.chi(end) - pf.chi(end-1))/(pf.p(end) - pf.p(end-1));
  chie = pf.chi(end) + a*(pe - pf.p(end));
%   pf.p = [p0; pf.p(2:end-1); pe];
%   pf.chi = [chi0; pf.chi(2:end-1); chie];
  
%   pf.pp = interp1(pf.p, pf.chi, 'pchip', 'pp');

  % editted Jan 22, 2018
  pf.p = [0.1; p0; pf.p(2:end-1); pe];
  pf.chi = [chi0; chi0; pf.chi(2:end-1); chie];
  pf.pp = pchip(pf.p, pf.chi);
end

function chi = pfl_fit_line (p, chi)
  n = numel(p);
  A = [ones(n, 1), p(:)];
  c = A \ chi(:);
  chi = c(1) + c(2)*p;
end

function is = pfl_select_is (p, chi, dr)
  n = numel(p);
  % Fit lines to 3-point pieces and measure residual. This measures
  % curvature.
  A = ones(3, 2);
  r = zeros(n-2, 1);
  for (i = 2:n-1)
    iis = i-1:i+1;
    A(:,2) = p(iis);
    c = A \ chi(iis);
    r(i) = norm(chi(iis) - (c(1) + c(2)*p(iis)), 2);
  end
  r = cumsum(r);
  r = [0; r; r(end)];
  % Choose is to be ~ equally spaced in r.
  is = 1;
  k = 2;
  while (k < n)
    i = find(r(k:end) - r(k-1) >= dr, 1);
    if (isempty(i)) break; end
    is(end+1) = k + i - 1;
    k = k + i;
  end
  is(end+1) = n;
  is = unique(is);
end

function chi = pfl_val (pf, p)
  chi = zeros(size(p));
  chi(p <= pf.p(1)) = pf.chi(1);
  chi(p >= pf.p(end)) = pf.chi(end);
  mask = p > pf.p(1) & p < pf.p(end);
  chi(mask) = ppval(pf.pp, p(mask));
  assert(all(chi <= 1));
end

function [x f F zeta nit] = mb_solve (a, b, c, d, g, h, k, p, tolx, tolf, x0)
% Solve the momentum balance equation for dp/dz (here x):
%     zeta(x) = c/(a + b x);
%     F(x) = 0 if zeta(x) >= 1
%            1 if zeta(x) <= 0
%            1 + (zeta^4 - 4 zeta)/3 else;
%     d (a + b x) F(x) + g sinh(h (a + b x)) (1 - F(x))^p + k = 0.
% Signs:
%     a, k < 0;
%     b, d, g, h, p > 0
%     c >= 0.
% The conditions assure a solution x >= -a/b exists.
%   Outputs: x is dp/dz. f is the value of the scalar equation being solved
% (ideally 0 on output). F and zeta are quantities in the science
% formulas. nit is the number of iterations in the Newton iteration.
  if (nargin < 11) x0 = []; end
  f = 0;
  nit = 1;
  
  % Type 1:
  %     zeta <= 0 (F = 1)
  %     d (b x + a) + k = 0.
  % The condition on zeta implies
  %     (i) c = 0 or (ii) x < -a/b
  % and the equation gives
  %     x = -(k/d + a)/b.
  % As a, k < 0,
  %     x = -(k/d + a)/b = -k/(d b) - a/b > -a/b,
  % which shows that only condition (i) can give a type-1 solution. From the
  % other direction, if c = 0, then all is well.
  if (c == 0)
    x = -(k/d + a)/b;
    zeta = c/(a + b*x);
    F = 1;
    return;
  end
    
  % Type 2:
  %     zeta >= 1 (F = 0)
  %     g sinh(h (a + b x)) + k = 0.
  % The condition on zeta implies c > 0 and
  %     0 < a + b x <= c => -a/b < x <= (c - a)/b.
  % The equation gives
  %     x = (asinh(-k/g)/h - a)/b.
  % First, x > -a/b, so the left inequality is trivially satisfied. The right
  % implies that a type-2 solution holds only if
  %     c >= asinh(-k/g)/h.
  akgh = asinh(-k/g)/h;
  if (c >= akgh)
    x = (akgh - a)/b;
    zeta = c/(a + b*x);
    F = 0;
    return;
  end
  
  % Type 3.
  %   Now we know 0 < F < 1, so we have a smooth function that we can
  % solve using Newton's method. For there always to exist a solution, we
  % must have
  %     (i) 0 < c < asinh(-k/g)/h.
  % The condition on zeta is
  %     0 < zeta < 1,
  % which implies
  %     0 < c/(a + b x) < 1
  %     0 < c/(a + b x) => a + b x > 0 => x > -a/b
  %     c/(a + b x) < 1 => a + b x > c => x > (c - a)/b.
  % As c > 0, we get the condition
  %     (ii) x > (c - a)/b.
  % So now it remains to determine whether the full equation always has such a
  % solution. Define
  %     f(x) = d F(x) (a + b x) + g sinh(h (a + b x)) (1 - F(x))^p + k.
  % At one extreme,
  %     lim_{x -> oo} {zeta, F, f(x)} = {0, 1, inf}.
  % At the other, x = (c - a)/b, with zeta = 1 and F = 0. Hence
  %     f(x) = g sinh(h c) + k.
  % From (i), we have an upper bound on c; substituting this, we get
  %     f(x) < 0.
  % Hence we have a zero crossing in the domain given by (ii).

  da = d*a;
  db = d*b;
  % We know this to be true, and we need it to be to proceed safely:
  x_lb = (c - a)/b;
  if (isempty(x0))
    F = 0.5;
    x = -(k + da*F + g*h*a*(1 - F)^p)/(db*F + g*h*b*(1 - F)^p);
  else
    x = x0;
  end
  x_lb_f = 1.1;
  if (x < x_lb) x = x_lb_f*x_lb; end
  x_ub = inf; % [x_lb, x_ub] is the safeguarding bracket.
  dx = inf;
  nit = 0;
  done = 0;
  while (1)
    %pr('%1.5e [%1.5e %1.5e] %1.2e %d\n', x, x_lb, x_ub, f, nit);
    % Update f.
    F = mb_solve_calc_F(a, b, c, x, 0);
    dbxda = db*x + da;
    habx = h*(a + b*x);
    gs = g*sinh(habx);
    omFp = (1 - F)^p;
    f = dbxda*F + gs*omFp + k;
    if (done) break; end
    is_bad = isnan(f) || isinf(f);
    if (is_bad && nit == 0 && ~isempty(x0))
      % Input value gave a bad f, so try a default value.
      nit = nit + 1;
      x = x_lb_f*x_lb;
      continue;
    end
    % Done?
    if (abs(f) <= tolf)
      % We don't need to recalc anything, so just break.
      break;
    end
    % Update the bracket.
    if (f > 0 || isnan(f))
      x_ub = x;
    else
      x_lb = x;
    end
    % Done?
    if (x_ub - x_lb <= tolx*x_lb)
      done = 1;
      continue;
    end
    % Get the gradient.
    F_x = mb_solve_calc_F_x(a, b, c, x);
    f_x = dbxda*F_x + db*F - gs*p*(1 - F)^(p - 1)*F_x + g*cosh(habx)*h*b*omFp;
    % Proposal step.
    dx_prev = dx;
    dx = -f/f_x;
    x = x + dx;
    % Done?
    if (abs(dx) <= tolx*x_lb)
      done = 1;
      continue;
    end
    % Safeguard it. The condition on the decrease of f is used in the Numerical
    % Recipes routine 'rtsafe'.
    if (isinf(x) || isnan(x) || ...
        x <= x_lb || x >= x_ub || ...
        abs(2*f) > abs(f_x*dx_prev))
      if (isinf(x_ub))
        x_lb_f = 2*x_lb_f;
        x = x_lb_f*x_lb;
      else
        x = 0.5*(x_lb + x_ub);
      end
    end
    nit = nit + 1;
    if (nit > 10000)
      % For safety, keep this loop bounded.
      error('nit 10000');
    end
  end
  zeta = c/(a + b*x);
end

function F = mb_solve_calc_F (a, b, c, x, clamp)
  if (nargin < 5) clamp = 1; end
  zeta = c/(a + b*x);
  zeta4 = zeta*zeta; zeta4 = zeta4*zeta4;
  if (clamp)
    if (zeta >= 1)
      F = 0;
      return;
    elseif (zeta <= 0)
      F = 1;
      return;
    end
  end
  F = 1 + (zeta4 - 4*zeta)/3;
end

function F_x = mb_solve_calc_F_x (a, b, c, x)
%assume 0 < F < 1
  ooab = 1/(a + b*x);
  zeta = c*ooab;
  zeta_x = -b*zeta*ooab;
  zeta3 = zeta*zeta*zeta;
  F_x = (4*zeta3 - 4)*zeta_x/3;
end

function [x f F zeta nit] = mbv_solve (a, b, d, g, h, k, tolx, tolf, x0)
% Solve the momentum balance equation for dp/dz (here x):
%     d (a + b x) + g sinh(h (a + b x)) + k = 0.
% Signs:
%     a, k < 0;
%     b, d, g, h, p > 0.
% No tau-yield. Frictional slip everywhere.
  if (nargin < 9) x0 = []; end
  f = 0;
  nit = 1;
  % Placeholder values:
  F = 1; zeta = 0;
    
  % We know this to be true, and we need it to be to proceed safely:
  x_lb = 0;
  if (isempty(x0))
    x = -(k + d*a + g*h*a)/(d*b + g*h*b);
  else
    x = x0;
  end
  x_lb_f = 1.1;
  if (x < x_lb) x = x_lb; end
  x_ub = inf; % [x_lb, x_ub] is the safeguarding bracket.
  dx = inf;
  nit = 0;
  done = 0;
  while (1)
    % Update f.
    abx = a + b*x;
    f = d*abx + g*sinh(h*abx) + k;
    %pr('%1.5e [%1.5e %1.5e] %1.2e %d\n', x, x_lb, x_ub, f, nit);
    if (done) break; end
    is_bad = isnan(f) || isinf(f);
    if (is_bad && nit == 0 && ~isempty(x0))
      % Input value gave a bad f, so try a default value.
      nit = nit + 1;
      x = x_lb_f*x_lb;
      continue;
    end
    % Done?
    if (abs(f) <= tolf)
      % We don't need to recalc anything, so just break.
      break;
    end
    % Update the bracket.
    if (f > 0 || isnan(f))
      x_ub = x;
    else
      x_lb = x;
    end
    % Done?
    if (x_ub - x_lb <= tolx*x_lb)
      done = 1;
      continue;
    end
    % Get the gradient.
    f_x = d*b + g*cosh(h*abx)*h*b;
    % Proposal step.
    dx_prev = dx;
    dx = -f/f_x;
    x = x + dx;
    % Done?
    if (abs(dx) <= tolx*x_lb)
      done = 1;
      continue;
    end
    % Safeguard it. The condition on the decrease of f is used in the Numerical
    % Recipes routine 'rtsafe'.
    if (isinf(x) || isnan(x) || ...
        x <= x_lb || x >= x_ub || ...
        abs(2*f) > abs(f_x*dx_prev))
      if (isinf(x_ub))
        assert(x > 0);
        x_lb_f = 2*x_lb_f;
        x = x_lb_f*x_lb;
      else
        x = 0.5*(x_lb + x_ub);
      end
    end
    nit = nit + 1;
    if (nit > 10000)
      % For safety, keep this loop bounded.
      warning('nit 10000');
      break;
    end
  end
end

function [mw phi_g] = solve_ch_bcs (m, options)
% Given p_ch, calculate mw and phi_g at the chamber according to the two BC
% equations.
  
  % Solve r == 0.
  function r = calc_r_mw (mw)
    Gamma = (1 - mw)./(m.B*mw);
    [Chi_hd Chi_cd] = solubility_liu_explicit(m.p_ch, m.T, mw);
    r = Gamma*1e-2*(m.chi_ch.total.h2o - Chi_hd) -...
        1e-6*(m.chi_ch.total.co2 - Chi_cd);
  end
  mw = fzero(@calc_r_mw, [1e-2 1]);

  % Solve rho_g(1)*phi_g(1) - total_exsolved*phi_l(1)*m.rho_l == 0.
  [Chi_hd Chi_cd] = solubility_liu_explicit(m.p_ch, m.T, mw);
  Chi_hd = 1e-2*Chi_hd;
  Chi_cd = 1e-6*Chi_cd;
  total_exsolved = 1e-2*m.chi_ch.total.h2o + ...
      1e-6*m.chi_ch.total.co2 - Chi_hd - Chi_cd;
  c1 = 1./(1 - Chi_hd - Chi_cd);    %editted c1,c2 Oct 2, 2017
  c2 = 1./(1 + Chi_hd.*c1*(m.rho_l/m.rho_hd) + Chi_cd.*c1*(m.rho_l/m.rho_cd));
  c12 = c1.*c2;
  phi_s_of_p = calc_phi_s_of_p(m, m.p_ch);
  rho_g = mw.*m.p_ch/(m.Rw*m.T) + (1 - mw).*m.p_ch/(m.Rc*m.T);
  d = (total_exsolved.*m.rho_l.*c12).*...
      (1 - phi_s_of_p.*m.rho_l.*c12./...
       (m.rho_s + phi_s_of_p.*(c12.*m.rho_l - m.rho_s)));
  phi_g = d/(rho_g + d);
end

function phi_s_calc = calc_phi_s_of_p (m, p)
  phi_s_calc = pfl_val(m.pf, p*1e-6);
end

function phi_s_calc = calc_phi_s_of_p1 (m, p)
% Old version of calc_phi_s_of_p.
  global xtl_lo;
  if (isempty(xtl_lo) || ~isfield(xtl_lo, 'Ti'))
    xtl_lo = load('xtl_tables.mat');
    xtl_lo.Ti = -1;
  end
  if (m.Ti ~= xtl_lo.Ti)
    xtl_lo.Ti = m.Ti;
    xtl_lo.tsi = flipud(xtl_lo.ts{m.Ti});
    xtl_lo.pp = interp1(xtl_lo.tsi(:,2), xtl_lo.tsi(:,1), 'pchip', 'pp');
  end
  ts = xtl_lo.tsi;
  p = p*1e-6;
  if (p >= ts(end,2))
    phi_s_calc = ts(end,1);
  elseif (p <= ts(1,2))
    phi_s_calc = ts(1,1);
  else
    if (0)
      i = find(p >= ts(:,2), 1, 'last');
      a = (p - ts(i,2))/(ts(i+1,2) - ts(i,2));
      phi_s_calc = (1 - a)*ts(i,1) + a*ts(i+1,1);
    else
      phi_s_calc = ppval(xtl_lo.pp, p);
    end
  end
end

function k_lat = calc_k_lat (m, z, phi_g, state)
  k_lat = calc_k(m, phi_g, state);
  % What was computed so far is k_lat_magma. If ~isinf(k_lat_wall), then the
  % wall has to be accounted for.
  if (m.klw.use && k_lat > 0)
%     k_lat_magma = k_lat;
    k_lat_magma = m.k_lat.*phi_g.^3;
    k_lat_wall = m.klw.top./(1e-3*z).^m.klw.mi;
    k_lat = (m.R+m.klw.L)./(m.R./k_lat_magma + m.klw.L./k_lat_wall);
  end
end

function k = calc_k (m, phi_g, state)
  if (is_thr_perm(m.k_lat_model) && ~state.k)
    k = 0;
    return;
  end
  switch (m.k_lat_model)
    case {0 3}
      k = m.k_lat;
    case {1 2}
        k = m.k_lat.*phi_g.^3;
%       k = m.k_lat.*(phi_g-m.phi_gc).^3;
    otherwise
      error('not a k_lat_model');
  end
end

function klw = klw_init (oklw)
  klw.use = ~isinf(oklw.top);
  if (~klw.use) 
      klw.L = 0;
      return; end
  klw.top = oklw.top;
  klw.mi = oklw.mi;
  klw.L = oklw.L;
end

function b = is_thr_perm (k_lat_model)
  b = any(k_lat_model == [2 3]);
end

function tauY = calc_yield (ys, phi_s, phi_l)
% Yield strength as a function of xtal volume fraction, after Zhou et al
% (1995, Rheologica Acta, equation 8).
%   Modified so that xtal fraction is taken relative to solid + liquid
% fraction. This change made so that yield strength doesn't vanish if the
% porosity gets too high.
%   Also modified so that the yield strength does not become singular at
% phim, by introducing a yield strength of the solid constituent.
  phi0 = ys.phi0;
  phim = ys.phim;
  tauc = ys.tauc;   % constant factor [Pa]
  tauys = ys.tauys; % yield strength of solid
  gamma = ys.gamma; % exponent

  % Just for this subroutine define:
  phi_s = phi_s./(phi_s + phi_l);

  tauY = tauc*((phi_s/phi0 - 1)./(1 - phi_s/phim)).^gamma;
  I1 = phi_s < phi0;
  tauY(I1) = 0;

  % Determine the solid fraction that gives solid yield strength .
  R = (tauys/tauc)^(1/gamma);
  phiys = phi0*(1 + R)/(1 + R*phi0/phim);

  I2 = phi_s > phiys;
  tauY(I2) = tauys;
end

function eta = calc_eta (c, phi, T, meltmodel, gdot)
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
        eta = 10.^((-3.545 + 0.833*lc) + (9601 - 2368*lc)./(T - (195.7 + 32.25*lc)));
    case {2}
        % Whittington et al 2009
        c = 100*c;
        eta = 10.^(-4.43 + (7618.3-17.25*log10(c+0.26))./(T-406.1+292.6*log10(c+0.26)));
    otherwise
        error('not a viscosity model');
end

if gdot==0
    % From Costa 2005b and 2007b
    kappa = 0.999916;
    phistar = 0.673;
    gamma = 3.98937;
%     delta = 1;
    delta = 16.9386;
else
    % From Costa 2005 and Caricchi 2007
    phistar = -0.066499*tanh(0.913424*log10(gdot)+3.850623) + 0.591806;
    delta   = -6.301095*tanh(0.818496*log10(gdot)+2.86) + 7.462405;
    kappa   = -0.000378*tanh(1.148101*log10(gdot)+3.92) + 0.999572;
    gamma   = 3.987815*tanh(0.8908*log10(gdot)+3.24) + 5.099645;
end

pih = 1.772453850905515882;
B = 2.5;
eta_phi = (1 + (phi./phistar).^delta) ./ ...
    (1 - kappa*erf((pih.*phi)./(2*kappa*phistar).* ...
    (1 + (phi./phistar).^gamma)) ...
    ).^(B*phistar);

% Brouwers model
% phi_m = 0.8;
% eta_phi = ((1-phi)./(1-phi/phi_m)).^(5/2*phi_m/(1-phi_m));
% eta_phi(phi>phi_m) = inf;


eta = eta.*eta_phi;
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
%     if (sum(Ih + Ic) > 0) warning('Melt is undersaturated!'); end
  end
end
