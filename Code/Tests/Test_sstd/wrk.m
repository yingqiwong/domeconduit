function varargout = wrk (varargin)
  [varargout{1:nargout}] = feval(varargin{:});
end

function des = run (o, Niter)
  if nargin == 1, Niter = 20; end

  o.slv.de.reltol = 1e-10;
  o.slv.de.abstol = 1e-20;
  o.slv.zero.tolx = 1e-10;
  o.slv.max_time  = 1000;
  des = plot_fs(o, Niter);
end

function des = plot_fs (o, nva)
% o is an options struct for smf. nva is the number of v_ch points at which to
% evaluate the ODE. Given these, this routine finds a bracket for the zero
% crossing of
%     f_ch(v_ch) = p_top(v_ch) - m.p_top.
% Then it evaluates f_ch(v_ch) at a bunch of v_ch within this bracket and plots
% f_ch(v_ch).  
  m = smf_rad_dz('m_init', o);
  init_gsm();
  global gsm;
  va = m.v_ch_guess; 
  [vi fi flag msg] = smf_rad_dz('solve1_get_initial_interval', @pr, m, va);
  if (flag)
    sprintf('When bracketing: exiting with flag %d (%s).\n', flag, msg);
    return;
  end
  if (numel(vi) == 2)
    sprintf('Interval: [%1.2e %1.2e] with f [%1.2e %.1e2] %d evals\n', ...
       vi(1), vi(2), fi(1), fi(2), gsm.cnt);
  else
    sprintf('v_ch %1.2e %d evals\n', vi, gsm.cnt);
  end
 
  vas = logspace(log10(1e-10*vi(1)), log10(vi(2)), nva);
  fs = zeros(1, nva);
  fs([1 end]) = fi;
  flag = zeros(1, nva);
  des = {};
  
  for (i = 1:nva)
    [fs(i) flag(i) des{i}] = smf_rad_dz('sm_calc_fsm', m, vas(i));
    des{i}.o = o;
    des{i}.m = m;
    des{i}.v_ch = vas(i);
    des{i}.fs = fs(i);
    
    % calculate rho and therefore q = rho*v at conduit inlet
    des_ch = des{i}; des_ch.z = des{i}.z(1);
    var_ch = smf_rad_dz('add_fields',des_ch);
    des{i}.q_ch = var_ch.rho*vas(i);
    
    clf; semilogx(vas, fs, 'o-'); drawnow;
  end
  
  hold on; plot(xlim, [0,0], 'k-'); hold off;
  xlabel('v_{ch}'); ylabel('p_{top}(v_{ch}) - m.p_{top}');
 
end

function init_gsm ()
% Keep a global struct.
  clear global gsm;
  global gsm;
  gsm.f = inf; %
  gsm.de = []; % Reuse the two final solutions from the bracketing step.
  gsm.vi = []; %
  gsm.cnt = 0;       % Number of calls to the fsm(v_ch) function.
  gsm.tmr.t = tic(); % Timer to enforce a timeout (m.slv.max_time).
end

function check_soln (s)
% Check solutions using finite differences. For d in des from plot_fs:
%   d.o.slv.use_mex=0; s=smf('add_fields',d);
%   wrk('check_soln',s);
  water_mass = -central_diff(s.water_mass.*s.v,s.z);
  water_vert = -central_diff(s.water_vert,s.z);
  water_loss = s.water_loss;
  C_mass = -central_diff(s.C_mass.*s.v,s.z);
  C_vert = -central_diff(s.C_vert,s.z);
  C_loss = s.C_loss;
  figure(1);clf;
  subplot(411); plot(s.z, water_mass + water_vert, '.-', s.z, -water_loss, 'o-');
  xlim(s.z([end 1]));
  ylim([-1e-3 0]);
  subplot(412); plot(s.z, C_mass + C_vert, '.-', s.z, -C_loss, 'o-');
  ylim([-6e-5 0])
  try
    % I separated out this variable in des_calc_exprs just to check that it's ~0.
    non_volatile_continuity = central_diff(s.non_volatile_continuity, s.z);
    subplot(413); plot(s.z, non_volatile_continuity, '.-');
    ylim([-1 1]*1e-12);
  catch
  end
  subplot(414); plot(s.z, s.state, '.-');
 
end

function [o] = model_struct(model, o, is)

o.p_ch  = model(is.p_ch);
o.k_lat = 10^(model(is.k_lat));
o.phi_gc = model(is.phi_gc);
o.total_h2o = model(is.h2o);
o.fr.v_ref = 1e-5;  % units m/s
o.fr.a  = 10^model(is.a);  % dimensionless
o.fr.f0 = 10^model(is.f0);  % dimensionless
o.R = model(is.R);  % dimensionless

end
