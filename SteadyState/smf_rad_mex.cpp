/* Speed up the bottlenecks in smf.m. To use, set
       o.slv.use_mex = 1
   As long as you have mex configured with a C++ compiler, to build type in
   Matlab
       mex smfmex.cpp

   Formulas in this file are transcribed from smf.m and should not be considered
   as documentary; see smf.m for comments, units, etc.

   Oct 2013. AMB ambrad@cs.stanford.edu. Initial.
   Jul 2014. AMB. See smf.m.
 */

#include "string.h"
#include "stdio.h"
#include "math.h"
#include <vector>
#include <cmath>
#include <string>
#include <algorithm>
#include <limits>
using namespace std;

typedef struct mxArray_tag mxArray;

// -----------------------------------------------------------------------------
// Tools.

// Mimic ppval(pp, x) for cubic segments.
class CubicSpline {
public:
  CubicSpline () {}
  //spec size(breaks) == n + 1
  //spec size(coefs) == 4 n
  //spec coefs packed by segment, from x^0 to x^3.
  CubicSpline (const double* breaks, const double* coefs, const size_t n)
    : n(n)
  {
    this->breaks.resize(n + 1);
    memcpy(&this->breaks[0], breaks, (n + 1)*sizeof(*breaks));
    this->coefs.resize(4*n);
    memcpy(&this->coefs[0], coefs, 4*n*sizeof(*coefs));
  }

  double call (double x) const {
    const size_t bi = get_segment(x);
    const double* c = &coefs[0] + 4*bi;
    x -= breaks[bi];
    return c[3] + x*(c[2] + x*(c[1] + x*c[0]));
  }

  size_t get_segment (double x) const {
    if (x <= breaks[0]) return 0;
    else if (x >= breaks[n-1]) return n - 1;
    else {
      const double* b = lower_bound(&breaks[0], &breaks[n-1], x);
      return (size_t) (b - &breaks[0]) - 1;
    }
  }

  double get_x_low () const { return breaks[0]; }
  double get_y_low () const { return coefs[3]; }
  double get_x_high () const { return breaks[n]; }
  double get_y_high () const { return coefs[4*n - 1]; }

private:
  size_t n;
  vector<double> breaks, coefs;

  friend void CubicSpline_init (CubicSpline& cs, const mxArray* mpp);
};

void panic (const string& err);

// -----------------------------------------------------------------------------
// Model.

struct Model {
  enum ModelType { bingham, newtonian } model;
  struct Yield { double phi0, phim, tauc, tauys, gamma, phiys; } ys;
  struct Friction { double v_ref, f0, a, Psi; } fr;
  struct { double slope, offset; } sig;
  struct Permeability {
    double k, eta_g, phydro_slope, L_diff, phi_gc;
    struct {
      bool use;
      double top, mi, L;
    } klw;
    int model, plug_gas_loss;
    bool use_k_vert;
  } k;
  struct Viscosity { int melt, use_strain; } eta;
  struct {
    enum Yis { p = 0, v, phi_g, mw, w };
    enum Fis { mb = 0, nv, h2o, co2, fw };
    double lb[4], ub[4];
  } blk;
  struct {
    double fd;
    short fdo;
  } slv;
  CubicSpline pf;
  double p_top, R, g, conduit_length, expn, T, Rw, Rc, B, rho_l, rho_hd, rho_cd,
    rho_s, rho_c, sy[4], syp[4], sF[4];
  bool use_phi_s_ratio;
  bool inited;

  Model () : inited(false) {}

  void finalize () {
    fr.Psi = 2*fr.v_ref/exp(fr.f0/fr.a);

    const double R = pow(ys.tauys/ys.tauc, 1/ys.gamma);
    ys.phiys = ys.phi0*(1 + R)/(1 + R*ys.phi0/ys.phim);

    inited = true;
  }
};

struct State { int k, plug; };

inline double calc_phi_s_of_p (const CubicSpline& cs, double p) {
  p *= 1e-6;
  if (p <= cs.get_x_low()) return cs.get_y_low();
  if (p >= cs.get_x_high()) return cs.get_y_high();
  return cs.call(p);  
}

inline void calc_solubility_liu_explicit (
  const double B, double p, const double T, const double mw,
  double& Chi_hd, double& Chi_cd)
{
  p *= 1e-6;
  
  const double
    Pw = p*mw,
    Pc = p*(1 - mw),
    Pwh = sqrt(Pw),
    Pw15 = Pw*Pwh;
  Chi_hd = (354.94*Pwh + 9.623*Pw - 1.5223*Pw15)/T +
    0.0012439*Pw15 + Pc*(-1.084e-4*Pwh - 1.362e-5*Pw);
  Chi_cd = Pc*(5668 - 55.99*Pw)/T + Pc*(0.4133*Pwh + 2.041e-3*Pw15);

  Chi_hd *= 1e-2;
  Chi_cd *= 1e-6;
}

inline double
calc_yield (const Model::Yield& y, double phi_s, const double phi_l) {
  phi_s /= phi_s + phi_l;
  if (phi_s < y.phi0) return 0;
  if (phi_s > y.phiys) return y.tauys;
  return y.tauc*pow((phi_s/y.phi0 - 1)/(1 - phi_s/y.phim), y.gamma);
}

inline double calc_eta (const double c, const double phi, const double T, 
                        const int meltmodel, const double gdot) {
    
  // initialize
  double eta, phistar, delta, kappa, gamma;
  
  if (meltmodel==1) {
      const double lc = log(1e2*c);
      eta = pow(10, (-3.545 + 0.833*lc) + (9601 - 2368*lc)/(T - (195.7 + 32.25*lc)));
  }
  else if (meltmodel==2) {
      eta = pow(10, -4.43 + (7618.3-17.25*log10(100*c+0.26))/(T-406.1+292.6*log10(100*c+0.26)));
  }

  static const double B = 2.5, pih = 1.772453850905515882;
  
  if (gdot==0) {
      kappa = 0.999916, phistar = 0.673, gamma = 3.98937, delta = 16.9386;
  }
  else {
      phistar = -0.066499*tanh(0.913424*log10(gdot)+3.850623) + 0.591806,
      delta = -6.301095*tanh(0.818496*log10(gdot)+2.86) + 7.462405,
      kappa = -0.000378*tanh(1.148101*log10(gdot)+3.92) + 0.999572,
      gamma = 3.987815*tanh(0.8908*log10(gdot)+3.24) + 5.099645;
  }
          
  eta *= ((1 + pow(phi/phistar, delta)) /
          pow(1 - kappa*erf(pih*phi/(2*kappa*phistar)*(1 + pow(phi/phistar, gamma))),
              B*phistar));
  // eta = 0.1*eta;
  return eta;
}

inline double pow3 (const double x) { return x*x*x; }
inline double pow4 (double x) { x = x*x; return x*x; }

inline bool is_thr_perm (int model) { return model == 2 || model == 3; }

inline void 
calc_permeability (const Model& m, const double z,
                   const double phi_g, double& k_lat, double& k_vert) {
  switch (m.k.model) {
  case 0:
  case 3:
    k_vert = m.k.k;
    break;
  case 1:
  case 2:
    k_vert = m.k.k*pow3(phi_g);
    break;
  }
  k_lat = k_vert;
  if (m.k.klw.use && k_lat > 0) {
    const double k_lat_wall = m.k.klw.top / pow(1e-3*z, m.k.klw.mi);
    k_lat = (m.R+m.k.klw.L)/(m.R/k_lat + m.k.klw.L/k_lat_wall);
  }
  if (!m.k.use_k_vert) k_vert = 0;
}

namespace mb {
inline double mypow (const double x, const double p) {
  if (p == 0) return 1;
  if (p == 1) return x;
  return pow(x, p);
}

struct Calc {
  const double a, b, c, d, g, h, k, p;
  Calc (const double a, const double b, const double c, const double d,
        const double g, const double h, const double k, const double p)
    : a(a), b(b), c(c), d(d), g(g), h(h), k(k), p(p) {}

  double abx, ooabx, zeta, zeta3, gs, omFp;

  void calc_f (const double x, double& f, double& F) {
    abx = a + b*x;
    ooabx = 1/abx;
    zeta = c*ooabx;
    zeta3 = zeta*zeta*zeta;
    F = 1 + (zeta3 - 4)*zeta/3;
    gs = g*sinh(h*abx);
    omFp = mypow(1 - F, p);
    f = d*abx*F + gs*omFp + k;
  }

  double calc_f_x (const double x, const double F) const {
    const double
      zeta_x = -b*zeta*ooabx,
      F_x = (4*zeta3 - 4)*zeta_x/3;
    return
      d*abx*F_x + d*b*F - gs*p*mypow(1 - F, p - 1)*F_x + g*cosh(h*abx)*h*b*omFp;
  }
};

void solve (const double a, const double b, const double c, const double d,
            const double g, const double h, const double k, const double p,
            const double tolx, const double tolf, const double x0,
            double& x, double& F, double& zeta, double& f) {
  f = 0;

  // Type 1 solution.
  if (c == 0) {
    x = -(k/d + a)/b;
    zeta = c/(a + b*x);
    F = 1;
    return;
  }

  // Type 2.
  double akgh = asinh(-k/g)/h;
  if (c >= akgh) {
    x = (akgh - a)/b;
    zeta = c/(a + b*x);
    F = 0;
    return;
  }

  // Type 3.
  double
    x_lb = (c - a)/b,
    x_lb_f = 1.1,
    x_ub = numeric_limits<double>::infinity(),
    dx = numeric_limits<double>::infinity();
  bool done = false;
  const bool x0_input = x0 > x_lb;
  if (x0_input)
    x = x0;
  else {
    F = 0.5;
    const double omFp = mypow(1 - F, p);
    x = -(k + d*a*F + g*h*a*omFp) / (d*b*F + g*h*b*omFp);
  }
  if (x < x_lb) x = x_lb_f*x_lb;
  for (size_t nit = 0; ; nit++) {
    //printf("c %1.5e [%1.5e %1.5e] %1.2e %ld\n", x, x_lb, x_ub, f, nit);
    Calc calc(a, b, c, d, g, h, k, p);
    // Update f.
    calc.calc_f(x, f, F);
    if (done) break;
    const bool is_bad = std::isinf(f) || std::isnan(f);
    if (is_bad && nit == 0 && x0_input) {
      // Input x was not good, so try our own.
      x = x_lb_f*x_lb;
      continue;
    }
    // Done?
    if (fabs(f) <= tolf) break;
    // Update the bracket.
    if (f > 0 || std::isnan(f)) x_ub = x; else x_lb = x;
    // Done?
    if (x_ub - x_lb <= tolx*x_lb) {
      done = true;
      continue;
    }
    // Get the gradient.
    const double f_x = calc.calc_f_x(x, F);
    // Proposal step.
    const double dx_prev = dx;
    dx = -f / f_x;
    x += dx;
    // Done?
    if (fabs(dx) <= tolx*x_lb) {
      done = true;
      continue;
    }
    // Safeguard the step. The condition on the decrease of f is used in the
    // Numerical Recipes routine 'rtsafe'.
    if (std::isnan(x) || std::isinf(x) ||
        x <= x_lb || x >= x_ub ||
        fabs(2*f) > fabs(f_x*dx_prev)) {
      if (std::isinf(x_ub)) {
        x_lb_f = 2*x_lb_f;
        x = x_lb_f*x_lb;
      } else
        x = 0.5*(x_lb + x_ub);
    }
    if (nit == 10000) {
      // For safety, keep this loop bounded.
      panic("nit 10000");
    }
  }
  zeta = c/(a + b*x);
}
} // namespace mb

namespace mbv {
void solve (const double a, const double b, const double d, const double g,
            const double h, const double k, const double tolx,
            const double tolf, const double x0, double& x, double& F,
            double& zeta, double& f) {
  f = 0;
  F = 1;
  zeta = 0;

  double
    x_lb = 0,
    x_lb_f = 1.1,
    x_ub = numeric_limits<double>::infinity(),
    dx = numeric_limits<double>::infinity();
  bool done = false;
  const bool x0_input = x0 > x_lb;
  if (x0_input)
    x = x0;
  else {
    x = -(k + d*a + g*h*a)/(d*b + g*h*b);
  }
  if (x < x_lb) x = x_lb_f*x_lb;
  for (size_t nit = 0; ; nit++) {
    //printf("c %1.5e [%1.5e %1.5e] %1.2e %ld\n", x, x_lb, x_ub, f, nit);
    // Update f.
    const double abx = a + b*x;
    f = d*abx + g*sinh(h*abx) + k;
    if (done) break;
    const bool is_bad = std::isinf(f) || std::isnan(f);
    if (is_bad && nit == 0 && x0_input) {
      // Input x was not good, so try our own.
      x = x_lb_f*x_lb;
      continue;
    }
    // Done?
    if (fabs(f) <= tolf) break;
    // Update the bracket.
    if (f > 0 || std::isnan(f)) x_ub = x; else x_lb = x;
    // Done?
    if (x_ub - x_lb <= tolx*x_lb) {
      done = true;
      continue;
    }
    // Get the gradient.
    const double f_x = d*b + g*cosh(h*abx)*h*b;
    // Proposal step.
    const double dx_prev = dx;
    dx = -f / f_x;
    x += dx;
    // Done?
    if (fabs(dx) <= tolx*x_lb) {
      done = true;
      continue;
    }
    // Safeguard the step. The condition on the decrease of f is used in the
    // Numerical Recipes routine 'rtsafe'.
    if (std::isnan(x) || std::isinf(x) ||
        x <= x_lb || x >= x_ub ||
        fabs(2*f) > fabs(f_x*dx_prev)) {
      if (std::isinf(x_ub)) {
        if (x == 0) panic("x == 0 in mbv::solve safeguard; shouldn't happen.");
        x_lb_f = 2*x_lb_f;
        x = x_lb_f*x_lb;
      } else
        x = 0.5*(x_lb + x_ub);
    }
    if (nit == 10000) {
      // For safety, keep this loop bounded.
      panic("nit 10000");
    }
  }
}
} // namespace mbv

struct Exprs {
  enum Equations { emb = 1, envol, evol };
  double F, vvfrac, zeta, pp, dbg;
  struct Fg {
    double f, g;
    Fg () { f = g = 0; }
  } nvol, h2o, co2;
  struct Mbe {
    double vv, vf;
    Mbe () { vv = vf = 0; }
  } mbe;
  Exprs () { F = zeta = pp = dbg = 0; }
};

struct Measurements {
  double p, phi_g, v, phi_l, F;
  void set (const double ip, const double iphi_g, const double iv,
            const double iphi_l, const double iF) {
    p = ip; phi_g = iphi_g; v = iv; phi_l = iphi_l; F = iF;
  }
  void get (double& op, double& ophi_g, double& ov, double& ophi_l, double& oF)
    const { op = p; ophi_g = phi_g; ov = v; ophi_l = phi_l; oF = F; }
};

void calc_exprs (
  const Model& m, const State& state, const double z, const double* const y,
  const double* const yp, const vector<Exprs::Equations>& eqs, Exprs& e,
  Measurements* ms = NULL)
{
  if (!m.inited) panic("m was not inited.");

  const double p = y[m.blk.p], v = y[m.blk.v], phi_g = y[m.blk.phi_g],
    mw = y[m.blk.mw];

  double Chi_hd, Chi_cd;
  calc_solubility_liu_explicit(m.B, p, m.T, mw, Chi_hd, Chi_cd);

  const double
    Gamma = (1 - mw)/(mw*m.B),
    c1 = 1/(1 - Chi_hd - Chi_cd),
    c2 = 1/(1 + Chi_hd*c1*(m.rho_l/m.rho_hd) + Chi_cd*c1*(m.rho_l/m.rho_cd)),
    c12 = c1*c2,
    phi_s_of_p = calc_phi_s_of_p(m.pf, p),
    phi_s = (phi_s_of_p*m.rho_l*c12*(1 - phi_g)/
             (m.rho_s + phi_s_of_p*(c12*m.rho_l - m.rho_s))),
    phi_l = (1 - phi_s - phi_g)*c2,
    rho_g = p*(mw/(m.Rw*m.T) + (1 - mw)/(m.Rc*m.T)),
    rho = m.rho_l*phi_l*c1 + m.rho_s*phi_s + rho_g*phi_g;
  
  double Visc, Den;
  { // Solve the momentum balance equation for dp/dz.
    const double
      phi_s_eta = m.use_phi_s_ratio ? phi_s/(1 - phi_g) : phi_s,
      gdot = m.eta.use_strain ? 2*v/m.R : 0,
      eta = calc_eta(Chi_hd, phi_s_eta, m.T, m.eta.melt, gdot),
      tolx = 1e2*std::numeric_limits<double>::epsilon(),
      tolf = 1e2*v*std::numeric_limits<double>::epsilon();
    Visc = 0.25*m.R/eta;
    Den = m.fr.a*(m.sig.slope*z + m.sig.offset);
    double pp0 = yp ? yp[m.blk.p] : -1;
    double f;
    switch (m.model) {
    case Model::bingham: {
      const double tauY = calc_yield(m.ys, phi_s, phi_l);
      mb::solve(-0.5*m.R*rho*m.g, 0.5*m.R, tauY, Visc, m.fr.Psi, 1/Den, -v,
                m.expn, tolx, tolf, pp0, e.pp, e.F, e.zeta, f);
    } break;
    case Model::newtonian: {
      mbv::solve(-0.5*m.R*rho*m.g, 0.5*m.R, Visc, m.fr.Psi, 1/Den, -v,
                 tolx, tolf, pp0, e.pp, e.F, e.zeta, f);
    } break;
    }
  }

  for (size_t i = 0; i < eqs.size(); i++)
    switch (eqs[i]) {

    case Exprs::emb: {
      const double tauR = -0.5*m.R*(-e.pp + rho*m.g);
      switch (m.model) {
      case Model::bingham: {
        panic("Shouldn't be in MBE eq with model = bingham.\n");
      } break;
      case Model::newtonian: {
        e.mbe.vv = tauR*Visc;
        e.mbe.vf = m.fr.Psi*sinh(tauR/Den);
      } break;
      }
    } break;

    case Exprs::envol: {
      e.nvol.f = (m.rho_l*phi_l + m.rho_s*phi_s)*v;
      e.nvol.g = 0;
    } break;

    case Exprs::evol: {
      const double C_fac = Gamma/(1 + Gamma);

      if (state.plug == 1) {
          
        //  double drhogdz = yp ? (rho_g/p*(e.pp) + p*(1/(m.Rw*m.T) - 1/(m.Rc*m.T))*(-yp[m.blk.mw])) : (rho_g/p*(e.pp));
          
        //  e.h2o.f = (Chi_hd*m.rho_l*phi_l*c1 + phi_g*rho_g/(1 + Gamma))*v;
        //  e.h2o.g = 1/(1+Gamma)*phi_g*v*drhogdz;
          
        //  e.co2.f = (Chi_cd*m.rho_l*phi_l*c1 + phi_g*rho_g*C_fac)*v;
        //  e.co2.g = C_fac*phi_g*v*drhogdz;
          
        e.h2o.f = phi_g/(1 + Gamma)*v;
        e.h2o.g = 0;

        e.co2.f = phi_g*C_fac*v;
        e.co2.g = 0;
      } else {
        double u_lat = 0, vgdiff = 0;
        if (!is_thr_perm(m.k.model) || state.k == 1) {
          double k_lat, k_vert;
          calc_permeability(m, z, phi_g, k_lat, k_vert);
          u_lat = k_lat/m.k.eta_g*(p - m.k.phydro_slope*z)/(m.R+m.k.klw.L);
          if (m.k.use_k_vert) vgdiff = (k_vert/m.k.eta_g)*e.pp;
        //  if (state.plug == 1) vgdiff = m.Kplug/rho_g - v;
        }

        const double
          water_mass = Chi_hd*m.rho_l*phi_l*c1 + phi_g*rho_g/(1 + Gamma),
          C_mass = Chi_cd*m.rho_l*phi_l*c1 + phi_g*rho_g*C_fac;  
        double water_vert = 0, C_vert = 0;
        if (vgdiff != 0) {
          water_vert = phi_g*rho_g/(1 + Gamma)*vgdiff;
          C_vert = phi_g*rho_g*vgdiff*C_fac;
        }
        double water_loss = 0, C_loss = 0;
        if (u_lat != 0) {
          water_loss = 2/m.R*rho_g*phi_g*u_lat/(1 + Gamma);
          C_loss = 2/m.R*rho_g*phi_g*u_lat*C_fac;
        }
        
        e.h2o.f = water_mass*v + water_vert;
        e.h2o.g = water_loss;
        e.co2.f = C_mass*v + C_vert;
        e.co2.g = C_loss;
        
      //  e.rho_g = rho_g;
      //  e.vgdiff = vgdiff;
        
      //  if (state.plug == 1) {
      //      e.h2o.g = 0;
      //      e.co2.g = 0;
      //  }
      }

    } break;
    default:
      panic("Not an equation number.");
    }

  if (ms) ms->set(p, phi_g, v, phi_l, e.F);
}

void calc_f_y (
  const Model& m, const State& state, const double z, const double* const y,
  const double* const f, double* f_y)
{
  if (!m.inited) panic("m was not inited.");
  
  vector<Exprs::Equations> eqs(2);
  eqs[0] = Exprs::envol;
  eqs[1] = Exprs::evol;

  double yp[] = {-1, 0, 0, 0};
  for (size_t i = 0; i < 4; i++) {
    double delta = m.slv.fd*m.sy[i];

    double yd[4];
    memcpy(yd, y, 4*sizeof(double));
    yd[i] += delta;

    const bool on_bound = yd[i] > m.blk.ub[i];
    if (on_bound) {
      delta *= -1;
      yd[i] = y[i] + delta;
    }

    Exprs e;
    calc_exprs(m, state, z, yd, yp, eqs, e);
    if (i == 0) yp[i] = e.pp;

    bool done = false;
    if (m.slv.fdo == 2 && !on_bound) {
      yd[i] = y[i] - delta;
      if (yd[i] >= m.blk.lb[i]) {
        Exprs em;
        calc_exprs(m, state, z, yd, yp, eqs, em);
        delta *= 2;
        f_y[0] = (e.nvol.f - em.nvol.f)/delta;
        f_y[1] = (e.h2o.f  - em.h2o.f )/delta;
        f_y[2] = (e.co2.f  - em.co2.f )/delta;
        done = true;
      }
    }
    if (!done) {
      f_y[0] = (e.nvol.f - f[0])/delta;
      f_y[1] = (e.h2o.f  - f[1])/delta;
      f_y[2] = (e.co2.f  - f[2])/delta;
    }

    f_y += 3;
  }
}

void calc_f_z (
  const Model& m, const State& state, const double z, const double* const y,
  const double* const f, double* f_z)
{
  if (!m.inited) panic("m was not inited.");
  
  vector<Exprs::Equations> eqs(2);
  eqs[0] = Exprs::envol;
  eqs[1] = Exprs::evol;

  double yp[] = {-1, 0, 0, 0};
  double delta = m.slv.fd*m.conduit_length;

  double zd;
  zd = z + delta;

  const bool on_bound = z > m.conduit_length;
  if (on_bound) {
    delta *= -1;
    zd = z + delta;
  }

  Exprs e;
  calc_exprs(m, state, zd, y, yp, eqs, e);
  yp[0] = e.pp;

  bool done = false;
  if (m.slv.fdo == 2 && !on_bound) {
    zd = z - delta;
    if (z >= 0) {
      Exprs em;
      calc_exprs(m, state, zd, y, yp, eqs, em);
      delta *= 2;
      f_z[0] = (e.nvol.f - em.nvol.f)/delta;
      f_z[1] = (e.h2o.f  - em.h2o.f )/delta;
      f_z[2] = (e.co2.f  - em.co2.f )/delta;
      done = true;
    }
  }
  if (!done) {
    f_z[0] = (e.nvol.f - f[0])/delta;
    f_z[1] = (e.h2o.f  - f[1])/delta;
    f_z[2] = (e.co2.f  - f[2])/delta;
  }
}

inline void transform (const Model& m, double& z, double* y, double* yp) {
  z = m.conduit_length*(1 - z);
  if (y) for (size_t i = 0; i < 4; i++) y[i] *= m.sy[i];
  if (yp) for (size_t i = 0; i < 4; i++) yp[i] *= m.syp[i];
}

inline void put_in_bounds (const Model& m, double y[4]) {
  for (size_t i = 0; i < 4; i++)
    if (y[i] < m.blk.lb[i]) y[i] = m.blk.lb[i];
    else if (y[i] > m.blk.ub[i]) y[i] = m.blk.ub[i];
}

void calc_F (
  const Model& m, const State& state, double z, const double* const yi,
  const double* const ypi, double* F)
{
  if (!m.inited) panic("m was not inited.");
  double y[4], yp[4];
  memcpy(y, yi, 4*sizeof(double));
  memcpy(yp, ypi, 4*sizeof(double));
  transform(m, z, y, yp);
  put_in_bounds(m, y);

  vector<Exprs::Equations> eqs(2);
  eqs[0] = Exprs::envol;
  eqs[1] = Exprs::evol;

  Exprs e;
  calc_exprs(m, state, z, y, yp, eqs, e);

  double f[3];
  f[0] = e.nvol.f;
  f[1] = e.h2o.f;
  f[2] = e.co2.f;
  double f_y[12];
  calc_f_y(m, state, z, y, f, f_y);
  double f_z[3];
  calc_f_z(m, state, z, y, f, f_z);
  
  F[0] = yp[m.blk.p] - e.pp;
  F[1] = f_y[0]*yp[0] + f_y[3]*yp[1] + f_y[6]*yp[2] + f_y[ 9]*yp[3] + f_z[0] - e.nvol.g;
  F[2] = f_y[1]*yp[0] + f_y[4]*yp[1] + f_y[7]*yp[2] + f_y[10]*yp[3] + f_z[1] - e.h2o.g;
  F[3] = f_y[2]*yp[0] + f_y[5]*yp[1] + f_y[8]*yp[2] + f_y[11]*yp[3] + f_z[2] - e.co2.g;
  for (size_t i = 0; i < 4; i++) F[i] /= m.sF[i];
}

// -----------------------------------------------------------------------------
// Interface.

#ifndef MAIN
#include <mex.h>

inline void CubicSpline_init (CubicSpline& cs, const mxArray* mpp) {
  const mxArray* mbreaks = mxGetField(mpp, 0, "breaks");
  cs.n = mxGetNumberOfElements(mbreaks) - 1;
  cs.breaks.resize(cs.n + 1);
  cs.coefs.resize(4 * cs.n);

  memcpy(&cs.breaks[0], mxGetPr(mbreaks), (cs.n + 1)*sizeof(double));

  double* c = &cs.coefs[0];
  const double* mc = mxGetPr(mxGetField(mpp, 0, "coefs"));
  for (size_t i = 0; i < cs.n; i++) {
    for (size_t j = 0; j < 4; j++) c[j] = mc[cs.n*j + i];
    c += 4;
  }
}

inline void panic (const string& err) { mexErrMsgTxt(err.c_str()); }

Model::ModelType parse_model (const mxArray* mm) {
  char ms[2];
  mxGetString(mxGetField(mm, 0, "model"), ms, 2);
  if (ms[0] == 'n') return Model::newtonian;
  return Model::bingham;
}

inline double get_double (const mxArray* mm, const char* fld) {
  return mxGetScalar(mxGetField(mm, 0, fld));
}

inline void Model_init (Model& m, const mxArray* mm) {
#define set(s, ms, fld) s.fld = get_double(ms, #fld);
#define setm(fld) set(m, mm, fld)
  setm(p_top); setm(R); setm(g); setm(conduit_length); setm(expn); setm(T);
  setm(Rw); setm(Rc); setm(B); setm(rho_l); setm(rho_hd); setm(rho_cd);
  setm(rho_s); setm(rho_c);
#undef setm
  m.use_phi_s_ratio = static_cast<bool>(get_double(mm, "use_phi_s_ratio"));
  m.model = parse_model(mm);

  const mxArray* s;

  s = mxGetField(mm, 0, "ys");
  set(m.ys, s, phi0); set(m.ys, s, phim); set(m.ys, s, tauc);
  set(m.ys, s, tauys); set(m.ys, s, gamma);

  s = mxGetField(mm, 0, "fr");
  set(m.fr, s, v_ref); set(m.fr, s, f0); set(m.fr, s, a);

  s = mxGetField(mm, 0, "sig");
  set(m.sig, s, slope); set(m.sig, s, offset);

  CubicSpline_init(m.pf, mxGetField(mxGetField(mm, 0, "pf"), 0, "pp"));
  
  m.k.k = get_double(mm, "k_lat");
  set(m.k, mm, eta_g); set(m.k, mm, L_diff); set(m.k, mm, phi_gc);
  m.k.phydro_slope = get_double(mxGetField(mm, 0, "phydro"), "slope");
  m.k.model = static_cast<int>(get_double(mm, "k_lat_model"));
  m.k.plug_gas_loss = static_cast<int>(get_double(mm, "plug_gas_loss"));
  m.k.use_k_vert = static_cast<bool>(get_double(mm, "use_k_vert"));
  s = mxGetField(mm, 0, "klw");
  m.k.klw.use = static_cast<bool>(get_double(s, "use"));
  if (m.k.klw.use) { set(m.k.klw, s, top); set(m.k.klw, s, mi); set(m.k.klw, s, L); }
  else { m.k.klw.top = m.k.klw.mi = 0; }

  s = mxGetField(mm,0,"eta");
  set(m.eta, s, melt); set(m.eta, s, use_strain);

  s = mxGetField(mm, 0, "blk");
  memcpy(m.blk.lb, mxGetPr(mxGetField(s, 0, "lb")), 4*sizeof(double));
  memcpy(m.blk.ub, mxGetPr(mxGetField(s, 0, "ub")), 4*sizeof(double));

  s = mxGetField(mm, 0, "slv");
  set(m.slv, s, fd);
  m.slv.fdo = static_cast<short>(get_double(s, "fdo"));

  memcpy(m.sy, mxGetPr(mxGetField(mm, 0, "sy")), 4*sizeof(double));
  memcpy(m.syp, mxGetPr(mxGetField(mm, 0, "syp")), 4*sizeof(double));
  memcpy(m.sF, mxGetPr(mxGetField(mm, 0, "sF")), 4*sizeof(double));

  m.finalize();
#undef set
}

inline void mex_to_string (const mxArray* ms, string& fn) {
  int strlen = mxGetNumberOfElements(ms) + 1;
  vector<char> vfn(strlen);
  mxGetString(ms, &vfn[0], strlen);
  fn.append(&vfn[0], strlen - 1);
}

inline void State_init (State& state, const mxArray* ms) {
  state.k = static_cast<int>(get_double(ms, "k"));
  state.plug = static_cast<int>(get_double(ms, "plug"));
}

inline void Equations_init (vector<Exprs::Equations>& eqs, const mxArray* me) {
  size_t n = mxGetNumberOfElements(me);
  eqs.resize(n);
  const double* de = mxGetPr(me);
  for (size_t i = 0; i < n; i++)
    eqs[i] = static_cast<Exprs::Equations>(de[i]);
}

inline mxArray* Fg_make_mex (const Exprs::Fg& fg) {
  const char* flds[] = {"f", "g"};
  mxArray* m = mxCreateStructMatrix(1, 1, 2, flds);
  mxSetField(m, 0, flds[0], mxCreateDoubleScalar(fg.f));
  mxSetField(m, 0, flds[1], mxCreateDoubleScalar(fg.g));
  return m;
}

inline mxArray* Mbe_make_mex (const Exprs::Mbe& mbe) {
  const char* flds[] = {"vv", "vf"};
  mxArray* m = mxCreateStructMatrix(1, 1, 2, flds);
  mxSetField(m, 0, flds[0], mxCreateDoubleScalar(mbe.vv));
  mxSetField(m, 0, flds[1], mxCreateDoubleScalar(mbe.vf));
  return m;
}

inline mxArray* Exprs_make_mex (const Exprs& e) {
  const char* flds[] = {"F", "zeta", "pp", "nv", "h2o", "co2", "mbe", "dbg"};
  mxArray* me = mxCreateStructMatrix(1, 1, 8, flds);
  mxSetField(me, 0, flds[0], mxCreateDoubleScalar(e.F));
  mxSetField(me, 0, flds[1], mxCreateDoubleScalar(e.zeta));
  mxSetField(me, 0, flds[2], mxCreateDoubleScalar(e.pp));
  mxSetField(me, 0, flds[3], Fg_make_mex(e.nvol));
  mxSetField(me, 0, flds[4], Fg_make_mex(e.h2o));
  mxSetField(me, 0, flds[5], Fg_make_mex(e.co2));
  mxSetField(me, 0, flds[6], Mbe_make_mex(e.mbe));
  mxSetField(me, 0, flds[7], mxCreateDoubleScalar(e.dbg));
  return me;
}

inline const double* get_darray (const mxArray* m, const char* fld) {
  const mxArray* mf = mxGetField(m, 0, fld);
  if (!mf) return NULL;
  return mxGetPr(mf);
}

inline double* set_dfield(mxArray* m, const char* fld, mxArray* f) {
  mxSetField(m, 0, fld, f);
  return mxGetPr(f);
}

mxArray* add_measurement_fields (const Model& m, const mxArray* de) {
  if (!m.inited) panic("m was not inited.");

  mxArray* dem = mxDuplicateArray(de);

  const double* z = get_darray(de, "z");
  if (!z) return dem;
  const size_t nz = mxGetNumberOfElements(mxGetField(de, 0, "z"));
  const double* y = get_darray(de, "y");
  const double* state = get_darray(de, "state");

  const char* flds[] = {"p", "phi_g", "v", "phi_l", "F", "vvfrac"};
  for (size_t i = 0; i < sizeof(flds)/sizeof(*flds); i++)
    mxAddField(dem, flds[i]);
  double* p     = set_dfield(dem, flds[0], mxCreateDoubleMatrix(nz, 1, mxREAL));
  double* phi_g = set_dfield(dem, flds[1], mxCreateDoubleMatrix(nz, 1, mxREAL));
  double* v     = set_dfield(dem, flds[2], mxCreateDoubleMatrix(nz, 1, mxREAL));
  double* phi_l = set_dfield(dem, flds[3], mxCreateDoubleMatrix(nz, 1, mxREAL));
  double* F     = set_dfield(dem, flds[4], mxCreateDoubleMatrix(nz, 1, mxREAL));
  double* vvfrac = set_dfield(
    dem, flds[5], mxCreateDoubleMatrix(
      m.model == Model::newtonian ? nz : 0, 1, mxREAL));

  vector<Exprs::Equations> eqs;
  if (m.model == Model::newtonian) eqs.push_back(Exprs::emb);
  for (size_t i = 0; i < nz; i++) {
    double yi[4];
    memcpy(yi, y, 4*sizeof(double));
    put_in_bounds(m, yi);
    State s;
    s.k = static_cast<int>(state[0]);
    s.plug = static_cast<int>(state[1]);
    Measurements ms;
    Exprs e;
    calc_exprs(m, s, z[i], yi, NULL, eqs, e, &ms);
    ms.get(p[i], phi_g[i], v[i], phi_l[i], F[i]);
    if (m.model == Model::newtonian) vvfrac[i] = e.mbe.vv/(e.mbe.vv + e.mbe.vf);
    state += 2;
    y += 4;
  }

  return dem;
}

static Model g_m;

void mexFunction (int nlhs, mxArray** plhs, int nrhs, const mxArray** prhs) {
#define cmpfn(str) (fn == str)

  if (nrhs < 1) mexErrMsgTxt("smfmex(cmd, ...)");
  string fn;
  mex_to_string(prhs[0], fn);

  if (cmpfn("ppval")) {
    if (nlhs != 1 || nrhs != 3) mexErrMsgTxt("y = ppval(pp, x)");
    CubicSpline cs;
    CubicSpline_init(cs, prhs[1]);
    const size_t nx = mxGetNumberOfElements(prhs[2]);
    const double* x = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]), mxGetN(prhs[2]), mxREAL);
    double* y = mxGetPr(plhs[0]);
    for (size_t i = 0; i < nx; i++) y[i] = cs.call(x[i]);
  } else if (cmpfn("phi_s_of_p")) {
    if (nlhs != 1 || nrhs != 3)
      mexErrMsgTxt("y = phi_s_of_p(pp, p)");
    CubicSpline cs;
    CubicSpline_init(cs, prhs[1]);
    const size_t np = mxGetNumberOfElements(prhs[2]);
    const double* p = mxGetPr(prhs[2]);
    plhs[0] = mxCreateDoubleMatrix(mxGetM(prhs[2]), mxGetN(prhs[2]), mxREAL);
    double* chi = mxGetPr(plhs[0]);
    for (size_t i = 0; i < np; i++) chi[i] = calc_phi_s_of_p(cs, p[i]);
  } else if (cmpfn("m_init")) {
    if (nlhs != 0 || nrhs != 2) mexErrMsgTxt("smfmex('m_init', m)");
    Model_init(g_m, prhs[1]);
  } else if (cmpfn("calc_exprs")) {
    if (nlhs != 1 || nrhs != 6)
      mexErrMsgTxt("c = calc_exprs(state, z, y, yp, eqs)");
    State state;
    State_init(state, prhs[1]);
    const double z = mxGetScalar(prhs[2]);
    const double* y = mxGetPr(prhs[3]);
    const double* yp = NULL;
    if (mxGetNumberOfElements(prhs[4]) > 0) yp = mxGetPr(prhs[4]);
    vector<Exprs::Equations> eqs;
    Equations_init(eqs, prhs[5]);
    Exprs e;
    calc_exprs(g_m, state, z, y, yp, eqs, e);
    plhs[0] = Exprs_make_mex(e);
  } else if (cmpfn("calc_f_y")) {
    if (nlhs != 1 || nrhs != 6)
      mexErrMsgTxt("f_y = calc_f_y(state, z, y, yp, f)");
    State state;
    State_init(state, prhs[1]);
    const double z = mxGetScalar(prhs[2]);
    const double* y = mxGetPr(prhs[3]);
    const double* yp = NULL;
    const double* f = mxGetPr(prhs[5]);
    plhs[0] = mxCreateDoubleMatrix(3, 4, mxREAL);
    calc_f_y(g_m, state, z, y, f, mxGetPr(plhs[0]));
  } else if (cmpfn("calc_F")) {
    if (nlhs != 1 || nrhs != 5) mexErrMsgTxt("f_y = calc_F(state, z, y, yp)");
    State state;
    State_init(state, prhs[1]);
    const double z = mxGetScalar(prhs[2]);
    const double* y = mxGetPr(prhs[3]);
    const double* yp = mxGetPr(prhs[4]);
    plhs[0] = mxCreateDoubleMatrix(4, 1, mxREAL);
    calc_F(g_m, state, z, y, yp, mxGetPr(plhs[0]));
  } else if (cmpfn("mb_solve")) {
    if (nlhs != 4 || nrhs != 12)
      mexErrMsgTxt(
        "[x f F zeta] = mb_solve(a, b, c, d, g, h, k, p, tolx, tolf, x)");
    double c[11], x, F, zeta, f;
    for (size_t i = 1; i < 12; i++) c[i-1] = mxGetScalar(prhs[i]);
    mb::solve(c[0], c[1], c[2], c[3], c[4], c[5], c[6], c[7], c[8], c[9],
              c[10], x, F, zeta, f);
    plhs[0] = mxCreateDoubleScalar(x);
    plhs[1] = mxCreateDoubleScalar(f);
    plhs[2] = mxCreateDoubleScalar(F);
    plhs[3] = mxCreateDoubleScalar(zeta);
  } else if (cmpfn("add_measurement_fields")) {
    if (nlhs != 1 || nrhs != 2) mexErrMsgTxt("de = add_measurement_fields(de)");
    plhs[0] = add_measurement_fields(g_m, prhs[1]);
  } else {
    mexErrMsgTxt((string("Invalid function: ") + fn).c_str());
  }

#undef cmpfn
}
#else
// Tests.

void panic (const string& err) {
  fprintf(stderr, "%s\n", err.c_str());
  exit(-1);
}

int main (int argc, char** argv) {
  const double breaks[] = {-1.1564, -9.7921e-01, -8.3137e-01, -5.3356e-01};
  const double coefs[] = {
    -1.9039e+02, -2.7021e+01, 2.7509e+01, -2.0026, 1.6848e+02, -4.5231e+01,
    0, 9.6423e-01, 3.1914, 7.7132e-01, -2.3263, 5.2006e-01};
  CubicSpline cs(breaks, coefs, 3);
  double xl = cs.get_x_low(), yl = cs.get_y_low(), xh = cs.get_x_high(),
    yh = cs.get_y_high();
  double x = breaks[0] - 0.2, dx = 0.05;
  double y;
  while (x < breaks[3] + 0.2) {
    y = cs.call(x);
    x += dx;
  }
}
#endif
