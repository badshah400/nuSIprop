// vim: set et ai ts=2 sw=2 tw=80:
#include "aux.hpp"
#include <cmath>
#include <complex>
#include <array>

static inline double sqr(const double x) {return x*x;}

// used for mixing matrix U_PMNS
typedef std::array<std::array<std::complex<double>, 3> ,3> cmat_3x3;

// parameter central value and its 1/3 sigma uncertainties organised as:
// {central value, -1sigma, +1sigma, -3sigma, +3sigma}
enum ValPos {CENTRAL=0, MINUS_1SIG, PLUS_1SIG, MINUS_3SIG, PLUS_3SIG};
typedef std::array<double, 5> ParamErr;

// convert degrees to radians
inline double deg_to_rad(const double x) { return x * (M_PI/180.0); }

inline double cos_sqr(const double sin_sqr) { return 1.0 - sin_sqr; }

class mixing_params
{
  /**
   * This class stores and allows access to best-fit parameters describing
   * three-flavour neutrino oscillation and their one/three sigma
   * limits. Default values are taken from Nu-Fit v6.1 [JHEP 12 (2024) 216] and
   * apply to normal ordered masses.
   */
public:
  mixing_params():
      ssq_th12{0.3088, 0.3088 - 0.0066, 0.3088 + 0.0067, 0.2893, 0.3295},
      ssq_th23{0.470,  0.470 - 0.014, 0.470 + 0.017, 0.432, 0.587},
      ssq_th13{0.02249, 0.02249 - 5.7E-04, 0.02249 + 5.7E-04, 0.0207, 0.02420},
      dcp{
        deg_to_rad(207),
        deg_to_rad(207-20),
        deg_to_rad(207+23),
        deg_to_rad(114),
        deg_to_rad(405)
      },
      dmsq_21{
        7.537E-05,
        7.537E-05 - 1.0E-06,
        7.537E-05 + 9.4E-07,
        7.236E-05,
        7.823E-05
      },
      dmsq_31{
        2.521E-03,
        2.521E-03 - 1.8E-05,
        2.521E-03 + 2.6E-05,
        2.454E-03, 2.592E-03
      },
      normal_hierarchy{true}
  {
    //  Default constructor uses best-fit values for normal ordered masses
    set_init_masses();
    return;
  }

  mixing_params(
    const ParamErr & s2_12,
    const ParamErr & s2_23,
    const ParamErr & s2_13,
    const ParamErr & delcp,
    const ParamErr & dm2_21,
    const ParamErr & dm2_31
  ):
    ssq_th12 {s2_12}, ssq_th23 {s2_23}, ssq_th13 {s2_13}, dcp {delcp},
    dmsq_21 {dm2_21}, dmsq_31 {dm2_31},
    normal_hierarchy{dm2_31[ValPos::CENTRAL] > 0 ? true : false}
  {
    set_init_masses();
    return;
  }

  enum ParName {SSQ_T12=0, SSQ_T23, SSQ_T13, DEL_CP, DEL_MSQ_21, DEL_MSQ_31};

  inline bool is_normal_ordered() {return normal_hierarchy;}

  inline double operator()(const ParName & p, const ValPos & e) const
  {
    double r;
    switch(p)
    {
      case SSQ_T12: r    = ssq_th12[e]; break;
      case SSQ_T23: r    = ssq_th23[e]; break;
      case SSQ_T13: r    = ssq_th13[e]; break;
      case DEL_CP:  r    = dcp[e];      break;
      case DEL_MSQ_21: r = dmsq_21[e];  break;
      case DEL_MSQ_31: r = dmsq_31[e];  break;
    }
    return r;
  }

  inline double central(const ParName & p) const
  {
    return this->operator()(p, ValPos::CENTRAL);
  }

  inline std::array<double, 2> one_sigma_limit(const ParName & p) const
  {
    return {
      this->operator()(p, ValPos::MINUS_1SIG),
      this->operator()(p, ValPos::PLUS_1SIG)
    };

  }

  inline std::array<double, 2> three_sigma_limit(const ParName & p) const
  {
    return {
      this->operator()(p, ValPos::MINUS_3SIG),
      this->operator()(p, ValPos::PLUS_3SIG)
    };

  }

  inline cmat_3x3 u_pmns(const ValPos v = ValPos::CENTRAL) const
  {
    cmat_3x3 U;
    double s12=sqrt(this->operator()(mixing_params::ParName::SSQ_T12, v));
    double c12=sqrt(1.0 - sqr(s12));
    double s13=sqrt(this->operator()(mixing_params::ParName::SSQ_T13, v));
    double c13=sqrt(1.0 - sqr(s13));
    double s23=sqrt(this->operator()(mixing_params::ParName::SSQ_T23, v));
    double c23=sqrt(1.0 - sqr(s23));

    double delcp{this->operator()(mixing_params::ParName::DEL_CP, v)};
    std::complex<double> del(cos(delcp), sin(delcp));
    // Standard leptonic mixing matrix
    U[0][0]=c12*c13;
    U[0][1]=s12*c13;
    U[0][2]=s13*1.0/del;
    U[1][0]=-s12*c23-c12*s23*s13*del;
    U[1][1]=c12*c23-s12*s23*s13*del;
    U[1][2]=s23*c13;
    U[2][0]=s12*s23-c12*c23*s13*del;
    U[2][1]=-c12*s23-s12*c23*s13*del;
    U[2][2]=c23*c13;
    return U;
  }

  void set_masses_from_total(const double mn_tot)
  {
    const double mL{
        nuSIaux::getmL(
            mn_tot, dmsq_21[ValPos::CENTRAL], dmsq_31[ValPos::CENTRAL]
          )
    };

    if (normal_hierarchy)
    {
      mass[0] = mL;
      mass[1] = sqrt(sqr(mass[0]) + dmsq_21[ValPos::CENTRAL]);
      mass[2] = sqrt(sqr(mass[0]) + dmsq_31[ValPos::CENTRAL]);
    } else {
      mass[2] = mL;
      mass[1] = sqrt(sqr(mass[2]) - dmsq_31[ValPos::CENTRAL]);
      mass[0] = sqrt(sqr(mass[1]) + dmsq_21[ValPos::CENTRAL]);
    }
  }

  inline const std::array<double, 3> & masses() const { return mass; }

private:
  ParamErr ssq_th12;
  ParamErr ssq_th23;
  ParamErr ssq_th13;
  ParamErr dcp;     // in rad
  ParamErr dmsq_21; // in eV^2
  ParamErr dmsq_31; // in eV^2
  bool normal_hierarchy;
  std::array<double, 3> mass;

  void set_init_masses()
  {
    if (normal_hierarchy)
    {
      mass[0] = 0.0;
      mass[1] = sqrt(dmsq_21[ValPos::CENTRAL]);
      mass[2] = sqrt(dmsq_31[ValPos::CENTRAL]); 
    } else {
      mass[2] = 0.0;
      mass[1] = sqrt(-1.0 * dmsq_31[ValPos::CENTRAL]);
      mass[0] = sqrt(sqr(mass[1]) + dmsq_21[ValPos::CENTRAL]); 
    }
  }
};

const mixing_params NU_FIT61_NOR {
  // \sin^{2}(\theta_12)
  {0.3088, 0.3088 - 0.0066, 0.3088 + 0.0067, 0.2893, 0.3295},

  // \sin^{2}(\theta_23)
  {0.470, 0.470 - 0.014, 0.470 + 0.017, 0.432, 0.587},

  // \sin^{2}(\theta_13)
  {0.02249, 0.02249 - 5.7E-04, 0.02249 + 5.7E-04, 0.0207, 0.02420},

  // \delta_{CP}
  {deg_to_rad(207), deg_to_rad(207-20), deg_to_rad(207+23), deg_to_rad(114),
    deg_to_rad(405)},

  // \Delta m^{2}_{21} [eV^2]
  {7.537E-05, 7.537E-05 - 1.0E-06, 7.537E-05 + 9.4E-07, 7.236E-05, 7.823E-05},

  // \Delta m^{2}_{31} [eV^2]
  {2.521E-03, 2.521E-03 - 1.8E-05, 2.521E-03 + 2.6E-05, 2.454E-03, 2.592E-03}
};

const mixing_params NU_FIT61_INV {
  // \sin^{2}(\theta_12)
  {0.3088, 0.3088 - 0.0066, 0.3088 + 0.0067, 0.2893, 0.3295},

  // \sin^{2}(\theta_23)
  {0.555, 0.555 - 0.016, 0.555 + 0.013, 0.437, 0.590},

  // \sin^{2}(\theta_13)
  {0.2261, 0.2261 - 5.6E-04, 0.2261 + 5.6E-04, 0.02091, 0.02433},

  // \delta_{CP}
  {deg_to_rad(283), deg_to_rad(283-28), deg_to_rad(283+24), deg_to_rad(202),
    deg_to_rad(347)},

  // \Delta m^{2}_{21} [eV^2]
  {7.537E-05, 7.537E-05 - 1.0E-06, 7.537E-05 + 9.4E-07, 7.236E-05, 7.822E-05},

  // \Delta m^{2}_{31} [eV^2]
  {-2.500E-03, -2.500E-03 - 2.3E-05, -2.500E-03 + 2.4E-05, 2.430E-03, 2.569E-03}
};

