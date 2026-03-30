// vim: set ai et ts=4 sw=4 tw=80:
//

#pragma once

// Note on consistent units:
// system: PDG natural units
// mass: MeV

#include "misc.hpp"
#include <complex>
#include <fmt/format.h>

inline double s(const double E_nu,   // MeV
                const double mass_nu // MeV
) {
    // units: MeV^2
    return 2.0 * E_nu * mass_nu;
}

enum class MassHierarchy { NORMAL = 0, INVERTED = 1 };
enum struct ScatterChannels : int {
    NONE = -1,
    ALL = 0,
    S = 2,
    T = 3,
    U = 5,
};

class Model {
  public:
    explicit Model(const arr3 &mass_nu, // ordered by increasing masses
                   const MassHierarchy order = MassHierarchy::NORMAL)
        : m_nu{mass_nu}, mass_order{order} {}

    inline double mass_nu1() {
        return (mass_order == MassHierarchy::NORMAL ? m_nu[0] : m_nu[1]) /
               MeV_TO_eV;
    }

    inline double mass_nu2() {
        return (mass_order == MassHierarchy::NORMAL ? m_nu[1] : m_nu[2]) /
               MeV_TO_eV;
    }

    inline double mass_nu3() {
        return (mass_order == MassHierarchy::NORMAL ? m_nu[2] : m_nu[0]) /
               MeV_TO_eV;
    }

    inline MassHierarchy nu_mass_hierarchy() { return mass_order; }
    virtual double sigma(const double, const ScatterChannels) = 0;
    virtual double gamma_2(const double) = 0;

  private:
    arr3 m_nu;
    MassHierarchy mass_order;
};

class KKModel : public Model {
  public:
    KKModel(const double g, // coupling constant
            const double m_KK, const double y_sm, const arr3 &mass_nu,
            const MassHierarchy order = MassHierarchy::NORMAL)
        : Model(mass_nu, order), gD{g}, m_0{m_KK}, y_SM{y_sm} {}

    inline double mass_zero() { return m_0; }
    // Decay width to two final states
    double gamma_2(const double m_ini) override;
    double sigma(const double y,
                 const ScatterChannels chan = ScatterChannels::ALL) override;

  private:
    /* data */
    static const unsigned int N_MODES;
    const double gD;
    double m_0;
    double y_SM;

    double cos_sqr_factor(const double m_ini);
    // returns scaled mass_n = mass_Zn / mass_Z0
    inline double mass_Zn_factor(const unsigned int n) { return (2 * n - 1); }
    inline double mass_Zn(const unsigned int n) {
        return mass_Zn_factor(n) * m_0;
    }
    std::complex<double>
    amp_n(const ScatterChannels ch, // interaction channel: s, t, or u
          const double y,           // y = s/mZ_0^2
          const double cos_th,      // cos(theta)
          const unsigned int n      // n-th KK mode
    );
    std::complex<double>
    sum_amp(const double y, const double cos_th,
            const ScatterChannels schan = ScatterChannels::ALL);
};

[[maybe_unused]] const arr3 MASS_NU_NOR_V61 = {0.02984, 0.03109, 0.059056};
[[maybe_unused]] const arr3 MASS_NU_INV_V61 = {0.014668, 0.052308, 0.053024};
