// vim: set ai et ts=4 sw=4 tw=80:
//

// Note on consistent units:
// system: PDG natural units
// mass: MeV

#include "model.hpp"
#include <algorithm>
#include <gsl/gsl_integration.h>
#include <ranges>
#include <stdexcept>

double
KKModel::cos_sqr_factor(const double m_ini)
{
  return sqr(cos(m_ini * y_SM));
}

double
KKModel::gamma_2(const double m_ini)
{
  return sqr(gD) * cos_sqr_factor(m_ini) * m_ini / 12.0 / pi;
}

double
KKModel::dsigma_dcos(const double y, const double cos_th)
{
  ///
  /// Return the |M^2| as a function of scaled s (= `y`) and cos(theta) =
  /// `cos_th`, where theta is the scattering angle.
  ///
  /// |M^2| is plugged in from symbolic calculations as a product of
  /// amplitudes:
  /// |M^2| = |\sum_{i=1}^{N_MODES} \sum_{j=1}^{N_MODES} x
  ///          (S_i * conj(S_j) + T_i * conj(T_j) + S_i * conj(T_j)...) |,
  /// where S_i is the s-channel amplitude for the i-th KK mode and so on.
  ///
  using namespace std::complex_literals;
  auto scaled_t = -0.5 * y * (1 - cos_th); // t / m_0^2
  auto scaled_u = -0.5 * y * (1 + cos_th); // u / m_0^2
  auto amp_s_nth = [this, &y](const size_t n) {
    auto num = cos_sqr_factor(n) * sqrt(2.0);
    return num / (y - sqr(mass_Zn_factor(n)) -
                  1.0i * gamma_2(mass_Zn(n)) * mass_Zn(n) / sqr(m_0));
  };
  auto amp_t_nth = [this, &scaled_t](const size_t n) {
    auto num = cos_sqr_factor(n) * sqrt(2.0);
    return num / (scaled_t - sqr(mass_Zn_factor(n)));
  };
  auto amp_u_nth = [this, &scaled_u](const size_t n) {
    auto num = cos_sqr_factor(n) * sqrt(2.0);
    return num / (scaled_u - sqr(mass_Zn_factor(n)));
  };

  std::complex<double> nu_anu_ss{ 0.0i }, nu_anu_st{ 0.0i }, nu_anu_tt{ 0.0i };
  std::complex<double> nu_nu_uu{ 0.0i }, nu_nu_tu{ 0.0i };
  [[maybe_unused]] std::array<std::complex<double>, N_MODES> amp_s;
  [[maybe_unused]] std::array<std::complex<double>, N_MODES> amp_t;
  [[maybe_unused]] std::array<std::complex<double>, N_MODES> amp_u;
  std::generate(amp_s.begin(), amp_s.end(), [&amp_s_nth, N = 0u]() mutable {
    return amp_s_nth(++N);
  });
  std::generate(amp_t.begin(), amp_t.end(), [&amp_t_nth, N = 0u]() mutable {
    return amp_t_nth(++N);
  });
  std::generate(amp_u.begin(), amp_u.end(), [&amp_u_nth, N = 0u]() mutable {
    return amp_u_nth(++N);
  });

  for (auto i : std::views::iota(1u, N_MODES + 1)) {
    for (auto j : std::views::iota(1u, N_MODES + 1)) {
      nu_anu_ss += (amp_s[i - 1] * std::conj(amp_s[j - 1])); // S_i * S_j
      nu_anu_tt += (amp_t[i - 1] * std::conj(amp_t[j - 1])); // T_i * T_j
      nu_anu_st += (amp_s[i - 1] * std::conj(amp_t[j - 1])); // S_i * T_j
      nu_nu_uu += (amp_u[i - 1] * std::conj(amp_u[j - 1]));  // U_i * U_j
      nu_nu_tu += (amp_u[i - 1] * std::conj(amp_t[j - 1]));  // U_i * T_j
    }
  }

  // Additional factor of 0.5 for nu-nu identical final states
  return 2.0 * // \sum |U^2_{mu l} - U^2_{tau k}|
         std::abs((sqr(scaled_u) * nu_anu_ss) +       // nu-anu
                  (1.5 * sqr(scaled_u) * nu_anu_tt) + // nu-anu + nu-nu
                  (0.5 * sqr(y) * nu_nu_uu) +         // nu-nu
                  (2.0 * sqr(scaled_u) * std::real(nu_anu_st)) + // nu-anu
                  (2.0 * 0.5 * sqr(y) * std::real(nu_nu_tu))     // nu-nu
         );
}

double
KKModel::sigma(const double y, [[maybe_unused]] const ScatterChannels chan)
{
  ///
  /// Compute
  ///             /
  ///             | |M|^2 d(cos_theta)
  ///             /
  ///
  double res{ 0.0 }, err{ 0.0 };

  ///
  /// integrand: use capture-less lambda for gsl C-style function. `x` is the
  /// cos_theta and `y` is obtained from the par instead.
  ///
  /// Note: we need to pass everything via the `void * par`, which is awkward,
  /// but no way out since captured lambdas are a pain to transform to C-style
  /// functions used by GSL.
  ///
  auto ig_func = [](double x, void* par) {
    auto* mypars = static_cast<std::tuple<KKModel*, const double>*>(par);
    auto obj = std::get<0>(*mypars);
    return obj->dsigma_dcos(std::get<1>(*mypars), x);
  };
  auto gsl_params = std::make_tuple(this, y);

  // Integration bull-work starts
  gsl_function F;
  F.function = ig_func;
  F.params = &gsl_params;
  const auto cos_limits = std::make_tuple(-1.0f, 1.0f);
  const size_t GSL_WSPACE_LIMIT{ 1'000 };
  const double GSL_ABS_ERROR{ 0.0 };
  const double GSL_REL_ERROR{ 1.0E-06 };
  gsl_integration_workspace* w =
    gsl_integration_workspace_alloc(GSL_WSPACE_LIMIT);
  auto status = gsl_integration_qag(&F,
                                    std::get<0>(cos_limits),
                                    std::get<1>(cos_limits),
                                    GSL_ABS_ERROR,
                                    GSL_REL_ERROR,
                                    GSL_WSPACE_LIMIT,
                                    GSL_INTEG_GAUSS21,
                                    w,
                                    &res,
                                    &err);
  gsl_integration_workspace_free(w);
  if (status) {
    throw std::runtime_error{ "Invalid return status from gsl integration." };
  }
  ///                           /
  /// Integration ends; `res` = | |M|^2 d(cos_theta) in CM frame if successful
  ///                           /
  ///

  ///
  ///                  /
  ///    sigma = g^4 * | |M|^2 d(cos_theta) / (32*pi*s)
  ///                  /
  ///
  return sqr(sqr(gD)) * res / y / 32.0 / pi / sqr(m_0);
}
