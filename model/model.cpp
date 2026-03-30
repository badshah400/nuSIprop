// vim: set ai et ts=4 sw=4 tw=80:
//

// Note on consistent units:
// system: PDG natural units
// mass: MeV

#include "model/model.hpp"
#include <gsl/gsl_integration.h>
#include <stdexcept>

const unsigned int KKModel::N_MODES{10};

double KKModel::cos_sqr_factor(const double m_ini) {
    return sqr(cos(m_ini * y_SM));
}

double KKModel::gamma_2(const double m_ini) {
    return sqr(gD) * cos_sqr_factor(m_ini) * m_ini / 12.0 / pi;
}

std::complex<double> KKModel::amp_n(const ScatterChannels ch, const double y,
                                    const double cos_theta,
                                    const unsigned int n) {
    using namespace std::complex_literals;
    auto scaled_t = -0.5 * y * (1 - cos_theta);
    auto scaled_u = -0.5 * y * (1 + cos_theta);
    // |M|^2 numerator is 2*u^2 consistently for every term,
    // so we multiply (sqrt(2) * u) with the amplitude for each term
    auto res = cos_sqr_factor(n) * sqrt(2.0) * scaled_u;
    switch (ch) {
    case ScatterChannels::S:
        // mZ/m0 = 2n-1
        return res / (y - sqr(mass_Zn(n) / m_0) +
                      (1.0i * gamma_2(mass_Zn(n)) * mass_Zn_factor(n)));
        break;
    case ScatterChannels::T:
        return res / (scaled_t - sqr(mass_Zn_factor(n)));
        break;
    case ScatterChannels::U:
        return res / (scaled_u - sqr(mass_Zn_factor(n)));
        break;
    default:
        return 0.0;
    }
}

std::complex<double>
KKModel::sum_amp(const double y, const double cos_th,
                 [[maybe_unused]] const ScatterChannels chan) {
    std::complex<double> res{0.0, 0.0};
    for (unsigned int i = 1; i <= KKModel::N_MODES; ++i) {
        // std::complex<double> amp_i{0.0, 0.0};
        res += amp_n(ScatterChannels::S, y, cos_th, i) +
               amp_n(ScatterChannels::T, y, cos_th, i) +
               amp_n(ScatterChannels::U, y, cos_th, i);
    }
    return res;
}

double KKModel::sigma(const double y, const ScatterChannels chan) {
    // Compute
    //             /
    //             | |M|^2 d(cos_theta)
    //             /
    double res{0.0}, err{0.0};

    // integrand: use capture-less lambda for gsl C-type function
    auto ig_func = [](double x, void *par) {
        auto *mypars =
            static_cast<std::tuple<KKModel *, double, ScatterChannels> *>(par);
        auto obj = std::get<0>(*mypars);
        auto amp = obj->sum_amp(std::get<1>(*mypars), x, std::get<2>(*mypars));
        return std::norm(amp);
    };
    auto gsl_params = std::make_tuple(this, y, chan);

    gsl_function F;
    F.function = ig_func;
    F.params = &gsl_params;
    const auto cos_limits = std::make_tuple(-1.0f, 1.0f);
    const size_t GSL_WSPACE_LIMIT{1'000};
    const double GSL_ABS_ERROR{0.0};
    const double GSL_REL_ERROR{1.0E-06};
    gsl_integration_workspace *w =
        gsl_integration_workspace_alloc(GSL_WSPACE_LIMIT);
    auto status = gsl_integration_qag(&F,
                                      std::get<0>(cos_limits),
                                      std::get<1>(cos_limits),
                                      GSL_ABS_ERROR,
                                      GSL_REL_ERROR,
                                      GSL_WSPACE_LIMIT,
                                      GSL_INTEG_GAUSS41,
                                      w,
                                      &res,
                                      &err);
    gsl_integration_workspace_free(w);
    if (status) {
        throw std::runtime_error{
            "Invalid return status from gsl integration."};
    }

    //            /
    //    sigma = | |M|^2 d(cos_theta) / (32*pi*s)
    //            /
    return sqr(sqr(gD)) * res / y / 32.0 / std::numbers::pi / sqr(m_0);
}
