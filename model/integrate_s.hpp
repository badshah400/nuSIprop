// vim: set ai et ts=4 sw=4 tw=80:
//

#pragma once
/*
 * Helper class for integrating cross sections between specified initial and
 * final energies
 * 
 */

#include "model/misc.hpp"
#include <fmt/format.h>
#include <gsl/gsl_integration.h>
#include <tuple>

template <typename T> class IntegrateS {
  public:
    explicit IntegrateS(const T &t) : model{t}, err_rel{1.0E-06}, err_abs{0.0} {}
    explicit IntegrateS(const IntegrateS<T> & IG):
        model{IG.model}, err_rel{IG.err_rel}, err_abs{IG.err_abs}
    {}
    IntegrateS<T>& operator=(const IntegrateS<T> & IG) {
        model = IG.model;
        err_rel = IG.err_rel;
        err_abs = IG.err_abs;
        return *this;
    }

    ResultType operator()(const double s_minus, const double s_plus) const {
        // fmt::print("{:.5E}", s_minus);
        gsl_function F;
        F.function = [](double x, void *par) {
            // auto *ig_par = static_cast<
            //     std::tuple<T *, const double, const double, const double> *>(
            //     par);
            // auto [_mod, _s_min, _s_max, _eng_res] = *ig_par;
            // auto _s_mean = 0.5 * (_s_min + _s_max);
            // return _mod->sigma(x) *
            //        exp(-0.5 * sqr(x - _s_mean) / sqr(_eng_res * _s_mean));
            auto _mod = static_cast<T *>(par);
            return _mod->sigma(x);
            // auto ig_par = static_cast<
            //     std::tuple<T *, double, double, double> *>(par);
            // auto _mod = std::get<0>(*ig_par);
            // auto _s_mean = 0.5 * (std::get<1>(*ig_par) + std::get<2>(*ig_par));
            // auto _eng_res = std::get<3>(*ig_par);
            // auto res = _mod->model->sigma(x) /
            //        exp(-0.5 * sqr(x - _s_mean) / sqr(_eng_res * _s_mean));
            // fmt::print("{:.5E}\n", res);
            // return res;
        };
        // auto pars = std::make_tuple(this, s_minus, s_plus, ENG_RESOLUTION);
        auto pars = this->model;
        // if (abs(s_minus - 1.0) < 1.0E-03) {
        //     fmt::print("{:15.5E}\n", s_minus);
        // }
        F.params = &pars;
        double res{0.0}, err{0.0};
        // size_t n_eval{0};
        gsl_integration_workspace *wspace =
            gsl_integration_workspace_alloc(N_SPACE + 1);
        gsl_integration_qag(&F,
                             s_minus,
                             s_plus,
                             err_abs,
                             err_rel,
                             N_SPACE,
                             GSL_INTEG_GAUSS21,
                             wspace,
                             &res,
                             &err);
        gsl_integration_workspace_free(wspace);
        return std::make_tuple(res, err);
    }

  private:
    T model;
    static constexpr unsigned int N_SPACE = 1000;
    double err_rel;
    double err_abs;
    static constexpr double ENG_RESOLUTION = 0.15;
};
