// vim: set ai et ts=4 sw=4 tw=80:
//

#pragma once
/*
 * Helper class for integrating cross sections between specified initial and
 * final energies
 */

#include "model/misc.hpp"
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
    }

    ResultType operator()(const double s_minus, const double s_plus) const {
        gsl_function F;
        F.function = [](double x, void *par) {
            auto _mod = static_cast<T *>(par);
            return _mod->sigma(x);
        };
        auto pars = this->model;
        F.params = &pars;
        double res{0.0}, err{0.0};
        gsl_integration_workspace *wspace =
            gsl_integration_workspace_alloc(N_SPACE + 1);
        gsl_integration_qags(&F,
                             s_minus,
                             s_plus,
                             err_abs,
                             err_rel,
                             N_SPACE,
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
