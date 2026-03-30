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

typedef std::tuple<double, size_t> ResultType;

template <typename T> class IntegrateS {
  public:
    explicit IntegrateS(T &t)
        : model{t}, err_rel{1.0E-06}, err_abs{0.0} {

          };

    ResultType operator()(const double s_minus, const double s_plus) {
        gsl_function F;
        F.function = [](double x, void *par) {
            auto _mod = static_cast<T *>(par);
            return _mod->sigma(x);
        };
        F.params = this;
        double res{0.0}, err{0.0};
        size_t n_eval{1000};
        gsl_integration_workspace *wspace =
            gsl_integration_workspace_alloc(N_SPACE);
        gsl_integration_qag(&F,
                            s_minus,
                            s_plus,
                            err_abs,
                            err_rel,
                            n_eval,
                            GSL_INTEG_GAUSS15,
                            wspace,
                            &res,
                            &err);
        gsl_integration_workspace_free(wspace);
        return std::make_tuple(res, n_eval);
    }

    virtual ~IntegrateS() {};

  protected:
    gsl_function integrand;

  private:
    T model;
    static const unsigned int N_SPACE = 1001;
    double err_rel;
    double err_abs;
};
