// vim: set ai cin et sw=2 ts=2 tw=80:
#pragma once

#include <cmath>

#include "aux.hpp"

struct cross_sec {
  virtual double operator()(const double sminus, const double splus) const = 0;
  virtual ~cross_sec() {}
};

struct scalar_med_s : cross_sec {
  scalar_med_s(const double coupling, const double med_mass,
               const double dec_wid)
      : g{coupling}, mediator_mass{med_mass}, decay_width{dec_wid} {}

  double g;
  double mediator_mass;
  double decay_width;
  double operator()(const double sminus, const double splus) const override {
    double Gamma_s{0.0};
    if (splus < 1e-5) // We Taylor-expand atandiff to avoid roundoff errors
    {
      Gamma_s =
          SQR(SQR(g)) /
          (32 * M_PI * SQR(this->mediator_mass) * this->decay_width) *
          (2 * this->mediator_mass *
               ((this->decay_width / this->mediator_mass *
                 (1 + SQR(this->decay_width / this->mediator_mass) +
                  2 * sminus)) /
                    SQR(1 + SQR(this->decay_width / this->mediator_mass)) *
                    (splus - sminus) +
                (this->decay_width / this->mediator_mass) /
                    SQR(1 + SQR(this->decay_width / this->mediator_mass)) *
                    SQR(splus - sminus)) +
           this->decay_width *
               (log1p(SQR(this->mediator_mass) /
                      (SQR(this->mediator_mass) + SQR(this->decay_width)) *
                      splus * (splus - 2)) -
                log1p(SQR(this->mediator_mass) /
                      (SQR(this->mediator_mass) + SQR(this->decay_width)) *
                      sminus * (sminus - 2))));
    } else {
      Gamma_s =
          SQR(SQR(g)) /
          (32 * M_PI * SQR(this->mediator_mass) * this->decay_width) *
          (2 * this->mediator_mass *
               nuSIaux::atandiff(
                   this->mediator_mass * (splus - 1) / this->decay_width,
                   this->mediator_mass * (sminus - 1) / this->decay_width) +
           this->decay_width *
               (log1p(SQR(this->mediator_mass) /
                      (SQR(this->mediator_mass) + SQR(this->decay_width)) *
                      splus * (splus - 2)) -
                log1p(SQR(this->mediator_mass) /
                      (SQR(this->mediator_mass) + SQR(this->decay_width)) *
                      sminus * (sminus - 2))));
    }
    return Gamma_s;
  }
};
