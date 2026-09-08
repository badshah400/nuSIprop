// vim: set ai et ts=2 sw=2 tw=80:

#include "flux_integration.hpp"
#include "tabdata.hpp"
#include <algorithm>
#include <array>
#include <vector>

const double GeV_TO_eV{1.0E09};

VecDD rebin_and_integrate(const TabData &t, const VecD &bin_left) {
  VecDD integrated_flux{{}, {}, {}, {}};
  auto kkflux_energies = t.col_data(0); // energies [eV]...
  std::transform(kkflux_energies.begin(), kkflux_energies.end(),
                 kkflux_energies.begin(),
                 [](const auto &x) { return x / GeV_TO_eV; }); // ...now [GeV]
  double prev_bin_eng = bin_left.at(0);
  for (auto eng : bin_left) {
    std::array<double, 3> val{0.0, 0.0, 0.0}; // tmp integrated flavour fluxes

    for (auto nbin_theory : std::views::iota(1u, kkflux_energies.size())) {
      const double e_min{kkflux_energies.at(nbin_theory - 1)};
      const double e_max{kkflux_energies.at(nbin_theory)};
      if (e_min > eng) break;
      if ((e_min < prev_bin_eng) || (e_max > eng)) {
        // fmt::println("{:.5E} {:15.5E} {:15.5E} {:15.5E}",
        //     e_min, prev_bin_eng, e_max, eng);
        continue;
      }
      double del_eng{e_max - e_min};
      for (auto flav : std::views::iota(0u, val.size())) {
        val.at(flav) += 0.5 * del_eng *
                        (t.col_data(flav + 1).at(nbin_theory) +
                         t.col_data(flav + 1).at(nbin_theory - 1));
      }
    }
    integrated_flux.at(0).push_back(eng);
    // fmt::print("{:.5E}", integrated_flux.at(0).back());
    for (auto n : std::views::iota(1u, integrated_flux.size())) {
      integrated_flux.at(n).push_back(val.at(n-1));
      // fmt::print("{:15.5E}", val.at(n-1));
    }
    // fmt::println("");
    prev_bin_eng = eng;
  }
  return integrated_flux;
}
