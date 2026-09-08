#include "nuSIprop.hpp"
#include <fmt/format.h>
#include <fmt/os.h>
#include <ranges>
#include <string>

const double MEDIATOR_MASS{1.0E6}; // [eV]
const double COUPLING{1.0E-2};
const double SPECTRAL_INDEX{2.2};
const bool IS_NORMAL_ORDERED{true};
// "./results/res_inv_mKK001_g01_E2_70_abc.dat"
const std::string FLUX_OUTPUT_FNAME{
  fmt::format("./results/res_{:s}_mKK{:03d}_g{:02d}_E{:s}_abc.dat",
      IS_NORMAL_ORDERED ? "nor" : "inv",
      static_cast<int>(MEDIATOR_MASS / 1.0E5),
      static_cast<int>(COUPLING * 100),
      [](){
        std::string s = fmt::format("{:.2f}", SPECTRAL_INDEX);
        std::replace(s.begin(), s.end(), '.', '_');
        return s;
      }())
};
auto out = fmt::output_file(FLUX_OUTPUT_FNAME);

int main() {
  // Construct evolver object
  nuSIprop::calculate_flux evolver(
      MEDIATOR_MASS,  // Mediator mass [eV]
      COUPLING, // Coupling
      0.12,  // Sum of neutrino masses [eV]
      SPECTRAL_INDEX,  // Spectral index
      6,    // Normalization of the free-streaming flux at 100 TeV [Default = 1]
      false, // Majorana neutrinos? [Default = true]
      true, // Include non s-channel contributions? Relevant for couplings
            // g>~0.1 [Default = true]
      IS_NORMAL_ORDERED, // Normal neutrino mass ordering? [Default = true]
      1200, // Number of energy bins, uniformly distributed in log space [Default
           // = 300]
      13,  // log_10 (E_min/eV) [Default = 13]
      16,  // log_10 (E_max/eV) [Default = 17]
      5,   // Largest redshift at which neutrino sources are included [Default =
           // 5]
      2,   // Flavor of interacting neutrinos [0=e, 1=mu, 2=tau. Default = 2]
      false // Consider double-scalar production? If set to true, the files
            // xsec/alpha_phiphi.bin and xsec/alphatilde_phiphi.bin must exist
            // [Default = false]
  );

  // Evolve it
  KKMuTau xs{evolver.g, evolver.mphi, 0.0};
  evolver.set_cross_section(xs);
  evolver.evolve();

  // Output the result
  out.print("#Energy[eV]  nu_e flux   nu_mu flux  nu_tau flux\n");
  for (const auto & i: std::views::iota(0, evolver.get_N_bins_E()))
    out.print("{:.5e}  {:.4e}  {:.4e}  {:.4e}\n",
           evolver.get_energy_min(i),  // Energy bin minimum
           evolver.get_flux_fla(0, i), // nu_e flux
           evolver.get_flux_fla(1, i), // nu_mu flux
           evolver.get_flux_fla(2, i)  // nu_tau flux
    );

  /* Other usage tips:

     -If you want the flux in the mass basis, just call
      evolver.get_flux(i, j), where i={0,1,2} is the mass eigenstate
      and j is the energy bin. Neutrino mass eigenstates are ordered
      in the usual convention: nu_1 (nu_3) is the eigenstate with the smallest
      (largest) electron neutrino content.

     -You can check that the total neutrino energy is conserved by
     checking the output of evolver.check_energy_conservation(). This
     returns the relative difference in total energy with and without
     neutrino self-interactions.  It thus gives an idea of the
     size of numerical errors.

     -The public member variables
        evolver.mphi
        evolver.g
        evolver.mntot
        evolver.si
        evolver.norm
      can be changed between consecutive runs. They correspond to the
      mediator mass [eV], the coupling, the total neutrino mass [eV],
      the spectral index of the propagated flux, and its normalization.
      Remember calling evolve() after updating them to obtain the correct
      flux.
  */
}
