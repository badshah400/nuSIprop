#include "nuSIprop.hpp"
#include <fmt/format.h>
#include <fmt/os.h>
#include <limits>
#include <ranges>
#include <string>

const double COUPLING{ 1.0E-02 };
const bool IS_NORMAL_ORDERED{ true };

int
main(int argc, const char* argv[])
{
  double spectral_index{ 2.2 };
  double medi_mass{ 1.05E6 }; // [eV]
  if (argc > 3) {
    fmt::println(stderr, "Invalid number of arguments passed.");
    return -100;
  }
  if (argc >= 2) {
    try {
      spectral_index = std::stod(argv[1]);
    } catch (...) {
      fmt::println(stderr, "Invalid spectral index passed: {:s}", argv[1]);
      return -1;
    }
  }
  if (argc == 3) {
    try {
      medi_mass = std::stod(argv[2]);
    } catch (...) {
      fmt::println(stderr, "Invalid mediator mass passed: {:s}", argv[2]);
      return -2;
    }
  }

  // string for astro nu flux spectral index gamma: `2.7 -> "2_70"`
  const auto SPEC_IDX_STRING = [&spectral_index]() {
    std::string s = fmt::format("{:.2f}", spectral_index);
    std::replace(s.begin(), s.end(), '.', '_');
    return s;
  }();

  //
  // Output filename example for
  // mKK = 0.1 MeV; g = 0.01; and astro. nu flux spectral index, gamma = 2.70:
  // ./results/res_inv_mKK001_g01_E2_70_abc.dat
  //
  const std::string FLUX_OUTPUT_FNAME{
    std::abs(COUPLING) < std::numeric_limits<double>::epsilon()
      ? fmt::format("./results/res_std_E{:s}_abc.dat", SPEC_IDX_STRING)
      : fmt::format("./results/resnew_{:s}_mKK{:03d}_g{:.0E}_E{:s}_abc.dat",
                    IS_NORMAL_ORDERED ? "nor" : "inv",
                    static_cast<int>(medi_mass / 1.0E5),
                    COUPLING,
                    SPEC_IDX_STRING)
  };

  // Construct evolver object
  nuSIprop::calculate_flux evolver(
    medi_mass,      // Mediator mass [eV]
    COUPLING,       // Coupling
    0.12,           // Sum of neutrino masses [eV]
    spectral_index, // Spectral index
    6,     // Normalization of the free-streaming flux at 100 TeV [Default = 1]
    false, // Majorana neutrinos? [Default = true]
    true,  // Include non s-channel contributions? Relevant for couplings
           // g>~0.1 [Default = true]
    IS_NORMAL_ORDERED, // Normal neutrino mass ordering? [Default = true]
    1200, // Number of energy bins, uniformly distributed in log space
          // [Default = 300]
    13,   // log_10 (E_min/eV) [Default = 13]
    19,   // log_10 (E_max/eV) [Default = 17]
    5,    // Largest redshift at which neutrino sources are included [Default =
          // 5]
    2,    // Flavor of interacting neutrinos [0=e, 1=mu, 2=tau. Default = 2]
    false // Consider double-scalar production? If set to true, the files
          // xsec/alpha_phiphi.bin and xsec/alphatilde_phiphi.bin must exist
          // [Default = false]
  );

  // Evolve it
  KKMuTau xs{ evolver.g, evolver.mphi, 0.0 };
  evolver.set_cross_section(xs);
  evolver.evolve();

  // Output the result
  auto out = fmt::output_file(FLUX_OUTPUT_FNAME);
  out.print("#Energy[eV]  nu_e flux   nu_mu flux  nu_tau flux\n");
  for (const auto& i : std::views::iota(0, evolver.get_N_bins_E()))
    out.print("{:.5e}  {:.4e}  {:.4e}  {:.4e}\n",
              evolver.get_energy_min(i),  // Energy bin minimum
              evolver.get_flux_fla(0, i), // nu_e flux
              evolver.get_flux_fla(1, i), // nu_mu flux
              evolver.get_flux_fla(2, i)  // nu_tau flux
    );

  fmt::println("Neutrino flux output saved to file: {}", FLUX_OUTPUT_FNAME);
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
