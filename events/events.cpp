// vim: set ai et ts=2 sw=2 tw=80:
//

///
/// This computes the event rates for a specified amount of IceCube run-time
/// given bin-by-bin effective areas and flux data.
///

#include "flux_integration.hpp"
#include "tabdata.hpp"
#include <algorithm>
#include <execution>
#include <fmt/format.h>
#include <fstream>
#include <functional>
#include <map>
#include <numeric>
#include <ranges>
#include <regex>
#include <string>
#include <vector>

const double IC_ENG_THRES{ 6.0E04 }; // 60 TeV threshold
const std::string FLUX_FILENAME{
  "./results/res_nor_mKK010_g01_E2_20_abc.dat"
};
const std::string ATMOS_EVENTS(
  "../IceCube-data-analysis/results/atmospheric_numu.tsv");
const std::string EFF_AREA_FILENAME{
  "../../../IceCube/HESE-7-year-data-release-main/HESE-7-year-data-release/"
  "effective_area_cm2_all.dat"
};

void
print_data(const DataMap& M)
{
  auto iter = M.cbegin();
  auto data = iter->second;
  for (auto i : std::views::iota(0u, data.size())) {
    for (auto it = M.cbegin(); it != M.cend(); it++) {
      auto row_i = it->second;
      fmt::print("{:15.5E}", row_i.at(i));
    }
    fmt::println("");
  }
  return;
}

// add a final bin at 1.0E7 GeV with zero flux for ease when looping for
// integration with zero fluxes for all flavours
void
append_ten_TeV_bin(DataMap& d, const std::vector<std::string>& energy_headers)
{
  for (auto& [k, v] : d) {
    if (std::any_of(energy_headers.begin(),
                    energy_headers.end(),
                    [&k](const auto& s) { return s == k; })) {
      v.push_back(1.0E7);
    } else {
      v.push_back(0.0); // push zero to other columns
    }
  }
}

// get iterator for 60 TeV energy bin; this is the bin previous to where
// `energy > 6.0E4` first happens
inline auto
bin_at_thres_energy(const VecD& energies, const double& thres = IC_ENG_THRES)
{
  auto _it = std::find_if_not(energies.cbegin(),
                              energies.cend(),
                              [&thres](const auto& x) { return x < thres; });
  // unless iterator points to first element, return iterator for previous bin
  return (_it == energies.cbegin()) ? _it : _it - 1;
}

int
main(int argc, const char* argv[])
{
  std::string flux_filename{ FLUX_FILENAME };
  const std::string eff_area_filename{ EFF_AREA_FILENAME };
  if (argc == 2) {
    flux_filename = argv[1];
  }
  const std::string EVENT_OUTPUT_FILENAME{ [&flux_filename]() {
    auto res =
      std::regex_replace(flux_filename, std::regex("res(new)?_"), "Events$1_");
    res = std::regex_replace(res, std::regex("\\.dat"), ".txt");
    return res;
  }() };
  std::ifstream flux_file(flux_filename, std::ios::in);

  auto flux = TabData(flux_filename,
                      { "energy", "flux_nue", "flux_numu", "flux_nutau" });
  auto flux_data = flux.get_data();

  append_ten_TeV_bin(flux_data, { "energy" });
  auto flux_energies = flux_data.at("energy");

  std::ifstream eff_area_file(eff_area_filename, std::ios::in);
  auto eff_area = TabData(
    eff_area_filename, { "energy_min", "energy_max", "nue", "mumu", "nutau" });

  auto eff_area_energies = eff_area.col_data(0);

  auto res = rebin_and_integrate(flux, eff_area_energies);

  VecD events{};
  for (auto i : std::views::iota(0u, res.at(0).size())) {
    events.push_back([&]() {
      double _r{ 0.0 };
      for (auto c : std::views::iota(1u, res.size())) {
        _r += res.at(c).at(i) * eff_area.col_data(c + 1).at(i);
      }
      return _r;
    }());
  }

  // drop events and energies when the latter is less than 10 TeV
  const auto it_threshold =
    std::find_if_not(eff_area_energies.cbegin(),
                     eff_area_energies.cend(),
                     [](const auto& x) { return x < 1.0E4; });

  if (it_threshold != eff_area_energies.cend()) {
    events.erase(events.cbegin(),
                 events.cbegin() +
                   std::distance(eff_area_energies.cbegin(), it_threshold));
    eff_area_energies.erase(
      eff_area_energies.cbegin(),
      eff_area_energies.cbegin() +
        std::distance(eff_area_energies.cbegin(), it_threshold));
  }

  // Atmospheric data above 60 TeV
  auto atmos_data = TabData(ATMOS_EVENTS, { "E_min", "E_max", "Events" });
  auto atmos_energy = atmos_data.col_data("E_min");
  auto atmos_events = atmos_data.col_data("Events");
  atmos_events.push_back(1.0E-03); // manual value for final bin
  auto ATM_ENERGY_THRES_IT = bin_at_thres_energy(atmos_energy);
  auto total_atmos_events =
    std::reduce(std::execution::par,
                atmos_events.cbegin() +
                  std::distance(atmos_energy.cbegin(), ATM_ENERGY_THRES_IT),
                atmos_events.cend(),
                0.0);

  const double ICHESE12_TOTAL_EVENTS{ 92 - total_atmos_events };

  // fmt::println("{:.0f}", total_atmos_events);
  // Bin counting must include energy bin before the point where it exceeds 60
  // TeV. For example, for bins B: {11, 22, 33, 44, 55, 66...}, iterator
  // IC_ENG_THRES_IT will point to 66, but the 55-66 bin will also have events
  // exceeding 60 TeV and needs to be counted. Account for this using:
  // (`iterator` - 1) instead of `iterator`.
  const auto IC_ENG_THRES_IT = bin_at_thres_energy(eff_area_energies);
  const double EVENT_NORM =
    ICHESE12_TOTAL_EVENTS /
    std::reduce(std::execution::par,
                events.cbegin() +
                  std::distance(eff_area_energies.cbegin(), IC_ENG_THRES_IT),
                events.cend());

  eff_area_energies.push_back(1.0E7);
  // normalise events to total signal events (minus background) from IC data
  std::transform(std::execution::par,
                 events.cbegin(),
                 events.cend(),
                 events.begin(),
                 [&EVENT_NORM](const auto& x) { return x * EVENT_NORM; });
  VecD total_events(events.size(), 0.0);
  // add atmosheric events to signal normalised events
  std::transform(std::execution::par,
                 events.cbegin(),
                 events.cend(),
                 atmos_events.cbegin(),
                 total_events.begin(),
                 std::plus<>{});

  std::FILE* output = std::fopen(EVENT_OUTPUT_FILENAME.c_str(), "w");
  fmt::println(output,
               "{}",
               "E_min/GeV  E_max/GeV  Ast.Events  Atm.Events  Total.Events");
  for (auto i : std::views::iota(0u, events.size())) {
    fmt::println(output,
                 "{:.5E} {:15.5E} {:15.5E} {:15.5E} {:15.5E}",
                 eff_area_energies.at(i),
                 eff_area_energies.at(i + 1),
                 events.at(i),
                 atmos_events.at(i),
                 total_events.at(i));
  }
  std::fclose(output);
  fmt::println("Event output saved to file: {:s}", EVENT_OUTPUT_FILENAME);

  return 0;
}
