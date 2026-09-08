// vim: set ai et ts=2 sw=2 tw=80:
//

///
/// This codes computes the chisq difference between hypothesised (model) data
/// and experimental IC data with poissonian error bars
///

#include <fmt/format.h>
#include <fmt/ranges.h>

#include <algorithm>
#include <execution>
#include <numeric>
#include <regex>
#include <string>

#include "tabdata.hpp"

size_t
fact(const size_t n)
{
  size_t res{ 1 };
  if (n == 0 || n == 1) {
    return res;
  }
  for (size_t i{ 2 }; i <= n; ++i) {
    res *= i;
  };
  return res;
}

double
loglike(const double& lam, const double& k)
{
  return lam > 0.0 ? (k * log(lam) - lam - log(fact(std::round(k))))
                   : (k * log(1.0E-06) - log(fact(std::round(k))));
}

int
main(int argc, const char* argv[])
{
  // Model event rate file must be specified on the command line
  if (argc != 2) {
    fmt::println("No model event file specified");
    fmt::println("Usage: {:s} <MODEL_EVENT_FILE>", argv[0]);
    return 1;
  }

  const std::string MODEL_EVENT_FILENAME{ argv[1] };
  std::smatch matched_regex;
  auto match = std::regex_search(
    MODEL_EVENT_FILENAME, matched_regex, std::regex("_E\\d_\\d+_abc"));

  if (!match) {
    fmt::println("Could not find spectral index match in filename {}",
                 MODEL_EVENT_FILENAME);
    return 2;
  }
  std::string _spec_idx = matched_regex.begin()->str();
  auto gamma = std::stof(
    std::regex_replace(_spec_idx, std::regex("_E(\\d+)_(\\d+)_abc"), "$1.$2"));

  float mKK{};
  match = std::regex_search(
    MODEL_EVENT_FILENAME, matched_regex, std::regex("_mKK\\d+_"));
  if (!match) {
    mKK = 0.0;
  } else {
    mKK = std::stof(std::regex_replace(
      matched_regex.begin()->str(), std::regex("_mKK(\\d+)"), "$1"));
    // Convert to MeV
    mKK /= 10.0;
  }

  const TabData MODEL_EVENT_DATA(
    MODEL_EVENT_FILENAME,
    { "e_min", "e_max", "sig_events", "atm_events", "tot_events" },
    false);
  const auto MODEL_TOTAL_EVENTS = MODEL_EVENT_DATA.col_data("tot_events");
  const VecD BIN_EMIN = MODEL_EVENT_DATA.col_data("e_min");

  const std::string ICECUBE_EVENT_ERR_FILENAME{
    "../IceCube-data-analysis/results/IC_12_year_events_68pc_cl.tsv"
  };
  // const auto [IC_EVENTS, IC_EVENT_MINUS, IC_EVENT_PLUS]
  const TabData IC_EVENT_DATA(
    ICECUBE_EVENT_ERR_FILENAME, { "events", "minus", "plus" }, false);
  const VecD IC_EVENTS = IC_EVENT_DATA.col_data("events");
  const VecD IC_EVENT_ERRMINUS = IC_EVENT_DATA.col_data("minus");
  const VecD IC_EVENT_ERRPLUS = IC_EVENT_DATA.col_data("plus");

  ///
  /// Find location of bin where energy exceeds 60 TeV
  ///
  const auto THRESHOLD_BIN = std::distance(
    BIN_EMIN.cbegin(),
    std::find_if(
      BIN_EMIN.cbegin(),
      BIN_EMIN.cend(),
      [](const auto & x) { return x > 6.0E4; })
  );

  ///
  /// symmetrised errors:
  /// For datum with central value X (X > 0), errors dm and dp, new
  /// symmetrised data is given by: X +/- sigma, where. variance sigma =
  /// sqrt(dm * dp)...
  ///
  VecD ic_symm_err(IC_EVENTS.size(), 0.0);
  std::transform(
    std::execution::par,
    IC_EVENT_ERRMINUS.begin(),
    IC_EVENT_ERRMINUS.end(),
    IC_EVENT_ERRPLUS.begin(),
    ic_symm_err.begin(),
    [](const auto& x, const auto& y) { return x > 0? std::min(x, y): y; });
  ///
  /// ...If sigma > X, then we enforce sigma = X. This is the case, for
  /// example, if X=1, where dm=0.63 and dp=1.75 => sigma=1.05 > 1 (=X) =>
  /// sigma=1
  ///
  std::transform(std::execution::par,
                 ic_symm_err.begin(),
                 ic_symm_err.end(),
                 IC_EVENTS.begin(),
                 ic_symm_err.begin(),
                 [](const auto& x, const auto& y) { return std::min(x, y); });
  std::transform(std::execution::par,
                 ic_symm_err.begin(),
                 ic_symm_err.end(),
                 IC_EVENT_ERRPLUS.begin(),
                 ic_symm_err.begin(),
                 [](const auto& x, const auto& y) { return x > 0 ? x : y; });

  VecD chisq_each_data(IC_EVENTS.size(), 0.0);
  VecD event_midpoint(IC_EVENTS.size(), 0.0);
  std::transform(std::execution::par,
                 IC_EVENT_ERRPLUS.cbegin(),
                 IC_EVENT_ERRPLUS.cend(),
                 IC_EVENT_ERRMINUS.cbegin(),
                 event_midpoint.begin(), 
                 [](const auto &x, const auto &y) { return 0.5 * (x - y); });
  std::transform(std::execution::par,
                 IC_EVENTS.cbegin(),
                 IC_EVENTS.cend(),
                 event_midpoint.cbegin(),
                 event_midpoint.begin(), 
                 std::plus<>{});
  std::transform(std::execution::par,
                 event_midpoint.cbegin(),
                 event_midpoint.cend(),
                 MODEL_TOTAL_EVENTS.cbegin(),
                 chisq_each_data.begin(),
                 // loglike);
                 [](const auto& x, const auto& y) { return y - x; });
  std::transform(
    std::execution::par,
    chisq_each_data.cbegin(),
    chisq_each_data.cend(),
    ic_symm_err.cbegin(),
    // IC_EVENTS.cbegin(),
    chisq_each_data.begin(),
    [](const auto& x, const auto& y) { return y > 0 ? (x*x) / (y*y) : 0.0; });

  fmt::println("# {:s} {:>15s} {:>12s}", "gamma", "mKK / MeV", "chi2");
  fmt::print("{:.3f} {:15.3F}", gamma, mKK);
  fmt::println(
    "{:15.5F}",
    std::reduce(
      chisq_each_data.cbegin() + THRESHOLD_BIN, // omit until 60 TeV bin
      chisq_each_data.cend(),
      0.0) /
      (chisq_each_data.size() - 1 - 3));

  const std::string CHISQ_OUTPUT_FILENAME{ std::regex_replace(
    MODEL_EVENT_FILENAME, std::regex("Events_(.*).txt"), "ChiSqNew_$1.tsv") };
  std::FILE* fout = std::fopen(CHISQ_OUTPUT_FILENAME.c_str(), "w+");

  if (fout) {
    fmt::println(fout,
        "#{:s} {:>10s} {:>15s} {:>15s} {:>15s} {:>15s}",
        "IC_EVT", "IC_ERR-", "IC_ERR+", "IC_ESYM", "MODEL_EVT", "CHISQ"
    );
    for (size_t bin_id{}; bin_id < IC_EVENTS.size(); bin_id++) {
      fmt::println(fout,
          "{:2.0F} {:15.3F} {:15.3F} {:15.3F} {:15.3F} {:15.3F}",
          IC_EVENTS.at(bin_id),
          IC_EVENT_ERRMINUS.at(bin_id),
          IC_EVENT_ERRPLUS.at(bin_id),
          ic_symm_err.at(bin_id),
          MODEL_TOTAL_EVENTS.at(bin_id),
          chisq_each_data.at(bin_id));
    }
    std::fclose(fout);
  }

  return 0;
}
