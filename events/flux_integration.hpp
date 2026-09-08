// vim: set ai et ts=2 sw=2 tw=80:

#pragma once
#include "tabdata.hpp"
#include <vector>

using VecDD = std::vector<std::vector<double>>;

VecDD rebin_and_integrate(const TabData &, const VecD &);
