#pragma once

#include <numbers>
#include <array>

typedef std::array<double, 3> arr3;

using std::numbers::pi;

const double MeV_TO_eV{1.0E06};
const double GeV_TO_MeV{1.0E3};

template <typename T> inline T sqr(const T t) { return t * t; }

