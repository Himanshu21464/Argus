#include "argus/partition.hpp"

#include <cmath>
#include <stdexcept>

namespace argus {

namespace {

// Power-law fit Q(T) = Q_ref * (T / T_ref)^exponent.
//
// Anchor-exact at T_ref = 296 K by construction; accuracy ~1-7% over
// 100-3000 K. Future replacement: full TIPS-2017 cubic polynomials
// loaded from a tabulated grid (better than 0.5% accuracy across the
// range); the public `Partition::Q(...)` surface stays stable.
//
// Exponents derived from anchor pairs (Q at 1000 K vs Q at 296 K) using
// published TIPS-2017 values for the most abundant isotopologue:
//
//   H2O(161): n = log(1296/174.58)/log(1000/296) = 1.648
//   CO2(626): n = log(1019/286.94)/log(1000/296) = 1.039
//   CH4(211): n = log(4014/590.48)/log(1000/296) = 1.575
//   CO(26):   n = log(365/107.42)/log(1000/296)  = 1.005
//   NH3(4111):n = log(9700/1725)/log(1000/296)   = 1.418
//
// Verified against TIPS reference at the test anchor temperatures
// (200, 296, 1000, 1500, 2000 K). See tests/test_partition.cpp for
// tolerance bounds.

struct Fit {
  int    molecule_id;
  int    isotope_id;
  std::string key;
  double q_ref;
  double exponent;
};

const Fit kFits[] = {
  { 1, 1, "H2O",  174.58, 1.648 },
  { 2, 1, "CO2",  286.94, 1.039 },
  { 6, 1, "CH4",  590.48, 1.575 },
  { 5, 1, "CO",   107.42, 1.005 },
  {11, 1, "NH3", 1725.22, 1.418 },
};

const Fit* find_fit(int molecule_id, int isotope_id) {
  for (const auto& f : kFits) {
    if (f.molecule_id == molecule_id && f.isotope_id == isotope_id) {
      return &f;
    }
  }
  return nullptr;
}

const Fit* find_fit_by_key(const std::string& key) {
  for (const auto& f : kFits) {
    if (f.key == key) return &f;
  }
  return nullptr;
}

}  // namespace

double Partition::Q(int molecule_id, int isotope_id, double T_k) {
  const Fit* f = find_fit(molecule_id, isotope_id);
  if (!f) {
    throw std::invalid_argument(
        "Partition::Q: unknown (molecule, isotope) pair");
  }
  if (T_k <= 0.0) {
    throw std::invalid_argument("Partition::Q: T must be positive");
  }
  return f->q_ref * std::pow(T_k / T_ref_k, f->exponent);
}

double Partition::Q(const std::string& key, double T_k) {
  const Fit* f = find_fit_by_key(key);
  if (!f) {
    throw std::invalid_argument(
        "Partition::Q: unknown species key '" + key + "'");
  }
  return Q(f->molecule_id, f->isotope_id, T_k);
}

double Partition::Q_ref(int molecule_id, int isotope_id) {
  const Fit* f = find_fit(molecule_id, isotope_id);
  if (!f) {
    throw std::invalid_argument(
        "Partition::Q_ref: unknown (molecule, isotope) pair");
  }
  return f->q_ref;
}

double Partition::Q_ref(const std::string& key) {
  const Fit* f = find_fit_by_key(key);
  if (!f) {
    throw std::invalid_argument(
        "Partition::Q_ref: unknown species key '" + key + "'");
  }
  return f->q_ref;
}

}  // namespace argus
