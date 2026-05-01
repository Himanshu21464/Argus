#include "argus/partition.hpp"

#include <cmath>
#include <stdexcept>

namespace argus {

namespace {

// Power-law fit Q(T) = Q_ref * (T / T_ref)^exponent.
//
// Anchor-exact at T_ref = 296 K by construction; accurate to ~10% over
// 200-3000 K for non-vibrating molecules. M3 will replace these with the
// full TIPS-2017 polynomial (a0 + a1*T + a2*T² + a3*T³ + a4*T⁴) loaded
// from a tabulated grid.
//
// Exponents derived from anchor pairs (Q at 1000 K vs Q at 296 K) using
// published TIPS-2017 values:
//
//   H2O(161): Q(1000)/Q(296) = 1296/174.58 = 7.42  -> n = log(7.42)/log(3.378) = 1.648
//   CO2(626): Q(1000)/Q(296) = 1019/286.94 = 3.55  -> n = log(3.55)/log(3.378) = 1.039
//   CH4(211): Q(1000)/Q(296) = 4014/590.48 = 6.80  -> n = log(6.80)/log(3.378) = 1.575
//   CO(26):   Q(1000)/Q(296) = 365/107.42  = 3.40  -> n = log(3.40)/log(3.378) = 1.005
//   NH3:      Q(1000)/Q(296) = 9700/1725  = 5.62   -> n = log(5.62)/log(3.378) = 1.418

struct Fit {
  int    molecule_id;
  int    isotope_id;
  std::string key;
  double q_ref;       // Q(T_ref = 296 K)
  double exponent;    // Q(T) = q_ref * (T/T_ref)^exponent
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
