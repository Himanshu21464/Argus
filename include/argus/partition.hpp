#pragma once

#include <string>

namespace argus {

// Total Internal Partition Sum Q(T) for a HITRAN molecule isotopologue.
//
// Used to scale line intensities away from HITRAN's reference temperature
// of 296 K:
//
//     S(T) = S(296) * (Q(296) / Q(T))
//          * exp(-c2 * E_lower * (1/T - 1/Tref))
//          * (1 - exp(-c2 * nu / T)) / (1 - exp(-c2 * nu / Tref))
//
// The third factor (induced-emission correction) is small for nu >> kT/h
// and is collapsed into the simple Boltzmann form for M2.
//
// M2 ships analytic polynomial fits for the three workhorses (H2O, CO2,
// CH4) of exoplanet atmospheres. M3 will load the full TIPS-2017 tables.
class Partition {
 public:
  // Return Q(T) for (molecule, isotope). Throws std::invalid_argument if
  // the (molecule, isotope) pair is unknown.
  // molecule_id follows HITRAN convention: 1=H2O, 2=CO2, 6=CH4 ...
  // isotope_id : 1 for the most abundant isotopologue
  static double Q(int molecule_id, int isotope_id, double T_k);

  // Convenience: take the canonical species key ("H2O", "CO2", "CH4", ...)
  // and assume isotope = 1.
  static double Q(const std::string& species_key, double T_k);

  // Q(296 K) — the HITRAN reference value.
  static double Q_ref(int molecule_id, int isotope_id);
  static double Q_ref(const std::string& species_key);

  // HITRAN reference temperature.
  static constexpr double T_ref_k = 296.0;
};

}  // namespace argus
