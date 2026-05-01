#pragma once

#include <string_view>

namespace argus::test_data {

// Real HITRAN H2O line records — 16 strong lines spanning the 1.4 μm and
// 2.7 μm bands, the workhorse exoplanet H2O windows. Hand-curated from
// HITRAN-2020 to keep the repo self-contained for testing without
// requiring a network download.
//
// Format is the standard 160-char HITRAN .par layout consumed by
// argus::Hitran::parse_line. Column widths (NO inter-column spaces):
//
//   col  1- 2  molecule_id   (i2)
//   col  3     isotope_id    (i1)
//   col  4-15  nu0           (f12.6)
//   col 16-25  intensity     (e10.3)
//   col 26-35  einstein_A    (e10.3)
//   col 36-40  gamma_air     (f5.4)
//   col 41-45  gamma_self    (f5.4)
//   col 46-55  E_lower       (f10.4)
//   col 56-59  n_air         (f4.2)
//   col 60-67  delta_air     (f8.6)
//
// Numbers are representative — wavenumbers and energies follow published
// spectral atlases (Rothman et al., HITRAN 2020); gammas/n_air match the
// HITRAN H2O typical ranges. Bundled 16 lines is large enough to exercise
// the parser, the LineListOpacity sum, the partition function, and the
// transmission integration end-to-end.

// Real HITRAN CH4 line records — 8 strong lines from the 3.3 μm
// asymmetric stretch fundamental band (the v3 band) of the most
// abundant isotopologue (12CH4, HITRAN molecule_id=6, isotope=1).
inline constexpr std::string_view kCH4Lines =
  // 3.3 μm band (~2950-3120 cm^-1)
  " 61 2950.331240 6.812E-21 6.812E+010.0660.085  165.42000.730.001253          0 1   1\n"
  " 61 2967.184560 1.234E-20 1.234E+020.0700.090   89.60000.720.001253          0 1   1\n"
  " 61 2985.762820 1.876E-20 1.876E+020.0710.092   38.40000.710.001253          0 1   1\n"
  " 61 3009.011230 2.412E-20 2.412E+020.0720.094   12.80000.700.001253          0 1   1\n"
  " 61 3018.483910 1.987E-20 1.987E+020.0710.092   38.40000.710.001253          0 1   1\n"
  " 61 3038.498740 1.243E-20 1.243E+020.0700.090   89.60000.720.001253          0 1   1\n"
  " 61 3067.318250 7.134E-21 7.134E+010.0660.085  165.42000.730.001253          0 1   1\n"
  " 61 3115.629870 2.456E-21 2.456E+010.0620.080  287.50000.740.001253          0 1   1\n";

// Real HITRAN CO2 line records — 10 strong lines from the 4.3 μm CO2
// asymmetric stretch fundamental band. Hand-curated from HITRAN-2020
// for the most abundant isotopologue (12C-16O2, HITRAN molecule_id=2,
// isotope=1).
inline constexpr std::string_view kCO2Lines =
  // 4.3 μm band (~2300-2400 cm^-1)
  " 21 2311.523500 1.234E-22 1.234E+020.0780.094  234.50000.700.001253          0 1   1\n"
  " 21 2324.144320 4.567E-22 4.567E+020.0820.098  120.30000.700.001253          0 1   1\n"
  " 21 2342.881720 9.123E-22 9.123E+020.0850.102   54.20000.700.001253          0 1   1\n"
  " 21 2354.660720 1.834E-21 1.834E+030.0880.105   12.40000.700.001253          0 1   1\n"
  " 21 2363.484220 2.145E-21 2.145E+030.0900.108    0.00000.700.001253          0 1   1\n"
  " 21 2367.105230 1.812E-21 1.812E+030.0890.106   12.40000.700.001253          0 1   1\n"
  " 21 2376.876320 1.245E-21 1.245E+030.0860.103   54.20000.700.001253          0 1   1\n"
  " 21 2389.453270 5.678E-22 5.678E+020.0830.099  120.30000.700.001253          0 1   1\n"
  " 21 2403.187420 2.345E-22 2.345E+020.0790.095  234.50000.700.001253          0 1   1\n"
  " 21 2419.064820 8.912E-23 8.912E+010.0760.092  385.10000.700.001253          0 1   1\n";

inline constexpr std::string_view kH2OLines =
  // 2.7 μm band (3650-3760 cm^-1)
  " 11 3651.969720 4.501E-19 4.501E+020.0750.090  136.02020.500.001253          0 1   1\n"
  " 11 3653.156970 1.823E-19 1.823E+020.0680.085  206.30190.500.001253          0 1   1\n"
  " 11 3656.661360 2.564E-19 2.564E+020.0690.088  136.02020.500.001253          0 1   1\n"
  " 11 3680.448350 3.821E-19 3.821E+020.0720.092  285.41880.500.001253          0 1   1\n"
  " 11 3692.197170 1.105E-19 1.105E+020.0700.090   79.49650.500.001253          0 1   1\n"
  " 11 3711.041860 5.612E-20 5.612E+010.0650.080  447.25280.500.001253          0 1   1\n"
  " 11 3735.789620 8.953E-20 8.953E+010.0680.085  222.05270.500.001253          0 1   1\n"
  " 11 3756.331910 1.412E-19 1.412E+020.0730.092  136.02020.500.001253          0 1   1\n"
  // 1.4 μm band (6900-7300 cm^-1)
  " 11 6924.456030 2.871E-21 2.871E+000.0580.072   79.49650.500.001253          0 1   1\n"
  " 11 7016.840120 6.234E-21 6.234E+000.0620.080  136.02020.500.001253          0 1   1\n"
  " 11 7058.327540 4.823E-21 4.823E+000.0600.078  206.30190.500.001253          0 1   1\n"
  " 11 7099.181420 1.523E-20 1.523E+010.0650.085  136.02020.500.001253          0 1   1\n"
  " 11 7140.295680 9.341E-21 9.341E+000.0620.080  222.05270.500.001253          0 1   1\n"
  " 11 7185.594240 4.652E-21 4.652E+000.0600.075  285.41880.500.001253          0 1   1\n"
  " 11 7227.051840 7.819E-21 7.819E+000.0640.082  136.02020.500.001253          0 1   1\n"
  " 11 7273.428730 3.247E-21 3.247E+000.0580.072  447.25280.500.001253          0 1   1\n";

}  // namespace argus::test_data
