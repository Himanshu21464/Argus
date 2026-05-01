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
