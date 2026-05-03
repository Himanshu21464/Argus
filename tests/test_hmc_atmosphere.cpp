// M3 substrate proof for the differentiable-physics claim:
// Hamiltonian Monte Carlo end-to-end through a real atmospheric
// forward model — a 3-H₂O-line isothermal layer with temperature-
// dependent Doppler + Lorentz broadening — to recover (T, log₁₀ VMR)
// from a noisy synthetic spectrum.
//
// HMC needs gradients, which Argus computes via forward-mode autograd
// (Dual<double>). The forward is templated on the scalar type T so
// HMC's `grad` machinery can call it with Dual-seeded inputs and read
// out ∂log p/∂(T_K, log₁₀ VMR) without finite differencing.
//
// This is the test that validates "the substrate is differentiable
// across the physics layer, not just the inference layer". Earlier
// HMC tests use analytic Gaussians; this one runs through Voigt line
// shapes + temperature-scaled widths + chi² likelihood.

#include <cassert>
#include <cmath>
#include <vector>

#include "argus/argus.hpp"

namespace {

// One H₂O-like line. Reference parameters at T_ref = 1500 K.
struct LineRef {
  double nu0;            // line centre, cm⁻¹
  double S0_ref;         // intensity at T_ref (in arbitrary cross-section units)
  double sigma_ref;      // Doppler width at T_ref, cm⁻¹
  double gamma_ref;      // Lorentz HWHM at T_ref, cm⁻¹
};

// Templated single-layer transmission depth (1 − exp(−τ)).
//   τ = N_col · VMR · Σ_lines S(T) · V(ν − ν₀; σ(T), γ(T))
//   σ(T) = σ_ref · √(T/T_ref)
//   γ(T) = γ_ref · √(T_ref/T)         (rough P/T scaling)
//   S(T) = S_ref                         (intensity T-correction omitted —
//                                         we want a clean substrate test)
template <typename T>
T transit_depth(T temperature, T log10_vmr, double nu_obs,
                const std::vector<LineRef>& lines) {
  using std::exp; using std::sqrt;
  using argus::exp; using argus::sqrt;
  constexpr double T_REF = 1500.0;
  constexpr double N_COL = 1.0e22;            // column density (arbitrary)
  constexpr double LN10  = 2.302585092994046;

  T tau = T{0.0};
  for (const auto& L : lines) {
    const T x       = T{nu_obs - L.nu0};
    const T sigma_g = T{L.sigma_ref} * sqrt(temperature / T{T_REF});
    const T gamma_l = T{L.gamma_ref} * sqrt(T{T_REF} / temperature);
    const T xsec    = T{L.S0_ref} * argus::voigt(x, sigma_g, gamma_l);
    tau = tau + xsec;
  }
  // 10^x = exp(x · ln 10) — Dual<T> has exp but only pow(Dual, scalar).
  const T vmr = exp(log10_vmr * T{LN10});
  tau = tau * T{N_COL} * vmr;
  return T{1.0} - exp(-tau);
}

}  // namespace

int main() {
  using namespace argus;
  using D = Dual<double>;

  // ─── 1. Truth: 3-line H₂O-like band, isothermal layer at 1500 K. ──
  std::vector<LineRef> lines{
    {1500.0, 1.0e-21, 0.04, 0.05},
    {1502.5, 0.7e-21, 0.04, 0.06},
    {1505.0, 1.2e-21, 0.04, 0.05},
  };
  std::vector<double> wn;
  for (int i = 0; i < 21; ++i) wn.push_back(1499.0 + i * 0.35);

  const double TRUE_T   = 1500.0;
  const double TRUE_LV  = -3.0;

  std::vector<double> truth_spec(wn.size());
  for (std::size_t i = 0; i < wn.size(); ++i) {
    truth_spec[i] = transit_depth<double>(TRUE_T, TRUE_LV, wn[i], lines);
  }

  // Add Gaussian noise (synthetic JWST-PRISM scale).
  std::mt19937_64 rng(2026);
  const double NOISE_SIGMA = 5.0e-4;
  std::normal_distribution<double> nz(0.0, NOISE_SIGMA);
  std::vector<double> observed(wn.size());
  for (std::size_t i = 0; i < wn.size(); ++i) {
    observed[i] = truth_spec[i] + nz(rng);
  }

  // ─── 2. Templated log-posterior, in STANDARDISED parameters. ──────
  // T_K       = 1500 + 200 · s[0]
  // log10_VMR =   -3 +   2 · s[1]
  // Both s[i] are O(1) so a single HMC step size works. Soft N(0,1)
  // priors on (s[0], s[1]) are broad enough to be non-restrictive at
  // truth.
  auto logp = [&](const std::vector<D>& s) -> D {
    const D temperature = D{1500.0, 0.0} + D{200.0, 0.0} * s[0];
    const D lvmr        = D{-3.0,   0.0} + D{2.0,   0.0} * s[1];

    D chi2{0.0, 0.0};
    for (std::size_t i = 0; i < wn.size(); ++i) {
      D model = transit_depth<D>(temperature, lvmr, wn[i], lines);
      D r = (D{observed[i], 0.0} - model) / D{NOISE_SIGMA, 0.0};
      chi2 = chi2 + r * r;
    }
    D log_prior = D{0.0, 0.0} - D{0.5, 0.0} * (s[0] * s[0] + s[1] * s[1]);
    D log_lik   = D{0.0, 0.0} - D{0.5, 0.0} * chi2;
    return log_prior + log_lik;
  };

  // Sanity: log p at truth is finite (truth is s = (0, 0)).
  std::vector<D> truth_d{D{0.0, 0.0}, D{0.0, 0.0}};
  assert(std::isfinite(logp(truth_d).v));

  // ─── 3. Run HMC starting from a deliberately wrong guess. ─────────
  // Initial guess: s = (-0.5, +0.1) → T = 1400 K, log10_VMR = -2.8.
  // Step size and leapfrog count tuned for the actual posterior
  // curvature: the chi² Hessian on the standardised parameters has
  // condition number ~10², so step ≈ 0.001 keeps the leapfrog
  // integrator stable end-to-end.
  HMC sampler(logp, /*step_size=*/0.001, /*n_leapfrog=*/200, /*seed=*/2026);
  std::vector<double> state{-0.5, 0.1};
  auto result = sampler.sample(state, 1500);

  const double acc = sampler.acceptance_rate();
  assert(acc > 0.50);

  // ─── 4. Recovery: unstandardise samples and check vs truth. ───────
  std::vector<double> Ts, LVs;
  for (std::size_t i = result.samples.size() / 2; i < result.samples.size();
       ++i) {
    Ts.push_back(1500.0 + 200.0 * result.samples[i][0]);
    LVs.push_back(-3.0  +   2.0 * result.samples[i][1]);
  }
  auto mean_std = [](const std::vector<double>& v) {
    double m = 0.0;
    for (double x : v) m += x;
    m /= static_cast<double>(v.size());
    double s = 0.0;
    for (double x : v) { const double d = x - m; s += d * d; }
    s = std::sqrt(s / static_cast<double>(v.size()));
    return std::pair<double, double>{m, s};
  };
  auto [T_mean, T_std]   = mean_std(Ts);
  auto [LV_mean, LV_std] = mean_std(LVs);

  assert(std::fabs(T_mean - TRUE_T) < 3.0 * T_std);
  assert(std::fabs(LV_mean - TRUE_LV) < 3.0 * LV_std);

  // ─── 5. Determinism: same seed → bit-equal samples. ───────────────
  HMC s2(logp, 0.001, 200, 2026);
  std::vector<double> st2{-0.5, 0.1};
  auto r2 = s2.sample(st2, 1500);
  assert(result.samples.size() == r2.samples.size());
  for (std::size_t i = 0; i < result.samples.size(); ++i) {
    for (std::size_t d = 0; d < 2; ++d) {
      assert(result.samples[i][d] == r2.samples[i][d]);
    }
  }

  return 0;
}
