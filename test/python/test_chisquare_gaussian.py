"""
Chi-square distribution test for a bivariate Gaussian PDF.

Generates 1e5 independent datasets from a known PDF, computes the
Pearson chi-square for each using equibin-defined bins and PDF-integrated
expected counts, then checks whether the resulting chi-square values
follow the theoretical chi2(dof) distribution.

A well-behaved sampler + PDF pair should produce:
  - mean chi2  ≈ dof
  - var  chi2  ≈ 2 * dof
  - KS p-value >> 0.05

Runtime note
------------
With default settings (n_events=1000, n_throws=100_000, n_bins=64):
  - Sampling:  ~100s   (1e8 total events at ~1M events/s)
  - Binning:   ~15s
  Total:       ~2 min

Increase n_events for a more sensitive test; increase n_bins to probe
finer structure (keep n_events / n_bins >= 10 for valid chi-square).

Run with:
  python test/python/test_chisquare_gaussian.py
"""

import sys
from pathlib import Path

import numpy as np
from scipy.stats import chi2 as chi2_dist, kstest
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
import sample_from_pdf as spdf
from chisquare_harness import run_chisquare_test

# ---------------------------------------------------------------------------
# PDF definition
# ---------------------------------------------------------------------------
mu_ep,  mu_eq  = 10.0, 5.0
sig_ep, sig_eq = 4.0,  1.0
rho = 0.6

def gaussian_pdf(ep_flat, eq_flat):
    ep = np.asarray(ep_flat, dtype=np.float64)
    eq = np.asarray(eq_flat, dtype=np.float64)
    z = (((ep - mu_ep) / sig_ep) ** 2
         - 2 * rho * (ep - mu_ep) / sig_ep * (eq - mu_eq) / sig_eq
         + ((eq - mu_eq) / sig_eq) ** 2)
    return np.exp(-z / (2 * (1 - rho**2)))

ep_range = (-5.0, 25.0)
eq_range = (1.0,   9.0)

# ---------------------------------------------------------------------------
# Test parameters
# ---------------------------------------------------------------------------
N_EVENTS  = 1_000      # events per chi-square throw
N_THROWS  = 100_000    # chi-square values to accumulate
N_BINS    = 64         # equibin target bins  (E_i ≈ N_EVENTS/N_BINS ≈ 15.6)
GRID_RES  = 200        # lintsampler grid resolution per axis

# ---------------------------------------------------------------------------
# Build sampler and run test
# ---------------------------------------------------------------------------
print("Building sampler...")
sampler = spdf.build_sampler(gaussian_pdf, ep_range, eq_range, GRID_RES, GRID_RES)

chi2_vals, dof, bins, expected = run_chisquare_test(
    gaussian_pdf, sampler, ep_range, eq_range,
    n_events   = N_EVENTS,
    n_throws   = N_THROWS,
    n_bins     = N_BINS,
    batch_size = 5_000,
    seed       = 42,
)

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
print(f"\n{'='*55}")
print(f"Chi-square test: bivariate Gaussian  (dof = {dof})")
print(f"{'='*55}")
print(f"  n_events per throw : {N_EVENTS:,}")
print(f"  n_throws           : {N_THROWS:,}")
print(f"  n_bins             : {len(bins)}  (target {N_BINS})")
print(f"  min expected count : {expected.min():.1f}")
print()
print(f"  Mean chi2   : {chi2_vals.mean():.3f}   (expect {dof})")
print(f"  Var  chi2   : {chi2_vals.var():.3f}   (expect {2*dof})")

ks_stat, ks_p = kstest(chi2_vals, chi2_dist(dof).cdf)
print(f"  KS statistic: {ks_stat:.4f}")
print(f"  KS p-value  : {ks_p:.4f}   (expect >> 0.05)")
print()
passed = ks_p > 0.05 and abs(chi2_vals.mean() - dof) / dof < 0.02
if passed:
    print("PASS — chi-square distribution consistent with chi2(dof).")
else:
    print("FAIL — chi-square distribution deviates from chi2(dof).")

# ---------------------------------------------------------------------------
# Reference plot
# ---------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))

x = np.linspace(chi2_dist(dof).ppf(1e-4), chi2_dist(dof).ppf(1 - 1e-4), 500)
ax.plot(x, chi2_dist(dof).pdf(x), color="C1", lw=2, label=rf"$\chi^2({dof})$ theory")
ax.hist(chi2_vals, bins=80, density=True, alpha=0.6, color="C0", label=f"{N_THROWS:,} throws")

ax.set_xlabel(r"$\chi^2$")
ax.set_ylabel("Probability density")
ax.set_title("Chi-square distribution: bivariate Gaussian")

stats_text = (
    f"mean = {chi2_vals.mean():.2f}  (expect {dof})\n"
    f"var  = {chi2_vals.var():.1f}  (expect {2*dof})\n"
    f"KS p = {ks_p:.3f}"
)
ax.text(0.97, 0.97, stats_text, transform=ax.transAxes,
        ha="right", va="top", fontsize=9,
        bbox=dict(boxstyle="round,pad=0.3", fc="white", alpha=0.8))

ax.legend()
fig.tight_layout()
outfile = REPO_ROOT / "figures" / "chisquare_gaussian.png"
fig.savefig(outfile, dpi=150)
print(f"\nPlot saved to {outfile}")
