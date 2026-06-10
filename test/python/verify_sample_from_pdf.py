"""
Verify sample_from_pdf with a skewed 2D Gaussian whose moments are known exactly.

The test PDF is a bivariate Gaussian with
  mean        = (mu_ep, mu_eq) = (10, 5)
  std         = (sig_ep, sig_eq) = (4, 1)
  correlation = rho = 0.6

All first- and second-order moments are checked against their analytic values.

Run with:
  python test/python/verify_sample_from_pdf.py
"""

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "python"))
import sample_from_pdf as spdf

# --- analytic parameters ---
mu_ep,  mu_eq  = 10.0, 5.0
sig_ep, sig_eq = 4.0,  1.0
rho = 0.6

def skewed_gauss(ep_flat, eq_flat):
    z = (((ep_flat - mu_ep) / sig_ep) ** 2
         - 2 * rho * (ep_flat - mu_ep) / sig_ep * (eq_flat - mu_eq) / sig_eq
         + ((eq_flat - mu_eq) / sig_eq) ** 2)
    return np.exp(-z / (2 * (1 - rho**2)))

# --- build the sampler ---
print("Building sampler (PDF evaluated on uniform grid)...")
sample = spdf.build_sampler(
    skewed_gauss,
    ep_range=(-5, 25),
    eq_range=(1, 9),
    n_ep=200,
    n_eq=200,
)

# --- draw samples ---
N = 500_000
print(f"Drawing {N:,} samples...")
Ep, Eq = sample(N, seed=42)

# --- moment checks ---
checks = [
    ("E[Ep]",        Ep.mean(),                   mu_ep),
    ("E[Eq]",        Eq.mean(),                   mu_eq),
    ("std(Ep)",      Ep.std(),                    sig_ep),
    ("std(Eq)",      Eq.std(),                    sig_eq),
    ("corr(Ep,Eq)",  np.corrcoef(Ep, Eq)[0, 1],  rho),
]

print(f"\n{'Quantity':<18} {'Sample':>10} {'Analytic':>10} {'Diff':>10}")
print("-" * 52)
all_pass = True
for name, sample_val, analytic_val in checks:
    diff = sample_val - analytic_val
    flag = "" if abs(diff) < 0.05 else "  <-- FAIL"
    if flag:
        all_pass = False
    print(f"{name:<18} {sample_val:>10.4f} {analytic_val:>10.4f} {diff:>+10.4f}{flag}")

print()
print("PASS" if all_pass else "FAIL")
