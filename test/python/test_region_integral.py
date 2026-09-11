"""
Validates the fast rectangular-region normalization integral three ways:

1. Fortran-native (ppqfort_pdf.ppqn_region/ppqg_region) vs. the Python
   batched-grid reimplementation (region_integral.ppqn_region_integral/
   ppqg_region_integral) -- same algorithm, two languages, should agree
   to near machine precision.
2. Both vs. the existing, algorithmically independent nested-quad
   reference (chisquare_harness.expected_counts_from_pdf, treating the
   whole region as one "bin", via band_breakpoints.make_ridge_breakpoints
   for the required ridge breakpoints) -- ground truth.
3. A timing comparison of all three on the same region.

Run from the repository root with:
  LD_LIBRARY_PATH=lib python test/python/test_region_integral.py
"""

import os
import sys
import time

import numpy as np

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(REPO_ROOT, "python"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from ppqfort_pdf import make_ppqg_pdf, make_ppqn_pdf, ppqg_region, ppqn_region
from region_integral import ppqg_region_integral, ppqn_region_integral
from band_breakpoints import make_ridge_breakpoints
from chisquare_harness import debug_bin

PARAMS_NR = dict(k=0.18, Z=32.0, F0=0.122, eps=3.0e-3, V=3.0,
                  p0=0.06421907, p10=0.48998486, q0=0.23718488, q10=0.27093151)
PARAMS_ER = dict(F0=0.122, eps=3.0e-3, V=3.0,
                  p0=0.06421907, p10=0.48998486, q0=0.23718488, q10=0.27093151)

# n_ep=300 converges to ~0.1% of the quad reference for a 247 keV-wide
# region (the widest tested here) -- trapezoid error drops ~4x per
# doubling of n_ep, so 1200 gets comfortably under 1e-4 for that region
# while staying fast; narrower regions converge faster still.  Wider
# regions than these need a correspondingly larger n_ep.
GRID = dict(n_ep=1200, n_eq_window=81, n_window_widths=10.0)

# (ep_min, ep_max, eq_min, eq_max, label) -- a few representative regions:
# a typical fit window, a narrow low-energy slice, and one intentionally
# off-band to confirm it correctly integrates to ~0.
REGIONS = [
    (3.0, 250.0, 0.5, 130.0, "typical fit region"),
    (3.0, 20.0, 0.5, 10.0, "low-energy slice"),
    (200.0, 250.0, 0.5, 5.0, "off-band (high Ep, tiny Eq)"),
]


def quad_reference(band, ep_min, ep_max, eq_min, eq_max, params):
    """The user's current method: nested scipy.integrate.quad with ridge
    breakpoints (chisquare_harness.debug_bin's "fixed" value -- the same
    computation expected_counts_from_pdf does, but unnormalized so it's
    the raw integral).  Also cross-checks against debug_bin's independent
    brute-force Simpson grid_ref."""
    if band == "NR":
        pdf_func = make_ppqn_pdf(**params, n_workers=os.cpu_count())
        ridge = make_ridge_breakpoints("NR", **{k: params[k] for k in
                                                 ("k", "Z", "eps", "V", "p0", "p10", "q0", "q10")},
                                        er_max=700.0, n_window_widths=10.0)
    else:
        pdf_func = make_ppqg_pdf(**params, n_workers=os.cpu_count())
        ridge = make_ridge_breakpoints("ER", k=None, Z=None,
                                        **{k: params[k] for k in
                                           ("eps", "V", "p0", "p10", "q0", "q10")},
                                        er_max=700.0, n_window_widths=10.0)
    import io
    from contextlib import redirect_stdout
    with redirect_stdout(io.StringIO()):
        result = debug_bin(pdf_func, (ep_min, ep_max, eq_min, eq_max), ridge)
    rel_fixed_grid = abs(result["fixed"] - result["grid_ref"]) / max(abs(result["grid_ref"]), 1e-300)
    if rel_fixed_grid > 1e-2:
        print(f"    NOTE: debug_bin's own fixed vs grid_ref disagree by {rel_fixed_grid:.3e}")
    return float(result["fixed"])


def run(band, params, region_integral_func, region_func):
    print(f"\n{'='*70}")
    print(f"Band: {band}")
    print(f"{'='*70}")
    all_ok = True
    for ep_min, ep_max, eq_min, eq_max, label in REGIONS:
        t0 = time.time()
        ref = quad_reference(band, ep_min, ep_max, eq_min, eq_max, params)
        t_quad = time.time() - t0

        t0 = time.time()
        py_val = region_integral_func(ep_min, ep_max, eq_min, eq_max, **GRID, **params)
        t_py = time.time() - t0

        t0 = time.time()
        fort_val = region_func(ep_min, ep_max, eq_min, eq_max, **GRID, **params)
        t_fort = time.time() - t0

        rel_py_fort = abs(py_val - fort_val) / max(abs(fort_val), 1e-300)
        rel_py_quad = abs(py_val - ref) / max(abs(ref), 1e-300)
        rel_fort_quad = abs(fort_val - ref) / max(abs(ref), 1e-300)

        print(f"\n  {label}: Ep in [{ep_min},{ep_max}], Eq in [{eq_min},{eq_max}]")
        print(f"    quad reference : {ref:.10e}   ({t_quad*1000:.1f} ms)")
        print(f"    python grid    : {py_val:.10e}   ({t_py*1000:.1f} ms)  rel vs quad: {rel_py_quad:.3e}")
        print(f"    fortran region : {fort_val:.10e}   ({t_fort*1000:.1f} ms)  rel vs quad: {rel_fort_quad:.3e}")
        print(f"    python vs fortran rel diff: {rel_py_fort:.3e}")

        ok = rel_py_fort < 1e-6 and rel_py_quad < 5e-4 and rel_fort_quad < 5e-4
        print(f"    {'PASS' if ok else 'FAIL'}")
        all_ok = all_ok and ok
    return all_ok


if __name__ == "__main__":
    ok_nr = run("NR", PARAMS_NR, ppqn_region_integral, ppqn_region)
    ok_er = run("ER", PARAMS_ER, ppqg_region_integral, ppqg_region)

    print(f"\n{'='*70}")
    if ok_nr and ok_er:
        print("ALL PASS")
    else:
        print("FAILURES ABOVE")
        sys.exit(1)
