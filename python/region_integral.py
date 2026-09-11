"""
Fast rectangular-region normalization integral for PpqN / PpqG.

Computes integral_region PpqN(Ep, Eq) dEp dEq (and the same for PpqG)
over a rectangle [ep_min, ep_max] x [eq_min, eq_max] -- e.g. to normalize
a likelihood to its fit region.

Why not scipy.integrate.quad
-----------------------------
Nested adaptive quad evaluates the PDF one (Ep, Eq) point at a time, and
even a single-point call into the compiled library pays its ~3-25 ms/point
*scalar* cost (see ppqfort_pdf.py's docstring) -- quad forced to
subdivide at the band's narrow ridge easily needs hundreds of such calls,
i.e. several seconds.

This module instead builds a ridge-aware fixed grid (same physics as
test/python/band_breakpoints.py's ridge/width derivation: a recoil of
energy Er noiselessly produces Eq = Y(Er)*Er, Ep = Er*(1 + Y(Er)*V /
(1000*eps)), so the band's centroid is a curve parameterized by Er, and
the curve's local width in Eq is sigq(Eq_ridge) and sigp(Ep) combined
through the ridge's local slope) and evaluates the *whole* grid in one
batched, thread-parallel call via make_ppqn_pdf/make_ppqg_pdf (~6 us/point
per the README's performance notes) instead of one Python/ctypes round
trip per quadrature node.

The outer (Ep) direction uses a plain uniform grid -- chisquare_harness.py
already treats plain adaptive quad as safe in that direction, only the
inner (Eq) direction needs ridge-tracking -- and the inner direction uses
a uniform grid over the local ridge +/- n_window_widths local band
widths, clipped to the region.  Both are trapezoid-summed; this is a
nested double sum (the inner window differs at every outer Ep point), not
a flattened tensor-product grid.

n_ep, n_eq_window, n_window_widths are required, explicit arguments (no
defaults) -- they directly control the accuracy of a number that feeds a
likelihood normalization, so silently hardcoding them would hide exactly
the choice a careful caller most needs to see and tune.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pq_dist_v9 as ppq
from ppqfort_pdf import make_ppqg_pdf, make_ppqn_pdf


def _ridge_table(band, *, k, Z, V, eps, er_hi, n_tab=4000):
    """Tabulate the noiseless band ridge Er -> (Ep(Er), Eq(Er))."""
    er_tab = np.geomspace(1e-3, er_hi, n_tab)
    if band == "NR":
        y_tab = ppq.Y(er_tab, k=k, Z=Z)
    elif band == "ER":
        y_tab = np.ones_like(er_tab)
    else:
        raise ValueError(f"band must be 'NR' or 'ER', got {band!r}")
    ep_tab = er_tab * (1.0 + y_tab * V / (1000.0 * eps))
    eq_tab = y_tab * er_tab
    return ep_tab, eq_tab


def _region_integral(band, ep_min, ep_max, eq_min, eq_max, *,
                      n_ep, n_eq_window, n_window_widths,
                      k, Z, V, eps, p0, p10, q0, q10, pdf_func):
    if n_ep < 2 or n_eq_window < 2:
        raise ValueError("n_ep and n_eq_window must each be >= 2")

    # Ep(Er) >= Er always (yield is non-negative), so er_hi need only
    # comfortably exceed the largest Ep or Eq of interest.
    er_hi = max(ep_max, eq_max) * 1.5 + 10.0
    ep_tab, eq_tab = _ridge_table(band, k=k, Z=Z, V=V, eps=eps, er_hi=er_hi)

    def eq_ridge(ep):
        return np.interp(ep, ep_tab, eq_tab)

    ep_grid = np.linspace(ep_min, ep_max, n_ep)
    dep = np.maximum(1e-3, 0.01 * ep_grid)
    slope = (eq_ridge(ep_grid + dep) - eq_ridge(ep_grid - dep)) / (2 * dep)
    ridge = eq_ridge(ep_grid)
    sigp_val = ppq.sigp(ep_grid, eps=eps, V=V, p0=p0, p10=p10)
    sigq_val = ppq.sigq(ridge, q0=q0, q10=q10)
    width = np.hypot(sigq_val, slope * sigp_val)

    eq_lo = np.maximum(eq_min, ridge - n_window_widths * width)
    eq_hi = np.minimum(eq_max, ridge + n_window_widths * width)

    # Build the nested grid: n_ep columns, each with its own (possibly
    # degenerate, if the window is clipped away entirely) Eq sub-grid.
    # eq_cols has shape (n_ep, n_eq_window); a degenerate column (eq_hi
    # <= eq_lo) gets all points pinned to eq_lo, contributing 0 to the
    # inner integral via h_eq below rather than needing a special case.
    valid = eq_hi > eq_lo
    h_eq = np.where(valid, (eq_hi - eq_lo) / (n_eq_window - 1), 0.0)
    j = np.arange(n_eq_window)
    eq_cols = eq_lo[:, None] + j[None, :] * h_eq[:, None]
    ep_cols = np.broadcast_to(ep_grid[:, None], eq_cols.shape)

    ep_flat = np.ascontiguousarray(ep_cols.ravel())
    eq_flat = np.ascontiguousarray(eq_cols.ravel())
    f_flat = pdf_func(ep_flat, eq_flat)
    f_cols = f_flat.reshape(eq_cols.shape)

    # Inner trapezoid sum over each column (Eq direction).
    col_w = np.ones(n_eq_window)
    col_w[0] = col_w[-1] = 0.5
    inner = (f_cols * col_w[None, :]).sum(axis=1) * h_eq

    # Outer trapezoid sum over Ep.
    ep_w = np.ones(n_ep)
    ep_w[0] = ep_w[-1] = 0.5
    h_ep = (ep_max - ep_min) / (n_ep - 1)
    return float((inner * ep_w).sum() * h_ep)


def ppqn_region_integral(ep_min, ep_max, eq_min, eq_max, *,
                          n_ep, n_eq_window, n_window_widths,
                          k, Z, F0, eps, V, p0, p10, q0, q10,
                          n_workers=None):
    """Integral of PpqN over [ep_min,ep_max] x [eq_min,eq_max]."""
    pdf_func = make_ppqn_pdf(k=k, Z=Z, F0=F0, eps=eps, V=V,
                              p0=p0, p10=p10, q0=q0, q10=q10,
                              n_workers=n_workers)
    return _region_integral("NR", ep_min, ep_max, eq_min, eq_max,
                             n_ep=n_ep, n_eq_window=n_eq_window,
                             n_window_widths=n_window_widths,
                             k=k, Z=Z, V=V, eps=eps, p0=p0, p10=p10,
                             q0=q0, q10=q10, pdf_func=pdf_func)


def ppqg_region_integral(ep_min, ep_max, eq_min, eq_max, *,
                          n_ep, n_eq_window, n_window_widths,
                          F0, eps, V, p0, p10, q0, q10,
                          n_workers=None):
    """Integral of PpqG over [ep_min,ep_max] x [eq_min,eq_max]."""
    pdf_func = make_ppqg_pdf(F0=F0, eps=eps, V=V,
                              p0=p0, p10=p10, q0=q0, q10=q10,
                              n_workers=n_workers)
    return _region_integral("ER", ep_min, ep_max, eq_min, eq_max,
                             n_ep=n_ep, n_eq_window=n_eq_window,
                             n_window_widths=n_window_widths,
                             k=None, Z=None, V=V, eps=eps, p0=p0, p10=p10,
                             q0=q0, q10=q10, pdf_func=pdf_func)
