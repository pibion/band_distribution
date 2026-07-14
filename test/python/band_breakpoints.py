"""
Ridge breakpoints for integrating the band PDFs over histogram bins.

The problem
-----------
The expected count in a bin is the integral of the PDF over the bin
rectangle, computed in chisquare_harness.expected_counts_from_pdf with
nested adaptive quadrature (scipy quad over Ep, inner quad over Eq).
Adaptive quadrature starts from a fixed set of sample points across the
integration interval: if the PDF is a narrow ridge inside a wide bin,
every one of those initial points can land in the ~zero region on
either side of the ridge, and quad confidently returns ~0 without ever
subdividing toward the peak.  See quad_breakpoints_demo.py in this
directory for a self-contained demonstration of the failure and the fix.

The fix is quad's `points=` argument: telling the inner (Eq) quadrature
where the ridge is forces a subdivision there, so the peak cannot be
missed.  This module computes those breakpoints from the band physics.

The physics
-----------
A recoil of true energy Er produces (noiselessly)

    N  = Y(Er) * Er / eps          e/h pairs, so
    Eq = Y(Er) * Er                (ionization channel)
    Ep = Er + (V/1000) * N
       = Er * (1 + Y(Er) * V / (1000 * eps))   (phonon + Neganov-Luke)

with Y(Er) = a * Er**b for nuclear recoils and Y = 1 for electron
recoils.  Eliminating Er gives the band centroid: the most probable Eq
at a given measured Ep.  Both Ep(Er) and Eq(Er) are monotone, so the
centroid is tabulated on an Er grid once and inverted by interpolation.

The band's local width in Eq at fixed Ep combines the charge resolution
with the phonon resolution mapped through the ridge slope:

    width(Ep) = sqrt( sigq(Eq_ridge)^2 + (dEq_ridge/dEp * sigp(Ep))^2 )

The breakpoints handed to the harness are the ridge center and a window
of +/- 10 local widths around it; the harness clips them to each bin.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "python"))
import pq_dist_v8 as ppq


def make_ridge_breakpoints(band, *, a=0.16, b=0.18, eps=3.0e-3, V=3.0,
                           p0=0.06421907, p10=0.48998486,
                           q0=0.23718488, q10=0.27093151,
                           er_max=700.0, n_window_widths=10.0):
    """
    Build an inner_points_func for chisquare_harness.expected_counts_from_pdf.

    Parameters
    ----------
    band : "NR" or "ER"
        Nuclear recoils use the yield model Y(Er) = a * |Er|^b; electron
        recoils use Y = 1 (the a and b arguments are then ignored).
    a, b, eps, V, p0, p10, q0, q10 : float
        Detector parameters, matching the PDF being integrated.
    er_max : float
        Upper end of the tabulated Er range (keV); must comfortably
        exceed the largest Ep of interest.
    n_window_widths : float
        Half-width of the breakpoint window in units of the local band
        width.

    Returns
    -------
    inner_points_func : callable
        inner_points_func(ep, ymin, ymax) -> [lo, ridge, hi], the Eq
        breakpoints for the inner quadrature at phonon energy ep.  The
        harness drops any that fall outside the bin's (ymin, ymax).
    """
    if band == "NR":
        y_a, y_b = a, b
    elif band == "ER":
        y_a, y_b = 1.0, 0.0
    else:
        raise ValueError(f"band must be 'NR' or 'ER', got {band!r}")

    # Tabulate the noiseless band centroid Er -> (Ep, Eq) and invert by
    # interpolation (both coordinates are monotone in Er).
    er_tab = np.geomspace(1e-3, er_max, 4000)
    y_tab = y_a * er_tab ** y_b
    ep_tab = er_tab * (1.0 + y_tab * V / (1000.0 * eps))
    eq_tab = y_tab * er_tab

    def eq_ridge(ep):
        return np.interp(ep, ep_tab, eq_tab)

    def inner_points_func(ep, ymin, ymax):
        eqr = float(eq_ridge(ep))
        # local band width in Eq: charge resolution plus the phonon
        # resolution mapped through the ridge slope
        dep = max(1e-3, 0.01 * ep)
        slope = (float(eq_ridge(ep + dep)) - float(eq_ridge(ep - dep))) / (2 * dep)
        sigp_val = float(ppq.sigp(ep, eps=eps, V=V, p0=p0, p10=p10))
        sigq_val = float(ppq.sigq(eqr, q0=q0, q10=q10))
        width = np.hypot(sigq_val, slope * sigp_val)
        return [eqr - n_window_widths * width, eqr,
                eqr + n_window_widths * width]

    return inner_points_func
