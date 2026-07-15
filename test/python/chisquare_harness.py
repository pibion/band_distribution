"""
Harness for testing whether chi-square values from a 2D PDF are
distributed as expected under the null hypothesis.

Workflow
--------
1. Generate a reference dataset; use equibin to define adaptive bins.
   The outer envelope of the bins is snapped to the declared region so
   every thrown event lands in exactly one bin.
2. Compute expected counts per bin by integrating pdf_func over each bin
   with adaptive quadrature, which is accurate regardless of bin size.
3. In batches, generate n_throws independent datasets of n_events each,
   count events in each bin, and accumulate the Pearson chi-square:
       chi2 = sum( (observed - expected)^2 / expected )
4. Return the array of chi-square values for downstream analysis.

The harness is PDF-agnostic: supply any pdf_func(Ep_flat, Eq_flat) ->
values_flat callable and a matching sampler(n_events, seed=None) ->
(Ep, Eq) callable.

Step 2's per-bin quadrature also requires an inner_points_func giving
the band-ridge breakpoints (see band_breakpoints.py); every PDF this
harness is used with in this repo is a narrow ridge, for which plain
adaptive quadrature over a wide bin can silently return ~0 (see
quad_breakpoints_demo.py), so there is no "skip the breakpoints"
option.
"""

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import equibin
from scipy.integrate import quad


# ---------------------------------------------------------------------------
# Expected counts from the PDF
# ---------------------------------------------------------------------------

def expected_counts_from_pdf(pdf_func, bins, n_events, inner_points_func, *,
                             n_workers=1):
    """
    Integrate pdf_func over each bin rectangle using adaptive quadrature
    (scipy quad over Ep, nested quad over Eq).

    Parameters
    ----------
    pdf_func : callable
        pdf_func(Ep_flat, Eq_flat) -> values_flat
    bins : list of (xmin, xmax, ymin, ymax)
        Bin definitions from equibin.BinningResult.bins.
    n_events : int
        Events per throw; scales the returned expected counts.
    inner_points_func : callable
        inner_points_func(ep, ymin, ymax) -> sequence of Eq breakpoints.
        The band PDFs in this repo are narrow ridges: plain adaptive
        quadrature starts from a fixed set of sample points, so a bin
        wider than the ridge can have every sample point land in the
        ~zero region on both sides and quad will confidently report a
        ~zero integral without ever subdividing toward the peak (see
        quad_breakpoints_demo.py).  Passing the ridge location as a
        breakpoint forces a subdivision there, so this argument is
        required rather than optional -- there is no PDF in this repo
        for which skipping it is safe.  Breakpoints outside (ymin, ymax)
        are dropped; an empty result for a bin is fine (scipy's
        points-based quadrature handles zero real breakpoints exactly
        like plain adaptive quadrature).
    n_workers : int
        Threads used to integrate bins concurrently.  Only helps when
        pdf_func releases the GIL (e.g. calls into a compiled library).

    Returns
    -------
    expected : ndarray, shape (n_bins,)
        Expected event count in each bin for a throw of n_events.
    """
    def integrate_bin(bin_rect):
        xmin, xmax, ymin, ymax = bin_rect

        def integrand(eq, ep):
            return float(pdf_func(np.array([ep]), np.array([eq]))[0])

        def inner(ep):
            pts = sorted(p for p in inner_points_func(ep, ymin, ymax)
                        if ymin < p < ymax)
            val, _ = quad(integrand, ymin, ymax, args=(ep,),
                          points=pts, limit=200)
            return val

        val, _ = quad(inner, xmin, xmax, limit=200)
        return max(val, 0.0)

    if n_workers > 1:
        with ThreadPoolExecutor(max_workers=n_workers) as pool:
            integrals = np.array(list(pool.map(integrate_bin, bins)))
    else:
        integrals = np.array([integrate_bin(b) for b in bins])

    total = integrals.sum()
    if total == 0.0:
        raise ValueError("PDF integrates to zero over all bins.")

    return n_events * integrals / total


# ---------------------------------------------------------------------------
# Main harness
# ---------------------------------------------------------------------------

def run_chisquare_test(pdf_func, sampler_func, ep_range, eq_range,
                       n_events, n_throws, n_bins, inner_points_func,
                       *, n_ref=None, n_workers=1,
                       batch_size=5_000, seed=0):
    """
    Run a chi-square distribution test for a 2D PDF.

    Parameters
    ----------
    pdf_func : callable
        pdf_func(Ep_flat, Eq_flat) -> values_flat
    sampler_func : callable
        sampler_func(n, seed=None) -> (Ep, Eq)  — from build_sampler()
    ep_range, eq_range : (float, float)
        Bounds of the sampling region, passed to equibin as hard limits.
    n_events : int
        Events per chi-square throw.  All bins should have
        expected count >> 1; aim for n_events / n_bins >= 10.
    n_throws : int
        Total number of chi-square values to compute (e.g. 100_000).
    n_bins : int
        Target number of equibin bins.
    inner_points_func : callable
        Passed to expected_counts_from_pdf; see its docstring for why
        this is required rather than optional.
    n_ref : int or None
        Reference dataset size for computing bin boundaries.
        Defaults to max(50 * n_bins, 5 * n_events).
    n_workers : int
        Threads for the expected-count integration.
    batch_size : int
        Throws per batch.  Each batch generates batch_size * n_events
        events at once; tune to fit available memory.
    seed : int
        Base random seed.

    Returns
    -------
    chi2_values : ndarray, shape (n_throws,)
    dof : int
        Degrees of freedom (n_bins - 1).
    bins : list of (xmin, xmax, ymin, ymax)
    expected : ndarray, shape (n_bins,)
        Expected counts per throw in each bin.
    """
    rng = np.random.default_rng(seed)

    if n_ref is None:
        n_ref = max(50 * n_bins, 5 * n_events)

    # --- step 1: reference data → equibin bins ---
    print(f"Generating {n_ref:,} reference events for binning...")
    Ep_ref, Eq_ref = sampler_func(n_ref, seed=int(rng.integers(0, 2**31)))
    print(f"Computing {n_bins} equibin bins...")
    result = equibin.bin_2d(
        Ep_ref, Eq_ref, n_bins=n_bins,
        xmin=ep_range[0], xmax=ep_range[1],
        ymin=eq_range[0], ymax=eq_range[1],
    )

    # equibin only tiles the extent of the reference data; snap the outer
    # envelope of the bins to the declared region so the bins tile it
    # fully and every thrown event lands in a bin.
    bins_arr = np.asarray(result.bins, dtype=np.float64)
    bins_arr[bins_arr[:, 0] == bins_arr[:, 0].min(), 0] = ep_range[0]
    bins_arr[bins_arr[:, 1] == bins_arr[:, 1].max(), 1] = ep_range[1]
    bins_arr[bins_arr[:, 2] == bins_arr[:, 2].min(), 2] = eq_range[0]
    bins_arr[bins_arr[:, 3] == bins_arr[:, 3].max(), 3] = eq_range[1]
    bins = [tuple(row) for row in bins_arr]
    n_bins_actual = len(bins)

    # --- step 2: expected counts from PDF ---
    print(f"Integrating PDF over {n_bins_actual} bins "
          f"(adaptive quadrature, {n_workers} workers)...")
    expected = expected_counts_from_pdf(
        pdf_func, bins, n_events, inner_points_func, n_workers=n_workers,
    )

    # --- step 3: batched chi-square throws ---
    print(f"Running {n_throws:,} throws of {n_events:,} events "
          f"(batch_size={batch_size:,})...")
    chi2_values = np.empty(n_throws)

    for batch_start in range(0, n_throws, batch_size):
        batch_end = min(batch_start + batch_size, n_throws)
        b = batch_end - batch_start                  # actual throws this batch

        # generate all events for this batch in one call
        seed_b = int(rng.integers(0, 2**31))
        Ep_b, Eq_b = sampler_func(b * n_events, seed=seed_b)
        Ep_b = Ep_b.reshape(b, n_events)            # (b, n_events)
        Eq_b = Eq_b.reshape(b, n_events)

        # count events per bin per throw: (b, n_bins)
        obs = np.zeros((b, n_bins_actual), dtype=np.intp)
        for j in range(n_bins_actual):
            xmin, xmax, ymin, ymax = bins_arr[j]
            mask = (
                (Ep_b >= xmin) & (Ep_b < xmax) &
                (Eq_b >= ymin) & (Eq_b < ymax)
            )
            obs[:, j] = mask.sum(axis=1)

        # Pearson chi-square: (b,)
        chi2_values[batch_start:batch_end] = (
            ((obs - expected) ** 2 / expected).sum(axis=1)
        )

        pct = 100 * batch_end / n_throws
        print(f"  {batch_end:>8,} / {n_throws:,}  ({pct:.0f}%)", end="\r")

    print()
    return chi2_values, n_bins_actual - 1, bins, expected
