"""
Sample (Ep, Eq) from an arbitrary numerically-defined 2D PDF.

The user specifies:
  - a PDF callable
  - a closed, bounded 2D region  (ep_range, eq_range)
  - a grid resolution            (n_ep, n_eq points per axis)

build_sampler() returns a lightweight callable that generates independent
(Ep, Eq) pairs using lintsampler's linear interpolant method.

Grid spacing
------------
lintsampler evaluates the PDF at every cell vertex and fits a bilinear
function within each cell, so accuracy is controlled by having cells
small enough that the PDF is approximately linear within them.

For PDFs whose feature width is roughly constant, uniform edges
(build_sampler) are appropriate.  For PDFs whose feature width varies
with position (e.g. detector bands whose resolution grows with energy),
use make_feature_scaled_edges() to space the edges proportionally to the
local feature width, then build_sampler_from_edges().  This keeps the
number of cells per feature width constant everywhere.

What does NOT work is quantile (equal-probability) spacing of the
marginals: it makes tail cells wide regardless of the local feature
width there, which causes the bilinear interpolant to over-sample the
tails and inflate variance estimates.  make_quantile_edges() is provided
only for constructing chi-square histogram bins, never for the sampling
grid.
"""

import numpy as np
from scipy.integrate import cumulative_trapezoid
from lintsampler import DensityGrid, LintSampler


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _evaluate_pdf_on_grid(pdf_func, ep_centers, eq_centers):
    """Evaluate pdf_func on a 2-D meshgrid; return (M, K) array."""
    Ep, Eq = np.meshgrid(
        np.asarray(ep_centers, dtype=np.float64),
        np.asarray(eq_centers, dtype=np.float64),
        indexing="ij",
    )
    vals = pdf_func(Ep.ravel(), Eq.ravel())
    return np.asarray(vals, dtype=np.float64).reshape(Ep.shape)


def _make_lint_pdf(pdf_func):
    """
    Wrap pdf_func(Ep_flat, Eq_flat) -> values_flat into the signature
    lintsampler expects when vectorizedpdf=True: f(x) where x has shape
    (..., 2) and the return has shape (...).
    """
    def lint_pdf(x):
        x = np.asarray(x, dtype=np.float64)
        ep = x[..., 0].ravel()
        eq = x[..., 1].ravel()
        vals = pdf_func(ep, eq)
        return np.asarray(vals, dtype=np.float64).reshape(x.shape[:-1])
    return lint_pdf


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def build_sampler_from_edges(pdf_func, ep_edges, eq_edges):
    """
    Evaluate the PDF on the given rectilinear grid and return a sampler.

    The edges need not be evenly spaced (lintsampler supports arbitrary
    monotonically increasing edges); see make_feature_scaled_edges() for
    constructing edges matched to a position-dependent feature width.

    The PDF is evaluated once at construction time at all grid vertices.
    The returned sampler draws independent samples cheaply using
    lintsampler's bilinear interpolant within each cell.

    Parameters
    ----------
    pdf_func : callable
        pdf_func(Ep_flat, Eq_flat) -> values_flat
        Accepts two 1-D float64 arrays and returns a 1-D float64 array.
    ep_edges, eq_edges : 1-D arrays
        Monotonically increasing cell edges along each axis.

    Returns
    -------
    sampler : callable
        sampler(n_events, *, seed=None) -> (Ep, Eq)
        Ep and Eq are independent 1-D float64 arrays of length n_events.
    """
    density_grid = DensityGrid(
        edges=(np.asarray(ep_edges, dtype=np.float64),
               np.asarray(eq_edges, dtype=np.float64)),
        pdf=_make_lint_pdf(pdf_func),
        vectorizedpdf=True,
    )

    def sampler(n_events, *, seed=None):
        ls = LintSampler(density_grid, seed=seed)
        samples = ls.sample(N=n_events)   # shape (n_events, 2)
        return samples[:, 0], samples[:, 1]

    return sampler


def build_sampler(pdf_func, ep_range, eq_range, n_ep, n_eq):
    """
    Evaluate the PDF on a uniform grid and return a sampler.

    Convenience wrapper around build_sampler_from_edges() with uniformly
    spaced edges.

    Parameters
    ----------
    pdf_func : callable
        pdf_func(Ep_flat, Eq_flat) -> values_flat
        Accepts two 1-D float64 arrays and returns a 1-D float64 array.
    ep_range : (float, float)
        Lower and upper bound of the sampling region along Ep.
        Should be wide enough to contain all significant PDF mass —
        probability density outside this region is silently ignored.
    eq_range : (float, float)
        Lower and upper bound of the sampling region along Eq.
    n_ep : int
        Number of grid points along the Ep axis.  Increase until the
        sampled moments are stable.
    n_eq : int
        Number of grid points along the Eq axis.

    Returns
    -------
    sampler : callable
        sampler(n_events, *, seed=None) -> (Ep, Eq)
    """
    return build_sampler_from_edges(
        pdf_func,
        np.linspace(ep_range[0], ep_range[1], n_ep),
        np.linspace(eq_range[0], eq_range[1], n_eq),
    )


def make_feature_scaled_edges(width_func, lo, hi, cells_per_width, *,
                              n_fine=20001):
    """
    Construct cell edges whose spacing tracks a position-dependent
    feature width, giving a constant number of cells per feature width.

    Solves ds = dE / width(E):  edges are uniform in s, so a feature of
    local width w(E) always spans ~cells_per_width cells regardless of
    where it sits.  For detector bands this means passing the local
    resolution function (e.g. sigp or sigq) as width_func.

    Parameters
    ----------
    width_func : callable
        width_func(E_array) -> local feature width at each E (same units
        as E, must be positive over [lo, hi]).
    lo, hi : float
        Bounds of the axis.
    cells_per_width : float
        Number of cells spanning one local feature width.  Larger values
        give a more accurate (and more expensive) sampling grid; ~5 or
        more is a reasonable starting point for the bilinear interpolant.
    n_fine : int
        Resolution of the internal grid used to integrate 1/width.

    Returns
    -------
    edges : 1-D array, length n_cells + 1
        Monotonically increasing edges from lo to hi.
    """
    if lo > 0:
        fine = np.geomspace(lo, hi, n_fine)
    else:
        fine = np.linspace(lo, hi, n_fine)

    w = np.asarray(width_func(fine), dtype=np.float64)
    if np.any(w <= 0):
        raise ValueError("width_func must be positive over [lo, hi]")

    # s(E) = integral of dE / w(E): "distance" measured in feature widths
    s = cumulative_trapezoid(1.0 / w, fine, initial=0.0)
    n_cells = int(np.ceil(cells_per_width * s[-1]))

    edges = np.interp(np.linspace(0.0, s[-1], n_cells + 1), s, fine)
    edges[0], edges[-1] = lo, hi
    return edges


def make_vertex_lookup_pdf(ep_edges, eq_edges, values):
    """
    Wrap precomputed vertex values as a pdf_func for build_sampler_from_edges.

    Useful when the PDF is expensive: evaluate it once on the grid
    vertices (possibly cached to disk), then hand the stored values to
    the sampler without re-evaluating.  The returned callable only
    accepts coordinates that exactly match grid vertices.

    Parameters
    ----------
    ep_edges, eq_edges : 1-D arrays
        The same edges that will be passed to build_sampler_from_edges.
    values : 2-D array, shape (len(ep_edges), len(eq_edges))
        PDF values at every vertex: values[i, j] = pdf(ep_edges[i], eq_edges[j]).

    Returns
    -------
    pdf_func : callable
        pdf_func(Ep_flat, Eq_flat) -> values_flat
    """
    ep_edges = np.asarray(ep_edges, dtype=np.float64)
    eq_edges = np.asarray(eq_edges, dtype=np.float64)
    values = np.asarray(values, dtype=np.float64)
    if values.shape != (ep_edges.size, eq_edges.size):
        raise ValueError(
            f"values shape {values.shape} != ({ep_edges.size}, {eq_edges.size})"
        )

    def pdf_func(ep_flat, eq_flat):
        ep = np.asarray(ep_flat, dtype=np.float64).ravel()
        eq = np.asarray(eq_flat, dtype=np.float64).ravel()
        i = np.searchsorted(ep_edges, ep)
        j = np.searchsorted(eq_edges, eq)
        i = np.clip(i, 0, ep_edges.size - 1)
        j = np.clip(j, 0, eq_edges.size - 1)
        if not (np.array_equal(ep_edges[i], ep) and np.array_equal(eq_edges[j], eq)):
            raise ValueError(
                "make_vertex_lookup_pdf was queried off-vertex; it only "
                "supports exact grid-vertex coordinates."
            )
        return values[i, j]

    return pdf_func


def make_quantile_edges(pdf_func, ep_range, eq_range, n_ep_bins, n_eq_bins, *,
                        n_coarse=500):
    """
    Determine histogram bin edges with approximately equal marginal
    probability per bin.

    Useful for constructing chi-square bins where expected counts per bin
    are roughly uniform.  Not suitable as a sampling grid (see module
    docstring).

    Parameters
    ----------
    pdf_func : callable
    ep_range, eq_range : (float, float)
    n_ep_bins, n_eq_bins : int   number of bins (edges = n_bins + 1)
    n_coarse : int               resolution of the internal estimation grid

    Returns
    -------
    ep_edges : ndarray, length n_ep_bins + 1
    eq_edges : ndarray, length n_eq_bins + 1
    """
    ep_fine = np.linspace(ep_range[0], ep_range[1], n_coarse)
    eq_fine = np.linspace(eq_range[0], eq_range[1], n_coarse)
    dep = ep_fine[1] - ep_fine[0]
    deq = eq_fine[1] - eq_fine[0]

    grid = np.clip(_evaluate_pdf_on_grid(pdf_func, ep_fine, eq_fine), 0.0, None)

    def _edges(centers, marginal, n_bins):
        cdf = np.cumsum(marginal)
        cdf /= cdf[-1]
        edges = np.interp(np.linspace(0.0, 1.0, n_bins + 1), cdf, centers)
        edges[0], edges[-1] = centers[0], centers[-1]
        return edges

    return (
        _edges(ep_fine, grid.sum(axis=1) * deq, n_ep_bins),
        _edges(eq_fine, grid.sum(axis=0) * dep, n_eq_bins),
    )
