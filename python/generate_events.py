"""
Generate (Ep, Eq) samples from the PpqN / PpqG probability densities using
the hierarchical generative model that underlies the analytic PDF.

The analytic PDF PpqN(Ep, Eq) is obtained by:
  1. Marginalizing N analytically (→ the erf factor in PpqFullN)
  2. Integrating over Er numerically

Generating samples instead follows the latent-variable chain:
  Er ~ PErN               (biexponential neutron recoil spectrum)
  N  | Er ~ TruncNormal   (Fano-factor fluctuation, truncated at N >= 0)
  Ep | Er,N ~ Normal      (phonon detector: Neganov-Luke + thermal)
  Eq | N    ~ Normal      (charge detector)

For gamma (ER) events, substitute PErG and Y=1 (a=1, b=0).
"""

import numpy as np
from scipy.stats import truncnorm
import sys, os

sys.path.insert(0, os.path.dirname(__file__))
import pq_dist_v6 as ppq


def _sample_PErN(n, rng, PNa=0.53693208, PNb=6.41515782, PNd=23.71789286):
    """Sample Er from the biexponential neutron spectrum PErN."""
    component = rng.uniform(size=n) < PNa
    return np.where(component,
                    rng.exponential(PNb, size=n),
                    rng.exponential(PNd, size=n))


def _sample_PErG(n, rng, PGa=0.573211975, PGb=0.169520023, PGd=279.552394):
    """Sample Er from the biexponential gamma spectrum PErG."""
    component = rng.uniform(size=n) < PGa
    return np.where(component,
                    rng.exponential(PGb, size=n),
                    rng.exponential(PGd, size=n))


def _sample_N(Nbar_arr, F_arr, rng):
    """
    Sample N (number of e/h pairs) for each event from a half-normal truncated
    at 0.  Mean = Nbar, Var = Nbar * F.

    Uses scipy.stats.truncnorm vectorised over events.
    """
    sigma_N = np.sqrt(np.maximum(Nbar_arr * F_arr, 0.0))

    # scalar path: avoid divide-by-zero when sigma is zero
    N = np.empty_like(Nbar_arr)
    zero_mask = sigma_N == 0.0
    N[zero_mask] = Nbar_arr[zero_mask]

    nonzero = ~zero_mask
    if nonzero.any():
        a_trunc = -Nbar_arr[nonzero] / sigma_N[nonzero]   # lower clip in std units
        b_trunc = np.inf
        N[nonzero] = truncnorm.rvs(a_trunc, b_trunc,
                                   loc=Nbar_arr[nonzero],
                                   scale=sigma_N[nonzero],
                                   random_state=rng)
    return N


def generate_NR_events(n_events, *,
                       a=0.16, b=0.18,
                       F0=0.122, s=0.0,
                       eps=3.0e-3,
                       V=3.0,
                       p0=0.06421907, p10=0.48998486,
                       q0=0.23718488, q10=0.27093151,
                       PNa=0.53693208, PNb=6.41515782, PNd=23.71789286,
                       seed=None):
    """
    Generate (Ep, Eq) pairs distributed according to PpqN.

    Returns
    -------
    Ep : ndarray, shape (n_events,)
    Eq : ndarray, shape (n_events,)
    Er : ndarray, shape (n_events,)   – latent true recoil energies
    N  : ndarray, shape (n_events,)   – latent e/h pair counts
    """
    rng = np.random.default_rng(seed)

    Er = _sample_PErN(n_events, rng, PNa=PNa, PNb=PNb, PNd=PNd)

    Nbar = ppq.Y(Er, a=a, b=b) * Er / eps          # mean e/h pairs
    F_val = F0 + s * Er                              # Fano factor

    N = _sample_N(Nbar, F_val, rng)

    Ep_true = Er + N * (V / 1e3)                    # Neganov-Luke + recoil
    Eq_true = N * eps                                # ionization

    # Detector smearing: evaluate resolution at the true (noise-free) signal.
    # sigp and sigq vary slowly, so this closely approximates the heteroscedastic
    # PDF where resolutions are evaluated at the measured values.
    sigp_val = ppq.sigp(Ep_true, eps=eps, V=V, p0=p0, p10=p10)
    sigq_val = ppq.sigq(Eq_true, q0=q0, q10=q10)

    Ep = rng.normal(Ep_true, sigp_val)
    Eq = rng.normal(Eq_true, sigq_val)

    return Ep, Eq, Er, N


def generate_ER_events(n_events, *,
                       F0=0.122, s=0.0,
                       eps=3.0e-3,
                       V=3.0,
                       p0=0.06421907, p10=0.48998486,
                       q0=0.23718488, q10=0.27093151,
                       PGa=0.573211975, PGb=0.169520023, PGd=279.552394,
                       seed=None):
    """
    Generate (Ep, Eq) pairs distributed according to PpqG.

    Gamma / electron-recoil events have ionization yield Y = 1 (a=1, b=0).

    Returns
    -------
    Ep : ndarray, shape (n_events,)
    Eq : ndarray, shape (n_events,)
    Er : ndarray, shape (n_events,)
    N  : ndarray, shape (n_events,)
    """
    return generate_NR_events(n_events,
                              a=1.0, b=0.0,
                              F0=F0, s=s,
                              eps=eps, V=V,
                              p0=p0, p10=p10,
                              q0=q0, q10=q10,
                              PNa=PGa, PNb=PGb, PNd=PGd,
                              seed=seed)
