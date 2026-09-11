"""
PpqPDF: the fixed context for one fit/MCMC run -- an (Ep, Eq) region, the
observed dataset, and the normalization-integral's quadrature grid
density -- exposing, per parameter point, the un-normalized PDF at the
bound dataset, the region's normalization integral, and their ratio
(normalized PDF values), for both the NR (PpqN) and ER (PpqG) bands.

All physics parameters (k, Z, F0, eps, V, p0, p10, q0, q10) are passed
fresh to every method call rather than bound at construction: those are
what an MCMC step actually varies, while the region, dataset, and grid
density don't change across a run.

Built entirely on the already-validated primitives in ppqfort_pdf.py --
ppqn_region/ppqg_region (-> Fortran PpqN_region/PpqG_region) and
make_ppqn_pdf/make_ppqg_pdf (-> Fortran PpqN_vector/PpqG_vector) -- no new
Fortran code and no new numerics.  region_integral.py is the independent
Python-side reimplementation used to cross-check the region-integral
algorithm; see test/python/test_region_integral.py.
"""

import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from ppqfort_pdf import make_ppqg_pdf, make_ppqn_pdf, ppqg_region, ppqn_region


class PpqPDF:
    """
    ep_min, ep_max, eq_min, eq_max : float
        The (Ep, Eq) region the PDF is normalized over.
    ep_data, eq_data : array_like
        The fixed observed dataset this fit is evaluated against.
    norm_n_ep, norm_n_eq_window, norm_n_window_widths :
        Quadrature grid density for the normalization integral (see
        PpqN_region's doc comment in src/PpqFort_m.f90) -- unrelated to
        the size of ep_data/eq_data, hence the norm_ prefix: n_ep is not
        the number of data points.
    n_workers : int or None
        Threads used for large batched PpqN_vector/PpqG_vector calls;
        see make_ppqn_pdf.
    """

    def __init__(self, ep_min, ep_max, eq_min, eq_max, ep_data, eq_data, *,
                 norm_n_ep, norm_n_eq_window, norm_n_window_widths,
                 n_workers=None):
        self.ep_min, self.ep_max = ep_min, ep_max
        self.eq_min, self.eq_max = eq_min, eq_max
        self.ep_data = np.ascontiguousarray(np.asarray(ep_data, dtype=np.float64))
        self.eq_data = np.ascontiguousarray(np.asarray(eq_data, dtype=np.float64))
        self.norm_n_ep = norm_n_ep
        self.norm_n_eq_window = norm_n_eq_window
        self.norm_n_window_widths = norm_n_window_widths
        self.n_workers = n_workers

    # ---- NR (PpqN) band ----

    def ppqn_integral(self, *, k, Z, F0, eps, V, p0, p10, q0, q10):
        """Normalization: integral of PpqN over this region."""
        return ppqn_region(self.ep_min, self.ep_max, self.eq_min, self.eq_max,
                            n_ep=self.norm_n_ep, n_eq_window=self.norm_n_eq_window,
                            n_window_widths=self.norm_n_window_widths,
                            k=k, Z=Z, F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)

    def ppqn_values(self, *, k, Z, F0, eps, V, p0, p10, q0, q10):
        """Un-normalized PpqN at the bound dataset."""
        pdf_func = make_ppqn_pdf(k=k, Z=Z, F0=F0, eps=eps, V=V, p0=p0, p10=p10,
                                  q0=q0, q10=q10, n_workers=self.n_workers)
        return pdf_func(self.ep_data, self.eq_data)

    def ppqn_normalized_values(self, *, k, Z, F0, eps, V, p0, p10, q0, q10):
        """ppqn_values(...) / ppqn_integral(...), elementwise."""
        return (self.ppqn_values(k=k, Z=Z, F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
                / self.ppqn_integral(k=k, Z=Z, F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10))

    # ---- ER (PpqG) band ----

    def ppqg_integral(self, *, F0, eps, V, p0, p10, q0, q10):
        """Normalization: integral of PpqG over this region."""
        return ppqg_region(self.ep_min, self.ep_max, self.eq_min, self.eq_max,
                            n_ep=self.norm_n_ep, n_eq_window=self.norm_n_eq_window,
                            n_window_widths=self.norm_n_window_widths,
                            F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)

    def ppqg_values(self, *, F0, eps, V, p0, p10, q0, q10):
        """Un-normalized PpqG at the bound dataset."""
        pdf_func = make_ppqg_pdf(F0=F0, eps=eps, V=V, p0=p0, p10=p10,
                                  q0=q0, q10=q10, n_workers=self.n_workers)
        return pdf_func(self.ep_data, self.eq_data)

    def ppqg_normalized_values(self, *, F0, eps, V, p0, p10, q0, q10):
        """ppqg_values(...) / ppqg_integral(...), elementwise."""
        return (self.ppqg_values(F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
                / self.ppqg_integral(F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10))
