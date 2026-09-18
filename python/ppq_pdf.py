"""
PpqPDF: the fixed context for one fit/MCMC run -- an (Ep, Eq) region, the
observed dataset, and the normalization-integral's convergence tolerance
-- exposing, per parameter point, the un-normalized PDF at the bound
dataset, the region's normalization integral, and their ratio (normalized
PDF values), for both the NR (PpqN) and ER (PpqG) bands.

All physics parameters (k, Z, F0, eps, V, p0, p10, q0, q10) are passed
fresh to every method call rather than bound at construction: those are
what an MCMC step actually varies, while the region, dataset, and
tolerance don't change across a run.

Built entirely on the already-validated primitives in ppqfort_pdf.py --
ppqn_region/ppqg_region (-> Fortran PpqN_region/PpqG_region, a doubling-
verified nested Gauss-Legendre quadrature) and make_ppqn_pdf/make_ppqg_pdf
(-> Fortran PpqN_vector/PpqG_vector) -- no new Fortran code and no new
numerics.

For MCMC, pass precomputed normalization tables (python/normgrid.py,
built once on a batch system) as ppqn_table/ppqg_table: the integral then
costs tens of microseconds instead of ~2 s, and raises rather than ever
extrapolating outside the table's parameter box.
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
    norm_epsrel, norm_epsabs : float or None
        Convergence tolerance for the normalization integral's doubling-
        verified quadrature (see PpqN_region's doc comment in
        src/PpqFort_m.f90) -- unrelated to the size of ep_data/eq_data,
        hence the norm_ prefix.  Required for any band without a table.
    ppqn_table, ppqg_table : path, normgrid.NormInterpolator, or None
        Precomputed normalization table for that band (built by
        python/normgrid.py for exactly this region).  If given, ppqn_integral/
        ppqg_integral interpolate it instead of running the quadrature.
    n_workers : int or None
        Threads used for large batched PpqN_vector/PpqG_vector calls;
        see make_ppqn_pdf.
    """

    def __init__(self, ep_min, ep_max, eq_min, eq_max, ep_data, eq_data, *,
                 norm_epsrel=None, norm_epsabs=None,
                 ppqn_table=None, ppqg_table=None,
                 n_workers=None):
        self.ep_min, self.ep_max = ep_min, ep_max
        self.eq_min, self.eq_max = eq_min, eq_max
        self.ep_data = np.ascontiguousarray(np.asarray(ep_data, dtype=np.float64))
        self.eq_data = np.ascontiguousarray(np.asarray(eq_data, dtype=np.float64))
        self.norm_epsrel = norm_epsrel
        self.norm_epsabs = norm_epsabs
        self.n_workers = n_workers
        self.ppqn_table = self._load_table(ppqn_table, "NR")
        self.ppqg_table = self._load_table(ppqg_table, "ER")

    def _load_table(self, table, band):
        if table is None:
            return None
        import normgrid
        if not isinstance(table, normgrid.NormInterpolator):
            table = normgrid.NormInterpolator.from_hdf5(table)
        if table.band != band:
            raise ValueError(f"table is for the {table.band} band, expected {band}")
        region = (self.ep_min, self.ep_max, self.eq_min, self.eq_max)
        if not np.allclose(table.region, region, rtol=0, atol=1e-12):
            raise ValueError(f"table was built for region {table.region}, not {region}")
        return table

    def _quadrature_tolerances(self):
        if self.norm_epsrel is None or self.norm_epsabs is None:
            raise ValueError("norm_epsrel and norm_epsabs are required to compute a normalization "
                             "integral without a precomputed table")
        return dict(epsrel=self.norm_epsrel, epsabs=self.norm_epsabs)

    # ---- NR (PpqN) band ----

    def ppqn_integral(self, *, k, Z, F0, eps, V, p0, p10, q0, q10):
        """Normalization: integral of PpqN over this region."""
        if self.ppqn_table is not None:
            return self.ppqn_table(k=k, Z=Z, F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
        return ppqn_region(self.ep_min, self.ep_max, self.eq_min, self.eq_max,
                            **self._quadrature_tolerances(),
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
        if self.ppqg_table is not None:
            return self.ppqg_table(F0=F0, eps=eps, V=V, p0=p0, p10=p10, q0=q0, q10=q10)
        return ppqg_region(self.ep_min, self.ep_max, self.eq_min, self.eq_max,
                            **self._quadrature_tolerances(),
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
