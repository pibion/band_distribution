"""
Python wrapper exposing the Fortran PpqN / PpqG PDFs as vectorized
pdf_func(Ep_flat, Eq_flat) -> values_flat callables.

All physics parameters are required keyword arguments — there are no
defaults, so the caller must specify the model fully.

The shared library evaluates points sequentially (~3-25 ms per point,
scaling with max(Ep, Eq) because the internal Er integration range grows
with energy).  Since ctypes releases the GIL during foreign calls, this
wrapper parallelizes large requests by splitting them into chunks and
dispatching PpqN_vector calls across a thread pool.  Points are
interleaved across chunks after sorting by max(Ep, Eq) so each chunk
carries a similar share of expensive (high-energy) points.
"""

import ctypes
import os
from concurrent.futures import ThreadPoolExecutor

import numpy as np

_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_LIB_DIR = os.path.join(_REPO_ROOT, "lib")

_DOUBLE_ARR = np.ctypeslib.ndpointer(dtype=np.float64, flags="C_CONTIGUOUS")

_api = None


def _load_library():
    """Load libband_distribution.so.

    The library has no SONAME and its RUNPATH points at nix store paths,
    so its dependencies (libassert.so, libjulienne.so) can only be found
    via LD_LIBRARY_PATH, which must be set before the Python process
    starts:  LD_LIBRARY_PATH=lib python your_script.py
    """
    global _api
    if _api is not None:
        return _api

    try:
        _api = ctypes.CDLL(os.path.join(_LIB_DIR, "libband_distribution.so"))
    except OSError as e:
        raise RuntimeError(
            f"Could not load libband_distribution.so ({e}). "
            f"Run with LD_LIBRARY_PATH={_LIB_DIR} set before Python starts, "
            f"e.g.: LD_LIBRARY_PATH=lib python your_script.py"
        ) from e

    _api.PpqN_vector.argtypes = [
        _DOUBLE_ARR, _DOUBLE_ARR, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double,
        _DOUBLE_ARR,
    ]
    _api.PpqN_vector.restype = None

    _api.PpqG_vector.argtypes = [
        _DOUBLE_ARR, _DOUBLE_ARR, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        _DOUBLE_ARR,
    ]
    _api.PpqG_vector.restype = None

    return _api


def make_ppqn_pdf(*, a, b, F0, s, eps, V, p0, p10, q0, q10, n_workers=None):
    """
    Build pdf_func(Ep_flat, Eq_flat) -> values_flat backed by the Fortran
    PpqN (nuclear recoil band PDF).

    All physics parameters are required:
      a, b           : ionization yield  Y(Er) = a * |Er|^b
      F0, s          : Fano factor       F(Er) = F0 + s * Er
      eps            : energy per e/h pair [keV]
      V              : bias voltage [V]
      p0, p10        : phonon resolution parameters
      q0, q10        : charge resolution parameters

    n_workers : int or None
        Threads used for large requests.  Defaults to os.cpu_count().
        Requests of <= 16 points bypass the pool (e.g. scalar calls from
        scipy quadrature).
    """
    api = _load_library()
    if n_workers is None:
        n_workers = os.cpu_count()
    scalars = (a, b, F0, s, eps, V, p0, p10, q0, q10)

    def _eval_chunk(ep, eq, out):
        api.PpqN_vector(ep, eq, ep.size, *scalars, out)

    return _make_threaded_pdf(_eval_chunk, n_workers)


def make_ppqg_pdf(*, F0, s, eps, V, p0, p10, q0, q10, n_workers=None):
    """Same as make_ppqn_pdf but for the Fortran PpqG (gamma / ER band).
    The yield is fixed at Y = 1 internally, so a and b are not taken."""
    api = _load_library()
    if n_workers is None:
        n_workers = os.cpu_count()
    scalars = (F0, s, eps, V, p0, p10, q0, q10)

    def _eval_chunk(ep, eq, out):
        api.PpqG_vector(ep, eq, ep.size, *scalars, out)

    return _make_threaded_pdf(_eval_chunk, n_workers)


def _make_threaded_pdf(eval_chunk, n_workers):
    def pdf(ep_flat, eq_flat):
        ep = np.ascontiguousarray(np.asarray(ep_flat, dtype=np.float64).ravel())
        eq = np.ascontiguousarray(np.asarray(eq_flat, dtype=np.float64).ravel())
        n = ep.size
        out = np.empty(n, dtype=np.float64)

        if n <= 16 or n_workers <= 1:
            eval_chunk(ep, eq, out)
        else:
            # Cost per point scales with max(Ep, Eq): sort by cost and deal
            # points round-robin so every chunk gets a balanced mix.
            order = np.argsort(-np.maximum(ep, eq), kind="stable")
            n_chunks = min(n_workers * 4, n)
            chunks = [order[i::n_chunks] for i in range(n_chunks)]

            def run(idx):
                ep_c = np.ascontiguousarray(ep[idx])
                eq_c = np.ascontiguousarray(eq[idx])
                out_c = np.empty(idx.size, dtype=np.float64)
                eval_chunk(ep_c, eq_c, out_c)
                out[idx] = out_c

            with ThreadPoolExecutor(max_workers=n_workers) as pool:
                list(pool.map(run, chunks))

        return np.nan_to_num(out, nan=0.0, posinf=0.0, neginf=0.0)

    return pdf
