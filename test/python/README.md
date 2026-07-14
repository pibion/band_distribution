# Python tests

Everything in this directory validates the band-distribution PDFs; the code
that *computes* the PDFs lives in `python/` at the repository root
(`pq_dist_v8.py` for the python reference implementation, `ppqfort_pdf.py`
for the wrapper around the compiled Fortran library).  Run all scripts from
the repository root; see the main README for expected timings and for how to
build `lib/libband_distribution.so` first.

## Test scripts

| Script | What it checks |
|--------|----------------|
| `test_PpqFort.py` | Fortran functions against the python reference implementation, one value at a time (machine-precision parity). |
| `test_PpqFort_vectorFuncs.py` | The vectorized `PpqN_vector` / `PpqG_vector` entry points. |
| `verify_sample_from_pdf.py` | Moments of the grid sampler used by the self-consistency tests. |
| `test_chisquare_gaussian.py` | The chi-square harness itself, against an analytic 2D Gaussian (no Fortran involved). |
| `test_chisquare_ppqn.py` | Self-consistency: events sampled *from* the Fortran `PpqN` follow `PpqN`. |
| `test_chisquare_nr_simulator.py` | Physics: events from the independent NR simulator follow the Fortran `PpqN`. |
| `test_chisquare_er_simulator.py` | Physics: events from the independent ER simulator follow the Fortran `PpqG`. |

The two simulator tests are the strongest statement: the event generator
(`generate_events.py`) knows nothing about the PDF's integrals — it draws
Er from the recoil spectrum, N from a truncated normal, and applies the
detector response — so agreement means the PDFs describe the physics.

## Supporting modules

| Module | Role |
|--------|------|
| `chisquare_harness.py` | PDF-agnostic machinery: equibin binning, expected counts by per-bin quadrature, batched chi-square throws. |
| `band_breakpoints.py` | Band physics for the per-bin quadrature: where the band ridge is and how wide it is, packaged as breakpoints for the harness. |
| `generate_events.py` | The independent physics simulator (truth generator) for the two simulator tests. |
| `sample_from_pdf.py` | Grid-based sampling from an arbitrary PDF (lintsampler), used by the self-consistency tests. |
| `quad_breakpoints_demo.py` | Self-contained demonstration (no repo dependencies) of why adaptive quadrature silently misses narrow peaks and how `points=` fixes it. |

## Integrating the PDF over bins

The expected count in a bin is the PDF integrated over the bin rectangle.
This is the subtlest part of the tests, because the band is a narrow ridge:
adaptive quadrature over a wide tail bin can sample only the flat ~zero
region, agree with itself that the integral is zero, and never subdivide
toward the peak — a silent wrong answer, not an error.

Suggested reading order:

1. `quad_breakpoints_demo.py` — run it and watch plain `quad` return ~0 for
   a narrow peak, then fix it with the `points=` argument.
2. `band_breakpoints.py` — where the ridge is physically: the noiseless
   band centroid Er → (Ep, Eq) from the yield model and Neganov-Luke
   amplification, inverted to give the most probable Eq at each Ep, plus
   the local band width setting the breakpoint window.
3. `chisquare_harness.py`, `expected_counts_from_pdf` — the nested
   quadrature that consumes those breakpoints (`inner_points_func`) and
   normalizes the per-bin integrals into expected counts.
