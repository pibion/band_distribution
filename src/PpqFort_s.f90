submodule(PpqFort_m) PpqFort_s
  implicit none

  ! Simpson 1/3 weights for 21 equally-spaced points
  real(c_double), parameter :: simps_w(21) = &
      [1.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, &
       4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 1.0d0]

  ! Lower edge of the Er integration domain (keV)
  real(c_double), parameter :: er_min = 1.0d-6

  ! Shortest decay length in the gamma recoil spectrum PErG (keV); must
  ! stay in sync with PGb inside the PErG module procedure.  The Er
  ! sampling step is capped at a fraction of this scale when the
  ! integration window reaches down to the low-Er spectrum spike.
  real(c_double), parameter :: PGb_scale = 0.169520023d0

  ! The Er integrand can be multimodal (typically bimodal, usually
  ! dominated by the first, sharp peak).  Rather than assume dominance,
  ! every local maximum of the closed-form log-integrand gets its own
  ! integration window; overlapping windows are merged.
  integer, parameter :: max_modes = 4

contains

  module procedure PpqN_vector
      integer :: i

      do concurrent(i = 1:n) default (none) shared(res_arr, Ep_arr, Eq_arr, a, b, F0, s, eps, V, p0, p10, q0, q10)
          res_arr(i) = PpqN(Ep_arr(i), Eq_arr(i), a, b, F0, s, eps, V, p0, p10, q0, q10)
      end do
  end procedure PpqN_vector

  module procedure PpqG_vector
      integer :: i

      do concurrent(i = 1:n) default (none) shared(res_arr, Ep_arr, Eq_arr, F0, s, eps, V, p0, p10, q0, q10)
          res_arr(i) = PpqG(Ep_arr(i), Eq_arr(i), F0, s, eps, V, p0, p10, q0, q10)
      end do
  end procedure PpqG_vector

  ! PpqN, PpqG: integrate P(Ep,Eq|Er)*P(Er) dEr.
  !
  ! The integrand is a narrow bump in Er: the phonon and ionization
  ! measurements each constrain Er to within a few sigma, so sampling the
  ! full physical range wastes almost every evaluation.  Instead, locate
  ! the peak of the closed-form log-integrand (Laplace approximation) and
  ! integrate only over peak +/- 8 sigma with a uniform trapezoid rule.
  ! For a smooth integrand that vanishes at both window edges the uniform
  ! trapezoid rule converges spectrally, so ~100 points reproduce the old
  ! full-range fine grid (~35,000 points) to near machine precision.

  module procedure PpqN
      if (F0 == 0.0d0 .and. s == 0.0d0) &
          error stop "Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero"
      res = integrate_band(Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, .false.)
  end procedure PpqN

  module procedure PpqG
      if (F0 == 0.0d0 .and. s == 0.0d0) &
          error stop "Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero"
      res = integrate_band(Ep, Eq, 1.0d0, 0.0d0, F0, s, eps, V, p0, p10, q0, q10, .true.)
  end procedure PpqG

  module procedure Y
      res = a * abs(Er)**b
  end procedure Y

  module procedure F
      res = F0 + s * Er
  end procedure F

  module procedure Nbar
      res = Y(Er, a, b) * Er / eps
  end procedure Nbar

  module procedure sigp
      real(c_double) :: e2pre, e2e

      e2pre = (p10**2 - p0**2)
      e2e = (Ep / (10.0d0 * (1.0d0 + (V / eps / 1000.0d0))))**2
      res = sqrt(p0**2 + e2pre * e2e)
  end procedure sigp

  module procedure sigq
      real(c_double) :: e2pre, e2e

      e2pre = (q10**2 - q0**2)
      e2e = (Eq / 10.0d0)**2
      res = sqrt(q0**2 + e2pre * e2e)
  end procedure sigq

  module procedure PErN
      real(c_double), parameter :: PNa = 0.53693208d0
      real(c_double), parameter :: PNb = 6.41515782d0
      real(c_double), parameter :: PNd = 23.71789286d0
      real(c_double), parameter :: PNa_over_PNb   = PNa / PNb
      real(c_double), parameter :: one_minus_PNa_over_PNd = (1.0d0 - PNa) / PNd

      res = PNa_over_PNb * exp(-Er / PNb) + one_minus_PNa_over_PNd * exp(-Er / PNd)
  end procedure PErN

  module procedure PErG
      real(c_double), parameter :: PGa = 0.573211975d0
      real(c_double), parameter :: PGb = 0.169520023d0
      real(c_double), parameter :: PGd = 279.552394d0
      real(c_double), parameter :: PGa_over_PGb   = PGa / PGb
      real(c_double), parameter :: one_minus_PGa_over_PGd = (1.0d0 - PGa) / PGd

      res = PGa_over_PGb * exp(-Er / PGb) + one_minus_PGa_over_PGd * exp(-Er / PGd)
  end procedure PErG

  ! PpqFullN: integrand for the Er integral, P(Ep,Eq|Er)*P(Er), NR case.
  ! Evaluates the N integral numerically via 21-point Simpson's rule.
  ! sigp and sigq are evaluated at exact noiseless energies (Er+V*N/1000,
  ! eps*N) per sample point, matching the data simulator exactly.
  module procedure PpqFullN
      res = PpqCondN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) * PErN(Er)
  end procedure PpqFullN

  ! ---------------------------------------------------------------------
  ! Private helpers
  ! ---------------------------------------------------------------------

  ! P(Ep,Eq|Er): the 21-point Simpson N integral, without the recoil
  ! spectrum factor.  The gamma / electron-recoil case is a=1, b=0.
  pure function PpqCondN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res)
      real(c_double), intent(in) :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
      real(c_double) :: res
      real(c_double) :: F_val, Nbar_val, sigma_N, norm_const
      real(c_double) :: sp0, sp1, sq0, sq1
      real(c_double) :: EP_mean, EQ_mean, sigp_mean_sq, sigq_mean_sq
      real(c_double) :: a_coeff, b_coeff, N_star, sigma_Neff
      real(c_double) :: N_lo, N_hi, h_N
      real(c_double) :: N_j, EP_nl, EQ_nl, sig_p, sig_q, exp_arg
      real(c_double) :: integrand_arr(21)
      integer :: j

      F_val    = F0 + s * Er
      Nbar_val = a * abs(Er)**b * Er / eps
      if (Nbar_val <= 0.0d0 .or. F_val <= 0.0d0) then
          res = 0.0d0
          return
      end if
      sigma_N    = sqrt(Nbar_val * F_val)
      norm_const = 1.0d0 / ((2.0d0 * pi)**1.5d0)
      sp0 = p0**2
      sp1 = (p10**2 - p0**2) / (10.0d0 * (1.0d0 + V / (eps * 1.0d3)))**2
      sq0 = q0**2
      sq1 = (q10**2 - q0**2) / 100.0d0

      EP_mean      = Er + (V * 1.0d-3) * Nbar_val
      EQ_mean      = eps * Nbar_val
      sigp_mean_sq = sp0 + sp1 * EP_mean**2
      sigq_mean_sq = sq0 + sq1 * EQ_mean**2
      a_coeff = (V * 1.0d-3) * (Ep - Er) / sigp_mean_sq &
              + eps * Eq / sigq_mean_sq &
              + 1.0d0 / F_val
      b_coeff = 1.0d0 / (2.0d0 * Nbar_val * F_val) &
              + eps**2 / (2.0d0 * sigq_mean_sq) &
              + V**2 / (2.0d0 * 1.0d6 * sigp_mean_sq)
      N_star     = a_coeff / (2.0d0 * b_coeff)
      sigma_Neff = 1.0d0 / sqrt(2.0d0 * b_coeff)
      if (N_star - 5.0d0*sigma_Neff < 0.0d0) then
          N_lo = 0.0d0
          N_hi = 10.0d0 * sigma_Neff
      else
          N_lo = N_star - 5.0d0 * sigma_Neff
          N_hi = N_star + 5.0d0 * sigma_Neff
      end if
      h_N = (N_hi - N_lo) / 20.0d0
      do j = 1, 21
          N_j   = N_lo + (j - 1) * h_N
          EP_nl = Er + (V * 1.0d-3) * N_j
          EQ_nl = eps * N_j
          sig_p = sqrt(sp0 + sp1 * EP_nl**2)
          sig_q = sqrt(sq0 + sq1 * EQ_nl**2)
          exp_arg = -0.5d0*((Ep - EP_nl)/sig_p)**2 &
                    - 0.5d0*((Eq - EQ_nl)/sig_q)**2 &
                    - 0.5d0*((N_j - Nbar_val)/sigma_N)**2
          if (exp_arg < -700.0d0) then
              integrand_arr(j) = 0.0d0
          else
              integrand_arr(j) = norm_const / (sigma_N * sig_p * sig_q) * exp(exp_arg)
          end if
      end do
      res = (h_N / 3.0d0) * dot_product(simps_w, integrand_arr)
  end function PpqCondN

  ! The full integrand P(Ep,Eq|Er)*P(Er) with the spectrum selected by
  ! is_gamma.  Pass a=1, b=0 together with is_gamma=.true. for the ER band.
  pure function band_integrand(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma) result(res)
      real(c_double), intent(in) :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
      logical, intent(in) :: is_gamma
      real(c_double) :: res

      res = PpqCondN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
      if (is_gamma) then
          res = res * PErG(Er)
      else
          res = res * PErN(Er)
      end if
  end function band_integrand

  ! Closed-form approximation to log of the Er integrand: complete the
  ! square in N analytically (with sigp/sigq frozen at the measured
  ! energies) and add the log recoil spectrum.  Not accurate enough for
  ! the PDF value itself -- the real integrand evaluates sigp/sigq at
  ! noiseless energies -- but the peak location and curvature agree to
  ! well within a fraction of the peak width, which is all the window
  ! placement needs.
  pure function log_integrand_approx(Er, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma) result(H)
      real(c_double), intent(in) :: Er, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2
      logical, intent(in) :: is_gamma
      real(c_double) :: H
      real(c_double) :: F_val, Nbar_val, aN, bN, cN, spec

      if (Er <= 0.0d0) then
          H = -1.0d30
          return
      end if
      F_val    = F0 + s * Er
      Nbar_val = a * Er**b * Er / eps
      if (Nbar_val <= 0.0d0 .or. F_val <= 0.0d0) then
          H = -1.0d30
          return
      end if
      aN = (V * 1.0d-3) * (Ep - Er) / sp2 + eps * Eq / sq2 + 1.0d0 / F_val
      bN = 1.0d0 / (2.0d0 * Nbar_val * F_val) + eps**2 / (2.0d0 * sq2) &
         + V**2 / (2.0d0 * 1.0d6 * sp2)
      cN = -(Ep - Er)**2 / (2.0d0 * sp2) - Eq**2 / (2.0d0 * sq2) &
         - Nbar_val / (2.0d0 * F_val)
      if (is_gamma) then
          spec = PErG(Er)
      else
          spec = PErN(Er)
      end if
      H = cN + aN**2 / (4.0d0 * bN) + log(spec)
  end function log_integrand_approx

  ! Locate every mode of the Er integrand and estimate their widths.
  !
  ! An unconditional log-spaced scan of the closed-form log-integrand
  ! (cheap: no N integral) finds all candidate local maxima -- log
  ! spacing so the sharp low-Er peak cannot fall between scan points.
  ! Each candidate is refined with Newton iterations; sigma = 1/sqrt(-H'')
  ! from the curvature at the refined peak.  A candidate at the low-Er
  ! boundary with the integrand monotonically decaying away from it is
  ! kept as a boundary mode whose sigma is the exponential decay length.
  ! Modes more than ~100 log units below the best one cannot contribute
  ! and are dropped; near-duplicate refinements are merged.
  !
  ! n_modes = 0 signals failure; the caller falls back to the original
  ! full-range grid.
  pure subroutine locate_modes(Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma, &
                               er_m, sig_m, H_m, bnd_m, n_modes)
      real(c_double), intent(in) :: Ep, Eq, a, b, F0, s, eps, V, sp2, sq2
      logical, intent(in) :: is_gamma
      real(c_double), intent(out) :: er_m(max_modes), sig_m(max_modes), H_m(max_modes)
      logical, intent(out) :: bnd_m(max_modes)
      integer, intent(out) :: n_modes

      integer, parameter :: n_scan = 512
      real(c_double), parameter :: x_scan_min = 1.0d-4
      real(c_double) :: xs(n_scan), Hs(n_scan)
      real(c_double) :: maxx_scan, dlx, x, d, H1, H2, sig, Hv
      integer :: k, m, j, ncand, cand(n_scan), best
      logical :: ok, duplicate

      n_modes = 0
      maxx_scan = max(1.1d0 * max(Ep, Eq), max(Ep, Eq) + 10.0d0)

      dlx = (log(maxx_scan) - log(x_scan_min)) / (n_scan - 1)
      do k = 1, n_scan
          xs(k) = exp(log(x_scan_min) + (k - 1) * dlx)
          Hs(k) = log_integrand_approx(xs(k), Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
      end do

      ! Candidate local maxima (including the low-Er boundary)
      ncand = 0
      do k = 1, n_scan
          if (Hs(k) < -1.0d29) cycle
          if ((k == 1      .or. Hs(k) >= Hs(k - 1)) .and. &
              (k == n_scan .or. Hs(k) >  Hs(k + 1))) then
              ncand = ncand + 1
              cand(ncand) = k
          end if
      end do
      if (ncand == 0) return

      ! Sort candidates by scan height, tallest first (ncand is small)
      do m = 1, ncand - 1
          best = m
          do j = m + 1, ncand
              if (Hs(cand(j)) > Hs(cand(best))) best = j
          end do
          k = cand(m); cand(m) = cand(best); cand(best) = k
      end do

      refine_candidates: do m = 1, ncand
          if (n_modes == max_modes) exit
          ! Anything this far below the tallest mode contributes nothing
          if (n_modes > 0) then
              if (Hs(cand(m)) < maxval(H_m(1:n_modes)) - 100.0d0) cycle
          end if

          k = cand(m)
          if (k == 1) then
              ! Boundary mode: integrand decays away from the low-Er edge
              x = xs(1)
              d = 1.0d-5
              H1 = (log_integrand_approx(x + d, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma) &
                  - Hs(1)) / d
              if (H1 >= 0.0d0) cycle
              n_modes = n_modes + 1
              er_m(n_modes)  = x
              sig_m(n_modes) = min(-1.0d0 / H1, maxx_scan)
              H_m(n_modes)   = Hs(1)
              bnd_m(n_modes) = .true.
              cycle
          end if

          x = xs(k)
          call newton_refine(x, ok, H1, H2)
          if (.not. ok) then
              ! Accept the scan point with finite-difference curvature if usable
              x = xs(k)
              d = max(1.0d-5, 1.0d-5 * x)
              H2 = (log_integrand_approx(x + d, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma) &
                  - 2.0d0 * Hs(k) &
                  + log_integrand_approx(x - d, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)) / d**2
              if (H2 >= 0.0d0) cycle
          end if
          sig = 1.0d0 / sqrt(-H2)
          Hv  = log_integrand_approx(x, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
          if (sig /= sig .or. sig <= 1.0d-6 .or. sig > maxx_scan) cycle

          ! Drop refinements that landed on an already-recorded mode
          duplicate = .false.
          do j = 1, n_modes
              if (abs(x - er_m(j)) < 0.5d0 * (sig + sig_m(j))) then
                  duplicate = .true.
                  exit
              end if
          end do
          if (duplicate) cycle

          n_modes = n_modes + 1
          er_m(n_modes)  = x
          sig_m(n_modes) = sig
          H_m(n_modes)   = Hv
          bnd_m(n_modes) = .false.
      end do refine_candidates

  contains

      pure subroutine newton_refine(xx, converged, H1_out, H2_out)
          real(c_double), intent(inout) :: xx
          logical, intent(out) :: converged
          real(c_double), intent(out) :: H1_out, H2_out
          integer :: it
          real(c_double) :: dd, Hm, H0, Hp, step

          converged = .false.
          H1_out = 0.0d0
          H2_out = 0.0d0
          do it = 1, 60
              dd = max(1.0d-4, 1.0d-4 * xx)
              Hm = log_integrand_approx(xx - dd, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
              H0 = log_integrand_approx(xx,      Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
              Hp = log_integrand_approx(xx + dd, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
              H1_out = (Hp - Hm) / (2.0d0 * dd)
              H2_out = (Hp - 2.0d0 * H0 + Hm) / dd**2
              if (H2_out >= 0.0d0) return
              step = -H1_out / H2_out
              step = sign(min(abs(step), 0.5d0 * xx + 1.0d0), step)
              if (xx + step <= 0.0d0) then
                  xx = 0.5d0 * xx
              else
                  xx = xx + step
              end if
              if (xx > 1.5d0 * maxx_scan) return
              if (abs(step) < 1.0d-4) then
                  converged = .true.
                  return
              end if
          end do
      end subroutine newton_refine

  end subroutine locate_modes

  ! Trapezoid rule over [lo, lo + (npts-1)*h].  For a smooth integrand
  ! vanishing at both edges this is spectrally accurate on a uniform grid.
  pure function integrate_window(lo, h, npts, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma) result(total)
      real(c_double), intent(in) :: lo, h, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
      integer, intent(in) :: npts
      logical, intent(in) :: is_gamma
      real(c_double) :: total
      integer :: i
      real(c_double) :: f_i, w_i

      total = 0.0d0
      do concurrent(i = 1:npts) default(none) reduce(+: total) &
          shared(lo, h, npts, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma) &
          local(f_i, w_i)
          f_i = band_integrand(lo + (i - 1) * h, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)
          w_i = merge(0.5d0, 1.0d0, i == 1 .or. i == npts)
          total = total + w_i * f_i
      end do
      total = total * h
  end function integrate_window

  ! Shared driver for PpqN / PpqG: locate every integrand mode, place an
  ! integration window around each, merge overlapping windows, and guard
  ! the window edges.
  pure function integrate_band(Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma) result(res)
      real(c_double), intent(in) :: Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
      logical, intent(in) :: is_gamma
      real(c_double) :: res

      integer, parameter :: max_pts = 200000
      real(c_double) :: sp2, sq2, maxx_scan, resolution
      real(c_double) :: er_m(max_modes), sig_m(max_modes), H_m(max_modes)
      logical :: bnd_m(max_modes)
      real(c_double) :: w_lo(max_modes), w_hi(max_modes), w_ht(max_modes), w_sig(max_modes)
      real(c_double) :: h, f_star, f_lo, f_hi, lo_limit, hi_limit, tmp
      integer :: n_modes, n_win, m, j, npts, guard
      logical :: need_lo, need_hi

      sp2 = sigp(Ep, eps, V, p0, p10)**2
      sq2 = sigq(Eq, q0, q10)**2

      call locate_modes(Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma, &
                        er_m, sig_m, H_m, bnd_m, n_modes)

      if (n_modes == 0) then
          ! Last resort: the original fine grid over the full range.
          ! Guaranteed correct, rarely taken.
          maxx_scan  = max(1.1d0 * max(Ep, Eq), max(Ep, Eq) + 10.0d0)
          resolution = merge(0.002d0, 0.01d0, maxx_scan < 15.0d0)
          npts = int((maxx_scan - er_min) / resolution) + 1
          res  = integrate_window(er_min, resolution, npts, Ep, Eq, a, b, F0, s, eps, V, &
                                  p0, p10, q0, q10, is_gamma)
          return
      end if

      ! Predicted peak log-value so small the integral underflows to zero
      if (maxval(H_m(1:n_modes)) < -600.0d0) then
          res = 0.0d0
          return
      end if

      ! One window per mode: +/- 8 sigma around an interior peak; for a
      ! boundary decay, 35 decay lengths puts the truncated tail below
      ! 1e-15 relative.
      do m = 1, n_modes
          if (bnd_m(m)) then
              w_lo(m) = er_min
              w_hi(m) = er_m(m) + 35.0d0 * sig_m(m)
          else
              w_lo(m) = max(er_min, er_m(m) - 8.0d0 * sig_m(m))
              w_hi(m) = er_m(m) + 8.0d0 * sig_m(m)
          end if
          w_ht(m)  = sig_m(m) / 6.0d0
          w_sig(m) = sig_m(m)
      end do

      ! Sort windows by lower edge (n_modes <= 4)
      do m = 1, n_modes - 1
          do j = 1, n_modes - m
              if (w_lo(j) > w_lo(j + 1)) then
                  tmp = w_lo(j);  w_lo(j)  = w_lo(j + 1);  w_lo(j + 1)  = tmp
                  tmp = w_hi(j);  w_hi(j)  = w_hi(j + 1);  w_hi(j + 1)  = tmp
                  tmp = w_ht(j);  w_ht(j)  = w_ht(j + 1);  w_ht(j + 1)  = tmp
                  tmp = w_sig(j); w_sig(j) = w_sig(j + 1); w_sig(j + 1) = tmp
              end if
          end do
      end do

      ! Merge overlapping windows, keeping the finer step
      n_win = 1
      do m = 2, n_modes
          if (w_lo(m) <= w_hi(n_win)) then
              w_hi(n_win)  = max(w_hi(n_win), w_hi(m))
              w_ht(n_win)  = min(w_ht(n_win), w_ht(m))
              w_sig(n_win) = min(w_sig(n_win), w_sig(m))
          else
              n_win = n_win + 1
              w_lo(n_win)  = w_lo(m)
              w_hi(n_win)  = w_hi(m)
              w_ht(n_win)  = w_ht(m)
              w_sig(n_win) = w_sig(m)
          end if
      end do

      ! Edge-guard threshold: the tallest mode sets the scale
      f_star = 0.0d0
      do m = 1, n_modes
          f_star = max(f_star, band_integrand(er_m(m), Ep, Eq, a, b, F0, s, eps, V, &
                                              p0, p10, q0, q10, is_gamma))
      end do

      res = 0.0d0
      do m = 1, n_win
          ! Resolve the sharp low-Er component of the gamma spectrum when
          ! the window reaches down into it
          if (is_gamma .and. w_lo(m) < 5.0d0 * PGb_scale) w_ht(m) = min(w_ht(m), PGb_scale / 6.0d0)

          ! Guard expansions must not cross into neighbouring windows:
          ! adjacent edges share a point, whose two half-weights sum to
          ! the full trapezoid weight of the union.
          lo_limit = er_min
          if (m > 1) lo_limit = w_hi(m - 1)
          hi_limit = huge(1.0d0)
          if (m < n_win) hi_limit = w_lo(m + 1)
          w_lo(m) = max(w_lo(m), lo_limit)

          do guard = 1, 4
              npts = min(max(int((w_hi(m) - w_lo(m)) / w_ht(m)) + 2, 25), max_pts)
              h    = (w_hi(m) - w_lo(m)) / (npts - 1)
              tmp  = integrate_window(w_lo(m), h, npts, Ep, Eq, a, b, F0, s, eps, V, &
                                      p0, p10, q0, q10, is_gamma)
              ! Edge guard: expand if the integrand has not died off at the
              ! window edges (protects against an underestimated sigma)
              f_lo = band_integrand(w_lo(m), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)
              f_hi = band_integrand(w_hi(m), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)
              need_lo = f_lo > 1.0d-10 * f_star .and. w_lo(m) > lo_limit
              need_hi = f_hi > 1.0d-10 * f_star .and. w_hi(m) < hi_limit
              if (.not. (need_lo .or. need_hi)) exit
              if (need_lo) w_lo(m) = max(lo_limit, w_lo(m) - 4.0d0 * w_sig(m))
              if (need_hi) w_hi(m) = min(hi_limit, w_hi(m) + 4.0d0 * w_sig(m))
          end do
          res = res + tmp
      end do
  end function integrate_band

end submodule PpqFort_s
