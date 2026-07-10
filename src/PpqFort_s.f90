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

  ! Peak-locator status codes
  integer, parameter :: PEAK_GAUSSIAN = 0   ! interior peak, sigma from curvature
  integer, parameter :: PEAK_BOUNDARY = 1   ! integrand decays from Er ~ 0; sigma is a decay length
  integer, parameter :: PEAK_FAILED   = 2   ! locator failed; integrate the full range

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

  ! Locate the peak of the Er integrand and estimate its width.
  !
  ! Seed: Er0 = Ep - (V/1000)*Eq/eps, i.e. measure N from the ionization
  ! channel and subtract the Luke phonons -- exact at the noiseless band
  ! centroid for any yield model.  Refine with Newton iterations on the
  ! closed-form log-integrand; fall back to a coarse scan if Newton
  ! wanders.  On success sigma_Er = 1/sqrt(-H''); if the integrand is
  ! monotonically decaying from Er ~ 0 (PEAK_BOUNDARY) sigma_Er is the
  ! exponential decay length instead.
  pure subroutine locate_peak(Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma, &
                              Er_star, sigma_Er, H_star, status)
      real(c_double), intent(in) :: Ep, Eq, a, b, F0, s, eps, V, sp2, sq2
      logical, intent(in) :: is_gamma
      real(c_double), intent(out) :: Er_star, sigma_Er, H_star
      integer, intent(out) :: status

      integer, parameter :: n_scan = 512
      real(c_double) :: maxx_scan, x, d, H1, H2, dx, xk, Hk, Hbest
      integer :: k, kbest
      logical :: ok

      maxx_scan = max(1.1d0 * max(Ep, Eq), max(Ep, Eq) + 10.0d0)

      ! Newton from the Luke-subtraction seed
      x = min(max(Ep - (V * 1.0d-3) * (Eq / eps), 1.0d-3), maxx_scan)
      call newton_refine(x, ok, H1, H2)

      if (.not. ok) then
          ! Coarse scan of the closed-form log-integrand (cheap: no N integral)
          dx = (maxx_scan - 1.0d-3) / (n_scan - 1)
          kbest = 1
          Hbest = -huge(1.0d0)
          do k = 1, n_scan
              xk = 1.0d-3 + (k - 1) * dx
              Hk = log_integrand_approx(xk, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
              if (Hk > Hbest) then
                  Hbest = Hk
                  kbest = k
              end if
          end do
          if (kbest == 1) then
              ! Monotonically decaying from the low-Er edge: window by decay length
              x = 1.0d-3
              d = 1.0d-4
              H1 = (log_integrand_approx(x + d, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma) &
                  - log_integrand_approx(x, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)) / d
              if (H1 < 0.0d0) then
                  Er_star  = x
                  sigma_Er = min(-1.0d0 / H1, maxx_scan)
                  H_star   = log_integrand_approx(x, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
                  status   = PEAK_BOUNDARY
              else
                  status = PEAK_FAILED
              end if
              return
          end if
          x = 1.0d-3 + (kbest - 1) * dx
          call newton_refine(x, ok, H1, H2)
          if (.not. ok) then
              status = PEAK_FAILED
              return
          end if
      end if

      Er_star  = x
      sigma_Er = 1.0d0 / sqrt(-H2)
      H_star   = log_integrand_approx(x, Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma)
      if (sigma_Er /= sigma_Er .or. sigma_Er <= 1.0d-6 .or. sigma_Er > maxx_scan) then
          status = PEAK_FAILED
      else
          status = PEAK_GAUSSIAN
      end if

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

  end subroutine locate_peak

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

  ! Shared driver for PpqN / PpqG: locate the integrand peak, place the
  ! integration window around it, and guard the window edges.
  pure function integrate_band(Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma) result(res)
      real(c_double), intent(in) :: Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
      logical, intent(in) :: is_gamma
      real(c_double) :: res

      integer, parameter :: max_pts = 200000
      real(c_double) :: sp2, sq2, Er_star, sigma_Er, H_star
      real(c_double) :: lo, hi, h_target, h, f_star, f_lo, f_hi, maxx_scan, resolution
      integer :: status, npts, guard
      logical :: need_lo, need_hi

      sp2 = sigp(Ep, eps, V, p0, p10)**2
      sq2 = sigq(Eq, q0, q10)**2

      call locate_peak(Ep, Eq, a, b, F0, s, eps, V, sp2, sq2, is_gamma, &
                       Er_star, sigma_Er, H_star, status)

      if (status == PEAK_FAILED) then
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
      if (H_star < -600.0d0) then
          res = 0.0d0
          return
      end if

      ! Window: +/- 8 sigma around an interior peak; for a boundary decay,
      ! 35 decay lengths puts the truncated tail below 1e-15 relative.
      lo = max(er_min, Er_star - 8.0d0 * sigma_Er)
      hi = Er_star + merge(35.0d0, 8.0d0, status == PEAK_BOUNDARY) * sigma_Er

      h_target = sigma_Er / 6.0d0
      ! Resolve the sharp low-Er component of the gamma spectrum when the
      ! window reaches down into it
      if (is_gamma .and. lo < 5.0d0 * PGb_scale) h_target = min(h_target, PGb_scale / 6.0d0)

      f_star = band_integrand(Er_star, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)

      do guard = 1, 4
          npts = min(max(int((hi - lo) / h_target) + 2, 25), max_pts)
          h    = (hi - lo) / (npts - 1)
          res  = integrate_window(lo, h, npts, Ep, Eq, a, b, F0, s, eps, V, &
                                  p0, p10, q0, q10, is_gamma)
          ! Edge guard: expand if the integrand has not died off at the
          ! window edges (protects against an underestimated sigma_Er)
          f_lo = band_integrand(lo, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)
          f_hi = band_integrand(hi, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10, is_gamma)
          need_lo = f_lo > 1.0d-10 * f_star .and. lo > er_min
          need_hi = f_hi > 1.0d-10 * f_star
          if (.not. (need_lo .or. need_hi)) exit
          if (need_lo) lo = max(er_min, lo - 4.0d0 * sigma_Er)
          if (need_hi) hi = hi + 4.0d0 * sigma_Er
      end do
  end function integrate_band

end submodule PpqFort_s
