submodule(PpqFort_m) PpqFort_s
  implicit none

  ! Simpson 1/3 weights for 21 equally-spaced points
  real(c_double), parameter :: simps_w(21) = &
      [1.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, &
       4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 2.0d0, 4.0d0, 1.0d0]

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

  module procedure PpqN
      real(c_double), parameter :: minx = 5E-21_c_double
      integer :: i, j
      real(c_double) :: integral, norm_const
      real(c_double) :: sp0, sp1, sq0, sq1   ! inlined sigp/sigq coefficients
      real(c_double) :: Er_i, F_val, Nbar_val, sigma_N
      real(c_double) :: EP_mean, EQ_mean, sigp_mean_sq, sigq_mean_sq
      real(c_double) :: a_coeff, b_coeff, N_star, sigma_Neff
      real(c_double) :: N_lo, N_hi, h_N
      real(c_double) :: N_j, EP_nl, EQ_nl, sig_p, sig_q, exp_arg
      real(c_double) :: integrand_arr(21), contrib_i

      if (F0 == 0.0d0 .and. s == 0.0d0) &
          error stop "Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero"

      norm_const = 1.0d0 / ((2.0d0 * pi)**1.5d0)
      ! Decompose sigp(EP) = sqrt(sp0 + sp1*EP^2), sigq(EQ) = sqrt(sq0 + sq1*EQ^2)
      sp0 = p0**2
      sp1 = (p10**2 - p0**2) / (10.0d0 * (1.0d0 + V / (eps * 1.0d3)))**2
      sq0 = q0**2
      sq1 = (q10**2 - q0**2) / 100.0d0

      associate(maxx => max(1.1_c_double * max(Ep, Eq), max(Ep, Eq) + 10.0_c_double))
        associate(resolution => merge(0.002_c_double, 0.01_c_double, maxx < 15.0_c_double))
          integrate_PpqFullN: &
          associate(npts => int((maxx - minx) / resolution) + 1)
            define_recoil_energies: &
            associate(er_arr => [(minx + (i - 1) * resolution, i = 1, npts)])
              integral = 0.0d0
              do concurrent(i = 1:npts) default(none) reduce(+: integral) &
                  shared(er_arr, a, b, Ep, Eq, F0, s, eps, V, norm_const, &
                         sp0, sp1, sq0, sq1) &
                  local(Er_i, F_val, Nbar_val, sigma_N, EP_mean, EQ_mean, &
                        sigp_mean_sq, sigq_mean_sq, a_coeff, b_coeff, N_star, sigma_Neff, &
                        N_lo, N_hi, h_N, N_j, EP_nl, EQ_nl, sig_p, sig_q, &
                        exp_arg, integrand_arr, contrib_i, j)
                  Er_i     = er_arr(i)
                  F_val    = F0 + s * Er_i
                  Nbar_val = a * abs(Er_i)**b * Er_i / eps
                  if (Nbar_val <= 0.0d0 .or. F_val <= 0.0d0) then
                      contrib_i = 0.0d0
                  else
                      sigma_N      = sqrt(Nbar_val * F_val)
                      EP_mean      = Er_i + (V * 1.0d-3) * Nbar_val
                      EQ_mean      = eps * Nbar_val
                      sigp_mean_sq = sp0 + sp1 * EP_mean**2
                      sigq_mean_sq = sq0 + sq1 * EQ_mean**2
                      a_coeff = (V * 1.0d-3) * (Ep - Er_i) / sigp_mean_sq &
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
                          EP_nl = Er_i + (V * 1.0d-3) * N_j
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
                      contrib_i = (h_N / 3.0d0) * dot_product(simps_w, integrand_arr) * PErN(Er_i)
                  end if
                  integral = integral + contrib_i
              end do
              res = integral * resolution
            end associate define_recoil_energies
          end associate integrate_PpqFullN
        end associate
      end associate
  end procedure PpqN

  module procedure PpqG
      real(c_double), parameter :: minx = 5E-21_c_double
      integer :: i, j
      real(c_double) :: integral, norm_const, resolution, maxx
      integer :: npts
      real(c_double), allocatable :: er_arr(:)
      real(c_double) :: sp0, sp1, sq0, sq1   ! inlined sigp/sigq coefficients
      real(c_double) :: Er_i, F_val, Nbar_val, sigma_N
      real(c_double) :: EP_mean, EQ_mean, sigp_mean_sq, sigq_mean_sq
      real(c_double) :: a_coeff, b_coeff, N_star, sigma_Neff
      real(c_double) :: N_lo, N_hi, h_N
      real(c_double) :: N_j, EP_nl, EQ_nl, sig_p, sig_q, exp_arg
      real(c_double) :: integrand_arr(21), contrib_i

      if (F0 == 0.0d0 .and. s == 0.0d0) &
          error stop "Fano factor F(Er) = F0 + s*Er is zero for all Er; F0 and s cannot both be zero"

      norm_const = 1.0d0 / ((2.0d0 * pi)**1.5d0)
      sp0 = p0**2
      sp1 = (p10**2 - p0**2) / (10.0d0 * (1.0d0 + V / (eps * 1.0d3)))**2
      sq0 = q0**2
      sq1 = (q10**2 - q0**2) / 100.0d0

      maxx       = max(1.1_c_double * max(Ep, Eq), max(Ep, Eq) + 10.0_c_double)
      resolution = merge(0.002d0, 0.01d0, maxx < 15.0d0)
      npts       = int((maxx - minx) / resolution) + 1

      allocate(er_arr(npts))
      do i = 1, npts
          er_arr(i) = minx + (i - 1) * resolution
      end do

      integral = 0.0d0

      do concurrent(i = 1:npts) default(none) reduce(+: integral) &
          shared(er_arr, Ep, Eq, F0, s, eps, V, norm_const, sp0, sp1, sq0, sq1) &
          local(Er_i, F_val, Nbar_val, sigma_N, EP_mean, EQ_mean, &
                sigp_mean_sq, sigq_mean_sq, a_coeff, b_coeff, N_star, sigma_Neff, &
                N_lo, N_hi, h_N, N_j, EP_nl, EQ_nl, sig_p, sig_q, &
                exp_arg, integrand_arr, contrib_i, j)
          Er_i     = er_arr(i)
          F_val    = F0 + s * Er_i
          Nbar_val = Er_i / eps          ! Y=1 for electron recoils
          if (Nbar_val <= 0.0d0 .or. F_val <= 0.0d0) then
              contrib_i = 0.0d0
          else
              sigma_N      = sqrt(Nbar_val * F_val)
              EP_mean      = Er_i + (V * 1.0d-3) * Nbar_val
              EQ_mean      = eps * Nbar_val          ! = Er_i for ER
              sigp_mean_sq = sp0 + sp1 * EP_mean**2
              sigq_mean_sq = sq0 + sq1 * EQ_mean**2
              a_coeff = (V * 1.0d-3) * (Ep - Er_i) / sigp_mean_sq &
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
                  EP_nl = Er_i + (V * 1.0d-3) * N_j
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
              contrib_i = (h_N / 3.0d0) * dot_product(simps_w, integrand_arr) * PErG(Er_i)
          end if
          integral = integral + contrib_i
      end do

      res = integral * resolution
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
  ! Evaluates the N integral numerically via 21-point simps from stdlib.
  ! sigp and sigq are evaluated at exact noiseless energies (Er+V*N/1000,
  ! eps*N) per sample point, matching the data simulator exactly.
  module procedure PpqFullN
      real(c_double) :: F_val, Nbar_val, sigma_N, norm_const
      real(c_double) :: sp0, sp1, sq0, sq1
      real(c_double) :: EP_mean, EQ_mean, sigp_mean_sq, sigq_mean_sq
      real(c_double) :: a_coeff, b_coeff, N_star, sigma_Neff
      real(c_double) :: N_lo, N_hi, h_N
      real(c_double) :: N_j, EP_nl, EQ_nl, sig_p, sig_q, exp_arg
      real(c_double) :: integrand_arr(21)
      integer :: j

      F_val    = F(Er, F0, s)
      Nbar_val = Nbar(Er, a, b, eps)
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
      res = (h_N / 3.0d0) * dot_product(simps_w, integrand_arr) * PErN(Er)
  end procedure PpqFullN

end submodule PpqFort_s
