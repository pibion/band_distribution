submodule(PpqFort_m) PpqFort_s
  implicit none

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
      integer :: i
      real(c_double) :: integral, sigp_val, sigq_val
      real(c_double) :: Er_i, F_val, Nbar_val, aN_val, bN_val, cN1_val, cN2_val
      real(c_double) :: exponent, prefac_i, erffac_i, contrib_i

      associate(maxx => max(1.1 * max(Ep, Eq), max(Ep, Eq) + 10))
        associate(resolution => merge(0.002_c_double, 0.01_c_double, maxx < 15))
          integrate_PpqFullN: &
          associate(npts => int((maxx - minx) / resolution) + 1)
            define_recoil_energies: &
            associate(er_arr => [(minx + (i - 1) * resolution, i = 1, npts)])
              sigp_val = sigp(Ep, eps, V, p0, p10)
              sigq_val = sigq(Eq, q0, q10)
#if ! CANNOT_DO_CONCURRENT
              integral = 0.0d0
              do concurrent(i = 1:npts) default(none) reduce(+: integral) &
                  shared(er_arr, a, b, Ep, Eq, F0, s, eps, V, sigp_val, sigq_val) &
                  local(Er_i, F_val, Nbar_val, aN_val, bN_val, cN1_val, cN2_val, exponent, prefac_i, erffac_i, contrib_i)
                  Er_i     = er_arr(i)
                  F_val    = F0 + s * Er_i
                  Nbar_val = a * abs(Er_i)**b * Er_i / eps
                  aN_val = (V/1.0d3)*(Ep - Er_i)/sigp_val**2 + eps*Eq/sigq_val**2 + 1.0d0/F_val
                  bN_val = 1.0d0/(2.0d0*Nbar_val*F_val) + eps**2/(2.0d0*sigq_val**2) &
                         + V**2/(2.0d0*(sigp_val*1.0d3)**2)
                  cN1_val = -((Ep - Er_i)**2)/(2.0d0*sigp_val**2) - (Eq**2)/(2.0d0*sigq_val**2) &
                          - Nbar_val/(2.0d0*F_val)
                  exponent = cN1_val + (aN_val**2)/(4.0d0*bN_val)
                  if (exponent < -700.0d0) then
                    contrib_i = 0.0d0
                  else
                    cN2_val = 1.0d0/(2.0d0*pi*sqrt(2.0d0*pi)*sigp_val*sigq_val*sqrt(Nbar_val*F_val))
                    prefac_i = 0.5d0*sqrt(pi)*cN2_val*exp(exponent)/sqrt(bN_val)
                    if (prefac_i > 0.0d0) then
                      erffac_i  = erf(aN_val/(2.0d0*sqrt(bN_val))) + 1.0d0
                      contrib_i = prefac_i * erffac_i * PErN(Er_i)
                    else
                      contrib_i = 0.0d0
                    end if
                  end if
                  integral = integral + contrib_i
              end do
#else
              integral = sum([(PpqFullN(er_arr(i), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10), i = 1, npts)])
#endif
              res = integral * resolution
            end associate define_recoil_energies
          end associate integrate_PpqFullN
        end associate
      end associate
  end procedure PpqN

  module procedure PpqG
      real(c_double), allocatable :: er_arr(:)
      real(c_double) :: minx, maxx, resolution, integral
      integer :: npts, i
      real(c_double) :: sigp_val, sigq_val
      real(c_double) :: Er_i, F_val, Nbar_val, aN_val, bN_val, cN1_val, cN2_val
      real(c_double) :: exponent, prefac_i, erffac_i, contrib_i

      minx = 5E-21
      maxx = max(1.1 * max(Ep, Eq), max(Ep, Eq) + 10)
      resolution = merge(0.002d0, 0.01d0, maxx < 15.0d0)
      npts = int((maxx - minx) / resolution) + 1

      allocate(er_arr(npts))
      do i = 1, npts
          er_arr(i) = minx + (i - 1) * resolution
      end do

      sigp_val = sigp(Ep, eps, V, p0, p10)
      sigq_val = sigq(Eq, q0, q10)

      integral = 0.0d0

      do concurrent(i = 1:npts) default(none) reduce(+: integral) &
          shared(er_arr, Ep, Eq, F0, s, eps, V, sigp_val, sigq_val) &
          local(Er_i, F_val, Nbar_val, aN_val, bN_val, cN1_val, cN2_val, exponent, prefac_i, erffac_i, contrib_i)
          Er_i     = er_arr(i)
          F_val    = F0 + s * Er_i
          Nbar_val = Er_i / eps                  ! Y=1 for gamma (a=1, b=0)
          aN_val = (V/1.0d3)*(Ep - Er_i)/sigp_val**2 + eps*Eq/sigq_val**2 + 1.0d0/F_val
          bN_val = 1.0d0/(2.0d0*Nbar_val*F_val) + eps**2/(2.0d0*sigq_val**2) &
                 + V**2/(2.0d0*(sigp_val*1.0d3)**2)
          cN1_val = -((Ep - Er_i)**2)/(2.0d0*sigp_val**2) - (Eq**2)/(2.0d0*sigq_val**2) &
                  - Nbar_val/(2.0d0*F_val)
          exponent = cN1_val + (aN_val**2)/(4.0d0*bN_val)
          if (exponent < -700.0d0) then
            contrib_i = 0.0d0
          else
            cN2_val = 1.0d0/(2.0d0*pi*sqrt(2.0d0*pi)*sigp_val*sigq_val*sqrt(Nbar_val*F_val))
            prefac_i = 0.5d0*sqrt(pi)*cN2_val*exp(exponent)/sqrt(bN_val)
            if (prefac_i > 0.0d0) then
              erffac_i  = erf(aN_val/(2.0d0*sqrt(bN_val))) + 1.0d0
              contrib_i = prefac_i * erffac_i * PErG(Er_i)
            else
              contrib_i = 0.0d0
            end if
          end if
          integral = integral + contrib_i
      end do

      res = integral * resolution
  end procedure PpqG

  module procedure Y
      ! Compute the ionization yield
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
      ! determined by fits to Geant output
      real(c_double), parameter :: PNa = 0.53693208d0
      real(c_double), parameter :: PNb = 6.41515782d0
      real(c_double), parameter :: PNd = 23.71789286d0
      real(c_double), parameter :: PNa_over_PNb   = PNa / PNb
      real(c_double), parameter :: one_minus_PNa_over_PNd = (1.0d0 - PNa) / PNd

      res = PNa_over_PNb * exp(-Er / PNb) + one_minus_PNa_over_PNd * exp(-Er / PNd)
  end procedure PErN

  module procedure PErG
      ! determined by fits to Geant output
      real(c_double), parameter :: PGa = 0.573211975d0
      real(c_double), parameter :: PGb = 0.169520023d0
      real(c_double), parameter :: PGd = 279.552394d0
      real(c_double), parameter :: PGa_over_PGb   = PGa / PGb
      real(c_double), parameter :: one_minus_PGa_over_PGd = (1.0d0 - PGa) / PGd

      res = PGa_over_PGb * exp(-Er / PGb) + one_minus_PGa_over_PGd * exp(-Er / PGd)
  end procedure PErG

  module procedure PpqFullN
      real(c_double) :: sigp_val, sigq_val, F_val, Nbar_val
      real(c_double) :: aN_val, bN_val, cN1_val, cN2_val
      real(c_double) :: exponent, prefac, erffac

      sigp_val = sigp(Ep, eps, V, p0, p10)
      sigq_val = sigq(Eq, q0, q10)
      F_val    = F(Er, F0, s)
      Nbar_val = Nbar(Er, a, b, eps)

      aN_val = (V/1.0d3)*(Ep - Er) / sigp_val**2 &
             + eps*Eq / sigq_val**2 &
             + 1.0d0 / F_val

      bN_val = 1.0d0 / (2.0d0 * Nbar_val * F_val) &
             + eps**2 / (2.0d0 * sigq_val**2) &
             + V**2 / (2.0d0 * (sigp_val * 1.0d3)**2)

      cN1_val = -((Ep - Er)**2) / (2.0d0 * sigp_val**2) &
              - (Eq**2) / (2.0d0 * sigq_val**2) &
              - Nbar_val / (2.0d0 * F_val)

      exponent = cN1_val + (aN_val**2) / (4.0d0 * bN_val)

      if (exponent < -700.0d0) then
        res = 0.0d0
        return
      end if

      cN2_val = 1.0d0 / (2.0d0 * pi * sqrt(2.0d0 * pi) &
              * sigp_val * sigq_val * sqrt(Nbar_val * F_val))

      prefac = 0.5d0 * sqrt(pi) * cN2_val * exp(exponent) / sqrt(bN_val)

      if (prefac > 0.0d0) then
        erffac = erf(aN_val / (2.0d0 * sqrt(bN_val))) + 1.0d0
        res = prefac * erffac * PErN(Er)
      else
        res = 0.0d0
      end if
  end procedure PpqFullN
end submodule PpqFort_s
