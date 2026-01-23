submodule(PpqFort_m) PpqFort_s
  implicit none

contains

  module procedure PpqN_vector
      integer :: i
  
      do concurrent(i = 1:n)
          res_arr(i) = PpqN(Ep_arr(i), Eq_arr(i), a, b, F0, s, eps, V, p0, p10, q0, q10)
      end do
  end procedure PpqN_vector

  module procedure PpqG_vector
      integer :: i
  
      do concurrent(i = 1:n)
          res_arr(i) = PpqG(Ep_arr(i), Eq_arr(i), F0, s, eps, V, p0, p10, q0, q10)
      end do
  end procedure PpqG_vector

  module procedure PpqN
      real(c_double), parameter :: minx = 5E-21_c_double
      integer :: i
      real(c_double) integral
  
      associate(maxx => max(1.1 * max(Ep, Eq), max(Ep, Eq) + 10))
        associate(resolution => merge(0.002_c_double, 0.01_c_double, maxx < 15))
          integrate_PpqFullN: &
          associate(npts => int((maxx - minx) / resolution) + 1)
            define_recoil_energies: &
            associate(er_arr => [(minx + (i - 1) * resolution, i = 1, npts)])
#if ! PREFER_DO_CONCURRENT
               integral = sum([(PpqFullN(er_arr(i), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10), i = 1, npts)])
#else
               integral = 0.
               do concurrent(i = 1:npts) default (none) reduce(+: integral) shared(er_arr, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
                   integral = integral + PpqFullN(er_arr(i), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
               end do
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
  
      minx = 5E-21
      maxx = max(1.1 * max(Ep, Eq), max(Ep, Eq) + 10) 
      if (maxx < 15) then
        resolution = 0.002d0
      else
        resolution = 0.01d0
      endif
      npts = int((maxx - minx) / resolution) + 1
  
      ! build the array of energies
      allocate(er_arr(npts))
      do i = 1, npts
          er_arr(i) = minx + (i - 1) * resolution
      end do
  
      integral = 0.0d0
  
      ! Optimizable loop
      do concurrent(i = 1:npts)
          integral = integral + PpqFullG(er_arr(i), Ep, Eq, F0, s, eps, V, p0, p10, q0, q10)
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
      real(c_double) PNa, PNb, PNd
  
      ! we don't ever change these parameters
      ! they were determined by fits to Geant output
      PNa = 0.53693208d0
      PNb = 6.41515782d0
      PNd = 23.71789286d0
  
      if (Er < 0.0d0) then
        res = 0.0d0
      else
        res = PNa * (1.0d0 / PNb) * exp(-Er / PNb) + &
              (1.0d0 - PNa) * (1.0d0 / PNd) * exp(-Er / PNd)
      end if
  end procedure PErN

  module procedure PErG
      real(c_double) :: PGa, PGb, PGd
  
      ! Set values for PGa, PGb, and PGd
      ! we don't ever change these parameters
      ! they were determined by fits to Geant output
      PGa = 0.573211975d0
      PGb = 0.169520023d0
      PGd = 279.552394d0
  
      if (Er < 0.0d0) then
        res = 0.0d0
      else
        res = PGa * (1.0d0 / PGb) * exp(-Er / PGb) + &
              (1.0d0 - PGa) * (1.0d0 / PGd) * exp(-Er / PGd)
      end if
  end procedure PErG

  module procedure aN
      real(c_double) :: t1, t2, t3
      
      t1 = (V/1.0d3)*(Ep - Er) / sigp(Ep, eps, V, p0, p10)**2
      t2 = eps*Eq / sigq(Eq, q0, q10)**2
      t3 = 1 / F(Er, F0, s)
      
      res = t1 + t2 + t3
  end procedure aN

  module procedure bN
      real(c_double) :: t1, t2, t3
  
      ! Compute the result
      t1 = 1.0d0 / (2.0d0 * Nbar(Er, a, b, eps) * F(Er, F0, s))
      t2 = eps**2 / (2.0d0 * sigq(Eq, q0, q10)**2)
      t3 = V**2 / (2.0d0 * (sigp(Ep, eps, V, p0, p10) * 1.0d3)**2)
  
      res = t1 + t2 + t3
  end procedure bN

  module procedure cN1
      real(c_double) :: t1, t2, t3
  
      ! Compute the result
      t1 = -((Ep - Er)**2) / (2.0d0 * (sigp(Ep, eps, V, p0, p10)**2))
      t2 = -(Eq**2) / (2.0d0 * (sigq(Eq, q0, q10)**2))
      t3 = -Nbar(Er, a, b, eps) / (2.0d0 * F(Er, F0, s))
  
      res = t1 + t2 + t3
  end procedure cN1

  module procedure cN2
      real(c_double) :: sigp_val, sigq_val, Nbar_val, F_val
  
      ! Compute intermediate values
      sigp_val = sigp(Ep, eps, V, p0, p10)
      sigq_val = sigq(Eq, q0, q10)
      Nbar_val = Nbar(Er, a, b, eps)
      F_val = F(Er, F0, s)
  
      ! Compute the result
      res = 1.0d0 / (2.0d0 * pi * sqrt(2.0d0 * pi) * sigp_val * sigq_val * sqrt(Nbar_val * F_val))
  end procedure cN2

  module procedure PpqExp
      real(c_double) :: cN1_val, aN_val, bN_val
  
      ! Compute the components of the exponent
      cN1_val = cN1(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
      aN_val = aN(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10)
      bN_val = bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
      ! Compute the result
      res = cN1_val + (aN_val**2) / (4.0d0 * bN_val)
  end procedure PpqExp

  module procedure PpqPreExp
      real(c_double) :: exponent, cN2_val, bN_val
  
      ! Compute the exponent
      exponent = PpqExp(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
     
      ! first we test that the exponent won't cause a machine underflow
      if (exponent < -700.0d0) then
        res = 0.0d0
      else
        ! Compute the result
        ! if the exponent is not too small
        cN2_val = cN2(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
        bN_val = bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
        res = 0.5d0 * sqrt(pi) * cN2_val * exp(exponent) * (1.0d0 / sqrt(bN_val))
      end if
  end procedure PpqPreExp

  module procedure PpqFullN
      real(c_double) :: prefac, erffac, aN_val, bN_val
  
      ! Compute the prefactor
      prefac = PpqPreExp(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
      if (prefac > 0) then
        ! Compute aN and bN values
        ! if it's worth it
        aN_val = aN(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10)
        bN_val = bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
        ! Compute the error function factor
        erffac = erf(aN_val / (2.0d0 * sqrt(bN_val))) + 1.0d0
  
        ! Compute the result
        res = prefac * erffac * PErN(Er)
      else
        res = 0
      end if
  end procedure PpqFullN

  module procedure PpqFullG
      real(c_double) :: a, b, prefac, erffac, aN_val, bN_val
  
      ! the yield is always 1.0 for gamma (electron recoil) interactions
      a = 1.0
      b = 0.0
  
      ! Compute the prefactor
      prefac = PpqPreExp(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
      if (prefac > 0) then
        ! Compute aN and bN values
        ! if it's worth it
        aN_val = aN(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10)
        bN_val = bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
  
        ! Compute the error function factor
        erffac = erf(aN_val / (2.0d0 * sqrt(bN_val))) + 1.0d0
  
        ! Compute the result
        res = prefac * erffac * PErG(Er)
      else
        res = 0
      end if
  end procedure PpqFullG
end submodule PpqFort_s
