module PpqFort
  use iso_c_binding
  implicit none

  ! Define pi
  real(c_double), parameter :: pi = acos(-1.0d0)

  contains
   !the ionization yield. 
   pure function Y(Er, a, b) result(res) bind(c, name="Y")
       real(c_double), value :: Er, a, b
       real(c_double) :: res
       
       ! Compute the ionization yield
       res = a * abs(Er)**b
    end function Y
  
    ! Model the Fano factor as a linear function
    pure function F(Er, F0, s) result(res) bind(c, name="F")
      real(c_double), value :: Er, F0, s
      real(c_double) :: res

      res = F0 + s * Er
   end function F

   ! average number of electron-hole pairs for a given Er
   pure function Nbar(Er, a, b, eps) result(res) bind(c, name="Nbar")
     real(c_double), value :: Er, a, b, eps
     real(c_double) :: res
   
     res = Y(Er, a, b) * Er / eps
   end function Nbar

   ! phonon sensor resolution
   pure function sigp(Ep, eps, V, p0, p10) result(res) bind(c, name="sigp")
      real(c_double), value :: Ep, eps, V, p0, p10
      real(c_double) :: res
      real(c_double) :: e2pre, e2e

      e2pre = (p10**2 - p0**2)
      e2e = (Ep / (10.0d0 * (1.0d0 + (V / eps / 1000.0d0))))**2
      res = sqrt(p0**2 + e2pre * e2e)
   end function sigp

   ! charge sensor resolution
   pure function sigq(Eq, q0, q10) result(res) bind(c, name="sigq")
      real(c_double), value :: Eq, q0, q10
      real(c_double) :: res
      real(c_double) :: e2pre, e2e

      e2pre = (q10**2 - q0**2)
      e2e = (Eq / 10.0d0)**2
      res = sqrt(q0**2 + e2pre * e2e)
   end function sigq

    ! the Er (energy) distribution of neutrons
    pure function PErN(Er) result(res) bind(c, name="PErN")
    real(c_double), value :: Er
    real(c_double) PNa, PNb, PNd
    real(c_double) :: res

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
  end function PErN

  ! the Er (energy) distribution of gamma events/electron recoils
  pure function PErG(Er) result(res) bind(c, name="PErG")
    real(c_double), value :: Er
    real(c_double) :: PGa, PGb, PGd
    real(c_double) :: res

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
  end function PErG

  pure function aN(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="aN")
    real(c_double), value :: Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10
    real(c_double) :: t1, t2, t3, res
    
    t1 = (V/1.0d3)*(Ep - Er) / sigp(Ep, eps, V, p0, p10)**2
    t2 = eps*Eq / sigq(Eq, q0, q10)**2
    t3 = 1 / F(Er, F0, s)
    
    res = t1 + t2 + t3
  end function aN

  pure function bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="bN")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: t1, t2, t3

    ! Output
    real(c_double) :: res

    ! Compute the result
    t1 = 1.0d0 / (2.0d0 * Nbar(Er, a, b, eps) * F(Er, F0, s))
    t2 = eps**2 / (2.0d0 * sigq(Eq, q0, q10)**2)
    t3 = V**2 / (2.0d0 * (sigp(Ep, eps, V, p0, p10) * 1.0d3)**2)

    res = t1 + t2 + t3
  end function bN

  pure function cN1(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="cN1")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: t1, t2, t3

    ! Output
    real(c_double) :: res

    ! Compute the result
    t1 = -((Ep - Er)**2) / (2.0d0 * (sigp(Ep, eps, V, p0, p10)**2))
    t2 = -(Eq**2) / (2.0d0 * (sigq(Eq, q0, q10)**2))
    t3 = -Nbar(Er, a, b, eps) / (2.0d0 * F(Er, F0, s))

    res = t1 + t2 + t3
  end function cN1

  pure function cN2(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="cN2")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: sigp_val, sigq_val, Nbar_val, F_val

    ! Output
    real(c_double) :: res

    ! Compute intermediate values
    sigp_val = sigp(Ep, eps, V, p0, p10)
    sigq_val = sigq(Eq, q0, q10)
    Nbar_val = Nbar(Er, a, b, eps)
    F_val = F(Er, F0, s)

    ! Compute the result
    res = 1.0d0 / (2.0d0 * pi * sqrt(2.0d0 * pi) * sigp_val * sigq_val * sqrt(Nbar_val * F_val))
  end function cN2

  pure function PpqExp(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqExp")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: cN1_val, aN_val, bN_val

    ! Output
    real(c_double) :: res

    ! Compute the components of the exponent
    cN1_val = cN1(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
    aN_val = aN(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10)
    bN_val = bN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)

    ! Compute the result
    res = cN1_val + (aN_val**2) / (4.0d0 * bN_val)
  end function PpqExp

  pure function PpqPreExp(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqPreExp")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: exponent, cN2_val, bN_val

    ! Output
    real(c_double) :: res

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
  end function PpqPreExp

  pure function PpqFullN(Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqFullN")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: prefac, erffac, aN_val, bN_val

    ! Output
    real(c_double) :: res

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
  end function PpqFullN

  pure function PpqFullG(Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqFullG")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, F0, s, eps, V, p0, p10, q0, q10

    ! Local variables
    real(c_double) :: a, b, prefac, erffac, aN_val, bN_val

    ! Output
    real(c_double) :: res

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
  end function PpqFullG

  pure function PpqN(Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqN")
      ! Inputs
      real(c_double), value :: Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10
  
      ! Locals
      real(c_double), allocatable :: er_arr(:)
      real(c_double) :: minx, maxx, resolution, res, integral
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
      ! gfortran may not support this in our Expanse version
      ! we want the latest llvm-flang (v. 20) compiler
      ! do concurrent(i = 1:npts) default (none)
      do concurrent(i = 1:npts)
          integral = integral + PpqFullN(er_arr(i), Ep, Eq, a, b, F0, s, eps, V, p0, p10, q0, q10)
      end do
  
      res = integral * resolution
  end function PpqN

  pure function PpqG(Ep, Eq, F0, s, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqG")
      ! Inputs
      real(c_double), value :: Ep, Eq, F0, s, eps, V, p0, p10, q0, q10
  
      ! Locals
      real(c_double), allocatable :: er_arr(:)
      real(c_double) :: minx, maxx, resolution, res, integral
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
  end function PpqG

  subroutine PpqN_vector(Ep_arr, Eq_arr, n, a, b, F0, s, eps, V, p0, p10, q0, q10, res_arr) bind(c, name="PpqN_vector")
    ! Inputs
      integer(c_int), value :: n
      real(c_double), intent(in) :: Ep_arr(n), Eq_arr(n)
      real(c_double), value :: a, b, F0, s, eps, V, p0, p10, q0, q10
  
      ! Output
      real(c_double), intent(out) :: res_arr(n)
  
      ! Locals
      integer :: i
  
      do concurrent(i = 1:n)
          res_arr(i) = PpqN(Ep_arr(i), Eq_arr(i), a, b, F0, s, eps, V, p0, p10, q0, q10)
      end do
  
  end subroutine PpqN_vector

    subroutine PpqG_vector(Ep_arr, Eq_arr, n, F0, s, eps, V, p0, p10, q0, q10, res_arr) bind(c, name="PpqG_vector")
    ! Inputs
      integer(c_int), value :: n
      real(c_double), intent(in) :: Ep_arr(n), Eq_arr(n)
      real(c_double), value :: F0, s, eps, V, p0, p10, q0, q10
  
      ! Output
      real(c_double), intent(out) :: res_arr(n)
  
      ! Locals
      integer :: i
  
      do concurrent(i = 1:n)
          res_arr(i) = PpqG(Ep_arr(i), Eq_arr(i), F0, s, eps, V, p0, p10, q0, q10)
      end do
  
  end subroutine PpqG_vector

end module

