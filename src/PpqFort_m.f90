module PpqFort_m
  !! TODO: Add comment for FORD describing the module's purpose
  use iso_c_binding, only : c_double, c_int
  implicit none

  interface

    pure module subroutine PpqN_vector(Ep_arr, Eq_arr, n, k, Z, F0, eps, V, p0, p10, q0, q10, res_arr) bind(c, name="PpqN_vector")
      !! TODO: Add comment for FORD describing the function's purpose
      integer(c_int), value :: n
      real(c_double), intent(in) :: Ep_arr(n), Eq_arr(n)
      real(c_double), value :: k, Z, F0, eps, V, p0, p10, q0, q10
      real(c_double), intent(out) :: res_arr(n)
    end subroutine PpqN_vector

    pure module subroutine PpqG_vector(Ep_arr, Eq_arr, n, F0, eps, V, p0, p10, q0, q10, res_arr) bind(c, name="PpqG_vector")
      !! TODO: Add comment for FORD describing the function's purpose
      integer(c_int), value :: n
      real(c_double), intent(in) :: Ep_arr(n), Eq_arr(n)
      real(c_double), value :: F0, eps, V, p0, p10, q0, q10
      real(c_double), intent(out) :: res_arr(n)
    end subroutine PpqG_vector

    pure module function PpqN(Ep, Eq, k, Z, F0, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqN")
      !! TODO: Add comment for FORD describing the function's purpose
      real(c_double), value :: Ep, Eq, k, Z, F0, eps, V, p0, p10, q0, q10
      real(c_double) :: res
    end function PpqN

    pure module function PpqG(Ep, Eq, F0, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqG")
      !! TODO: Add comment for FORD describing the function's purpose
      real(c_double), value :: Ep, Eq, F0, eps, V, p0, p10, q0, q10
      real(c_double) :: res
    end function PpqG

   pure module function Y(Er, k, Z) result(res) bind(c, name="Y")
       !! the ionization yield, using the Lindhard model
       !! Y(Er) = k*g(eps_L) / (1 + k*g(eps_L)), with
       !! eps_L = 11.5*Er*Z**(-7/3) and g(eps_L) = 3*eps_L**0.15 + 0.7*eps_L**0.29 + eps_L
       real(c_double), value :: Er, k, Z
       real(c_double) :: res
    end function Y

   pure module function Nbar(Er, k, Z, eps) result(res) bind(c, name="Nbar")
     !! average number of electron-hole pairs for a given Er
     real(c_double), value :: Er, k, Z, eps
     real(c_double) :: res
   end function Nbar

   pure module function sigp(Ep, eps, V, p0, p10) result(res) bind(c, name="sigp")
      !! phonon sensor resolution
      real(c_double), value :: Ep, eps, V, p0, p10
      real(c_double) :: res
   end function sigp

   pure module function sigq(Eq, q0, q10) result(res) bind(c, name="sigq")
   !! charge sensor resolution
      real(c_double), value :: Eq, q0, q10
      real(c_double) :: res
   end function sigq

  pure module function PErN(Er) result(res) bind(c, name="PErN")
    !! the Er (energy) distribution of neutrons
    real(c_double), value :: Er
    real(c_double) :: res
  end function PErN

  pure module function PErG(Er) result(res) bind(c, name="PErG")
    !! the Er (energy) distribution of gamma events/electron recoils
    real(c_double), value :: Er
    real(c_double) :: res
  end function PErG

  pure module function PpqFullN(Er, Ep, Eq, k, Z, F0, eps, V, p0, p10, q0, q10) result(res) bind(c, name="PpqFullN")
    ! Input arguments
    real(c_double), value :: Er, Ep, Eq, k, Z, F0, eps, V, p0, p10, q0, q10
    real(c_double) :: res
  end function PpqFullN

  pure module subroutine PpqFort_version(major, minor, patch) bind(c, name="PpqFort_version")
    !! Report the package version (see fpm.toml's version field, which
    !! must be kept in sync by hand)
    integer(c_int), intent(out) :: major, minor, patch
  end subroutine PpqFort_version

  end interface

  ! Define pi
  real(c_double), parameter :: pi = 3.14159265358979323846_c_double
    ! TODO: move the pi definition to just above the submodule's "contains" statement

  ! Package version (semver).  Keep in sync with fpm.toml's version field.
  integer(c_int), parameter :: version_major = 1
  integer(c_int), parameter :: version_minor = 0
  integer(c_int), parameter :: version_patch = 0

end module PpqFort_m
