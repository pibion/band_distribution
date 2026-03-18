program profile_driver
  use iso_c_binding, only : c_double, c_int
  use PpqFort_m, only : PpqN_vector, PpqG_vector
  implicit none

  ! Parameters from test suite
  real(c_double), parameter :: a   = 0.16d0
  real(c_double), parameter :: b   = 0.18d0
  real(c_double), parameter :: F0  = 0.122d0
  real(c_double), parameter :: s   = 0.0d0
  real(c_double), parameter :: eps = 3.0d-3
  real(c_double), parameter :: V   = 3.0d0
  real(c_double), parameter :: p0  = 0.06421907d0
  real(c_double), parameter :: p10 = 0.48998486d0
  real(c_double), parameter :: q0  = 0.06421907d0
  real(c_double), parameter :: q10 = 0.48998486d0

  integer, parameter :: npts = 10
  integer(c_int), parameter :: n = int(npts, c_int)

  real(c_double) :: Ep_arr(npts), Eq_arr(npts)
  real(c_double) :: resN(npts), resG(npts)
  real(c_double) :: dummy
  integer :: i, iter
  integer, parameter :: n_iter = 1

  ! Build a grid of Ep and Eq values spanning the physically relevant range
  do i = 1, npts
    Ep_arr(i) = 300.0d0 + (100.0d0 * (i - 1)) / (npts - 1)  ! 300 to 400
    Eq_arr(i) =  50.0d0 + (100.0d0 * (i - 1)) / (npts - 1)  !  50 to 150
  end do

  ! Repeat computation to get a long enough run for TAU sampling
  do iter = 1, n_iter
    call PpqN_vector(Ep_arr, Eq_arr, n, a, b, F0, s, eps, V, p0, p10, q0, q10, resN)
    call PpqG_vector(Ep_arr, Eq_arr, n,       F0, s, eps, V, p0, p10, q0, q10, resG)
  end do

  ! Use results to prevent the compiler from optimizing the loop away
  dummy = sum(resN) + sum(resG)
  print *, "sum(resN) + sum(resG) =", dummy

end program profile_driver
