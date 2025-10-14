program test_ppq
  use iso_c_binding
  use PpqFort
  implicit none

  integer, parameter :: n = 3
  real(c_double) :: a, b, F0, s, eps, V, p0, p10, q0, q10
  real(c_double) :: Eq_arr(n), Ep_arr(n), resN(n), resG(n)
  integer(c_int) :: ni

  ! Parameters
  a   = 0.16d0
  b   = 0.18d0
  F0  = 0.122d0
  s   = 0.0d0
  eps = 3.0d-3
  V   = 3.0d0
  p0  = 0.06421907d0
  p10 = 0.48998486d0
  q0  = 0.06421907d0
  q10 = 0.48998486d0

  ! Arrays
  Eq_arr = (/100.0d0, 100.0d0, 100.0d0/)
  Ep_arr = (/347.0d0, 346.0d0, 348.0d0/)

  ni = n

  ! Call PpqN_vector
  call PpqN_vector(Ep_arr, Eq_arr, ni, a, b, F0, s, eps, V, p0, p10, q0, q10, resN)

  ! Call PpqG_vector
  call PpqG_vector(Ep_arr, Eq_arr, ni, F0, s, eps, V, p0, p10, q0, q10, resG)

  ! Print results
  print *, "Results from PpqN_vector:"
  print '(3ES20.10E3)', resN

  print *, "Results from PpqG_vector:"
  print '(3ES20.10E3)', resG

end program test_ppq


