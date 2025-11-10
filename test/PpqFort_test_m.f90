module PpqFort_test_m
  use iso_c_binding, only : c_double, c_int
  use julienne_m, only : &
     test_t, test_description_t, test_diagnosis_t, test_result_t &
    ,operator(.approximates.), operator(.within.), operator(.all.), operator(//)
  use PpqFort, only :  PpqG_vector, PpqN_vector
  implicit none

  type, extends(test_t) :: PpqFort_test_t
  contains
    procedure, nopass :: subject
    procedure, nopass :: results
  end type

  ! Parameters
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

  ! Arrays
 real(c_double), parameter :: Eq_arr(*) = [100.0d0, 100.0d0, 100.0d0]
 real(c_double), parameter :: Ep_arr(*) = [347.0d0, 346.0d0, 348.0d0]
 integer(c_int), parameter :: ni = int(size(Eq_arr), c_int)

contains

  pure function subject() result(test_subject)
    character(len=:), allocatable :: test_subject
    test_subject = 'The PpqFort module'
  end function

  function results() result(test_results)
    type(test_result_t), allocatable :: test_results(:)
    type(PpqFort_test_t) :: PpqFort_test

    test_results = PpqFort_test%run([ & 
       test_description_t('computing the resN array', check_PpqG_vector) &
      ,test_description_t('computing the resG array', check_PpqN_vector) &
    ])
  end function

  function check_PpqN_vector() result(test_diagnosis)
    type(test_diagnosis_t) test_diagnosis
    real(c_double), parameter :: expected_resN(*) = [real(c_double) :: 2.7532559194D-008, 2.9389448480D-008,  2.5705989565D-008]
    real(c_double), parameter :: tolerance(*) = expected_resN/100
    real(c_double) resN(size(Eq_arr))

    call PpqN_vector(Ep_arr, Eq_arr, ni, a, b, F0, s, eps, V, p0, p10, q0, q10, resN)
    test_diagnosis =  .all. ( resN .approximates. expected_resN .within. tolerance)  // ' (PpqN)'
  end function

  function check_PpqG_vector() result(test_diagnosis)
    type(test_diagnosis_t) test_diagnosis
    real(c_double), parameter :: expected_resG(*) = [real(c_double) :: 1.2591349782D-033, 2.6066426223D-033, 6.0702183456D-034]
    real(c_double), parameter :: tolerance(*) = expected_resG/100
    real(c_double) resG(size(Eq_arr))

    call PpqG_vector(Ep_arr, Eq_arr, ni, F0, s, eps, V, p0, p10, q0, q10, resG)
    test_diagnosis =  .all. ( resG .approximates. expected_resG .within. tolerance)  // ' (PpqN)'
  end function

end module PpqFort_test_m
