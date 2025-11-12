submodule(PpqFort_m) PpqFort_s
  implicit none

contains

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

end submodule PpqFort_s
