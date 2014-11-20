
module planck

use shr_kind_mod,  only: r8 => shr_kind_r8

use shr_const_mod, only: k_b => shr_const_boltz, & ! Boltzman's constant (J/K/molecule)
                         pi  => shr_const_pi       ! 3.14...

implicit none

real(r8), parameter :: &
    h = 6.6260755e-34_r8, & ! Planck's constant (J.s)
    c = 2.99792458e8_r8     ! Speed of light in a vacuum (m/s)
integer, parameter :: nquad = 50 ! number of quadrature values

contains

subroutine integrated_planck(wavenumber_start, wavenumber_end, t, flux)
    real(r8), intent(in)  :: wavenumber_start
    real(r8), intent(in)  :: wavenumber_end ! cm^-1
    real(r8), intent(in)  :: t                                ! temperature in K
    real(r8), intent(out) :: flux                             ! W/(m^2*sr*cm^-1)
    real(r8) :: delta  ! delta(wavenumber)
    real(r8) :: cflux  ! cumulative integral value of flux during integration
    real(r8) :: quad_midpoint
    integer :: iquad

    delta = (wavenumber_end - wavenumber_start) / nquad
    cflux = 0._r8
    do iquad = 1, nquad
      quad_midpoint = wavenumber_start + delta * (iquad - 0.5_r8)
      cflux = cflux + planckf(quad_midpoint, t)
    enddo

    flux = cflux * delta * pi * 1.e8_r8

    return
end subroutine integrated_planck


function planckf(wavenumber, t)
    real(r8), intent(in)  :: wavenumber ! cm^-1
    real(r8), intent(in)  :: t          ! temperature in K
    real(r8) planckf                    ! W/(m^2*sr*cm^-1)
    real(r8) exponenta

    exponenta = 100._r8*h*c*wavenumber/(k_b*t)
    planckf = 2._r8 * h * c**2 * wavenumber**3 / &
        (exp(exponenta) - 1._r8)

    return
end function planckf

end  module
