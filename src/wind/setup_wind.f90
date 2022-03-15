SUBROUTINE setup_wind
!
! Routine takes inputs and sets up disc winds
! X Ray winds from Owen et al (2010)
! FUV/EUV winds from Pont et al, Alexander et al
!

use stardata
use winddata
use unitdata

implicit none

integer :: i
real :: xwind,sigdotxray

If(allocated(sigdot_wind)) deallocate(sigdot_wind)
allocate(sigdot_wind(nrannuli))

sigdot_wind(:) = 0.0

!print*, 'Setting up wind parameters', mstar/solarmass

mdot_wind = 0.0d0
windnorm = 0.0d0

If (mstar .gt. 0.25d0*solarmass) Then
  do i = isr, ier
    xwind = 0.85*(rz(i)/1.496d13)*(mstar/1.989d33)**(-1.0d0)
    call xraywind(sigdotxray,xwind)
    mdot_wind = mdot_wind + sigdotxray*2.0d0*3.14159*rz(i)/drzm1(i)
  enddo

!  Lx = (Lx - Lx*ran2(iseed)/4.0)
  If (Lx .lt. 1.0e29) Then
    Lx = 1.0e29
  EndIf

!  print*, 'Xray luminosity set to: ', Lx,'  erg/s'

  windnorm = 6.25d-9*(mstar/solarmass)**(-0.068)*(Lx/1.0d30)**1.14/mdot_wind

Else
  windnorm = 0.0d0
EndIf

return
END SUBROUTINE setup_wind
