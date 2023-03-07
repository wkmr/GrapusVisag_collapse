subroutine compute_wind
!
! Compute the surface density loss due to winds
! Interpolates between X Ray Winds (Owen et al 2010)
! and FUV/EUV winds (Pont et al, Alexander et al)
!

use winddata
use stardata
use unitdata

implicit none

integer :: i
real*8 sigdot_diff, sigdot_dir,rdir_min
real :: xwind, sigdotxray

sigdot_wind(:) = 0.0

rdir_min = 0.1d0*AU

If (mstar .gt. 0.25d0*solarmass) Then
  do i= isr,ier
   Call wind_profiles(sigdot_diff,sigdot_dir,rz(i), mstar, rdir_min)

   If (sigdot_diff.lt.0.0d0) sigdot_diff = 0.0d0
   If (sigdot_dir.lt.0.0d0) sigdot_dir = 0.0d0

   xwind = 0.85*(rz(i)/AU)*(mstar/solarmass)**(-1.0d0)

   call xraywind(sigdotxray,xwind)
   sigdotxray = sigdotxray*windnorm

   sigdotxray = sigdotxray*msolyr

   if(rz(i)/AU > rwind_xray) then
     sigdot_wind(i) = rwind_xray*(sigdot_diff + sigdot_dir)
!     sigdot_wind(i) = 0.0d0
   else
     sigdot_wind(i) = sigdotxray
   endif
  
!   sigdot_wind(i) = 10.0d0*(sigdot_diff + sigdot_dir)

  enddo
Else
  sigdot_wind(:) = 0.0d0
EndIf

end subroutine compute_wind
