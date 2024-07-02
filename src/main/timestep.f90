subroutine timestep(t)
!
!
! Determine maximum safe timestep for the simulation
! Must account for viscous evolution and planetary torques
! Timestep cannot exceed predefined maximum
!
!

use stardata
use eosdata
use embryodata
use planetdata

implicit none

integer :: i, iplanet
real(kind=8) :: dtmin,C0,C1,dr2
real :: t
real :: dtvisc, dtmin_visc, dttorq, dtmin_torque, dtmin_planet

dtmin_visc = 1.0e30
dtmin_torque = dtmin_visc
dtmin = dtmin_visc

If (t .lt. 1.0d8) Then
   C0 = 0.01d0
Else If (t .lt. 1.0d9) Then
   C0 = 0.025d0
Else
   C0 = 0.25d0
EndIf
C1 = 0.5d0

!C0 = C0/10.0

iplanet = 1

If (nembryo .gt. 0) Then
  adot(iplanet) = 0.0d0
EndIf  

do i = isr, ier

   ! Compute minimum viscous timestep
   dr  = (rf(i+1)-rf(i))
   dr2 = dr**2
   if ((nu_tc(i) .gt. 0.0d0) .and. (sigma_d(i) .gt. 0.0d0)) then
     dtvisc  = C0 * dr2 / (6.0d0 * nu_tc(i))
   else
     dtvisc = 1.0d30
   endif
   dtmin_visc = min(dtmin_visc,dtvisc)
!   dtmin = min(dtmin,0.1*tcool(i))

!   print*, i, t, dtvisc, dtmin_visc, dr2, nu_tc(i), alpha_d(i), cs_d(i), sigma_d(i)

   planetchoice = 'y'

  ! If planets present, compute minimum timestep due to torques
   if ((planetchoice=='y') .and. (nembryo .gt. 0)) then
      if(torque_term(i)>1.0e-30) then
         dttorq = C0*dr*rz(i)*sigma_d(i)/torque_term(i)
      else
         dttorq = 1.0e30
      endif

      dtmin_torque = min(dtmin_torque, dttorq)   

      adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma_d(i)/drzm1(i)

   endif

enddo

adot(iplanet) = abs(adot(iplanet)*(ap(iplanet)*G*mstar)**0.5*(4.0d0*pi/mp(iplanet)))

dtmin_planet = au/adot(iplanet)
If ((iplanetrad(iplanet) .gt. isr) .and. (iplanetrad(iplanet) .lt. ier)) Then
  dtmin_planet =  (rf(iplanetrad(iplanet)+1)-rf(iplanetrad(iplanet)-1))/4.0/adot(iplanet)
EndIf  

!dt = dtmin_visc
!dt = min(dtmin, dtmin_visc, dtmin_torque) 
dt = min(dtmin_visc, dtmin_torque, dtmin_planet)
dt = min(dt,maxstep*yr) ! timestep can not exceed maxstep
    
!print*, t/yr,dt/yr, dtmin_visc, dtmin_torque
!IF(dt>tdump) then
!   print*, 'Reducing dt for rapid snapshotting'
!   dt = tdump/10.0
!endif

end subroutine timestep
