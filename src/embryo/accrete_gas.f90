subroutine accrete_gas(j,t)
!****************************************
! Accrete gas from disc onto embryo j
! depends on if a gap is opened
! (Future development)
!****************************************

use stardata
use embryodata
use eosdata
use unitdata

implicit none

real,intent(in) :: t

real :: r_hill, r_acc, r_tot
real :: r_feed_inf, r_feed_sup
real :: M_feed
real :: sig_cross, v_rel
real :: Mdot_env_3D, Mdot_env_2D, accr_mass
real :: Mdot_env_max

integer :: i, k, l
integer, intent(in) :: j
integer :: numpts

sigdot_accr(:) = 0.0d0

!if(embryo(j)%migtype==1) then

call find_planets_in_disc

i = embryo(j)%icurrent

If (mstar .gt. 0.0d0) Then
  r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.3333d0
Else 
  r_hill = 0.0d0
EndIf  

If (cs_d(i) .gt. 0.0d0) then
  r_acc = G*embryo(j)%m/(cs_d(i))**2.0d0
else
  r_acc = 0.0d0
endif

If (r_hill .gt. 0.0d0) Then
  If (r_acc .gt. 0.0d0) then
    r_tot = 1.0/(1.0*r_acc) + 1.0/(0.25d0*r_hill)
    r_tot = 1.0/r_tot
  else
    r_tot = 1.0d0/(0.25d0*r_hill)
    r_tot = 1.0d0/r_tot
  endif
Else
  r_tot = 0.0d0
EndIf  

if (omega_d(i)*r_tot .lt. cs_d(i)) then
  v_rel = cs_d(i)
else
  v_rel = omega_d(i)*r_tot
endif

If (mstar .gt. 0.0d0) Then
  r_feed_inf = embryo(j)%a*(1.0d0-0.0d0)
  r_feed_inf = r_feed_inf*(1.0d0 - 0.5*(embryo(j)%m/(3.0d0*mstar))**0.5d0)

  r_feed_sup = embryo(j)%a*(1.0d0-0.0d0)
  r_feed_sup = r_feed_sup*(1.0d0 + 0.5*(embryo(j)%m/(3.0d0*mstar))**0.5d0)
Else
  r_feed_inf = embryo(j)%a
  r_feed_sup = embryo(j)%a
EndIf  

M_feed = 2.0d0*pi*embryo(j)%a*sigma_d(i)*(r_feed_sup - r_feed_inf)

Mdot_env_2D = 2.0d0*r_tot*sigma_d(i)*v_rel

Mdot_env_3D = 0.0d0
if (H_d(i) .gt. 0.0d0) then
  Mdot_env_3D = sigma_d(i)/H_d(i)*pi*r_tot**2.0d0*v_rel
endif

Mdot_env_max = Mdot_env_3D
If (Mdot_env_max .gt. Mdot_env_2D) Then
  Mdot_env_max = Mdot_env_2D
EndIf

mdotvisc = 0.0d0
numpts = 0
do k = isr, ier 
  if (sigma_d(k) .gt. 0.0d0) then
    mdotvisc = mdotvisc + 3.0*pi*alpha_d(k)*cs_d(k)*cs_d(k)*sigma_d(k)/omega_d(k)
    numpts = numpts + 1
  endif
enddo

mdotvisc = mdotvisc/numpts

if (Mdot_env_max*dt .gt. M_feed) then
  Mdot_env_max = M_feed/dt
endif

Mdot_env_max = Mdot_env_max*0.25d0

If (accr_on_disc == 'y') then
  Mdot_env_max = Mdot_env_max + 0.01*mdotvisc
EndIf

accr_mass = 0.0d0
do l = isr, ier
  sigdot_accr(l) = 0.0d0
  if (rz(l) .ge. r_feed_inf) then
    if (rz(l) .le. r_feed_sup) then
      sigdot_accr(l) = Mdot_env_max/(pi*(r_feed_sup**2.0d0 - r_feed_inf**2.0d0))
    endif
  endif

  accr_mass = accr_mass + sigdot_accr(l)*pi*(rz(l+1)**2.0d0-rz(l)**2.0d0)*dt
enddo
if (accr_mass .gt. 0.0d0) then
  sigdot_accr(:) = sigdot_accr(:)*Mdot_env_max*dt/accr_mass
endif

embryo(j)%M = embryo(j)%M + Mdot_env_max*dt 

embryo(j)%R = (embryo(j)%m/(4.0d0*theta_grad*embryo(j)%rhoc*pi))**(1.0d0/3.0d0)

if (embryo(j)%idiss == 1)  then
  IF(embryo(j)%M/Mjup < 70.0 ) THEN
      embryo(j)%R = -2.8e-3*embryo(j)%M/Mjup + 1.0
   ELSE
      embryo(j)%R = 0.01* embryo(j)%M/Mjup + 0.2
   ENDIF
   embryo(j)%R = embryo(j)%R*rjup
endif

embryo(j)%R0 = embryo(j)%R*(1.0 + 2.0*(t-embryo(j)%t_form)/embryo(j)%t_cool0)**0.5

end subroutine accrete_gas

