subroutine evolve_radius(j,t)
!****************************************************************
! Evolve radius of embryo j at time t
! Currently assumes a polytropic collapse (cf protostar formation)
!****************************************************************

use eosdata,only: pi
use embryodata

implicit none

integer, intent(in) :: j
real, intent(in) :: t
real :: rchoose


IF(embryo(j)%itidal==0) rchoose = embryo(j)%R0
IF(embryo(j)%itidal==1) rchoose = embryo(j)%R 

rchoose = embryo(j)%R0
!rchoose = embryo(j)%r

!embryo(j)%R = rchoose/(1.0 + 2.0*(t-embryo(j)%t_0_cool)/embryo(j)%t_cool0)**0.5
embryo(j)%R = rchoose/(1.0 + 2.0*t/embryo(j)%t_cool0)**0.5
!embryo(j)%R = rchoose/(1.0 + 2.0*dt/embryo(j)%t_cool0)**0.5
if ((embryo(j)%R .lt. rjup) .and. (embryo(j)%m .gt. 0.5d0*mjup)) then 
  embryo(j)%r = rjup
  embryo(j)%R0 = embryo(j)%R*(1.0d0 + 2.0d0*t/embryo(j)%t_cool0)**0.5d0
endif

embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)           

!print*, j, rchoose/1.5d13, embryo(j)%r0/1.5d13, embryo(j)%r/1.5d13, embryo(j)%t_cool0

end subroutine evolve_radius
