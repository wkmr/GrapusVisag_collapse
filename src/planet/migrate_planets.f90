subroutine migrate_planets
!
! Subroutine uses the torque exerted by each planet to find the
! migration timescale for each planet
!

use planetdata
use stardata
use unitdata
use embryodata

implicit none

integer :: iplanet,i

!
! Assume that planet torque contributes to its migration only
! (Detailed torque balance between disc and individual bodies)
!

do iplanet=1,nplanet

   ! Integrate torque*sigma over the entire disc
   adot(iplanet) = 0.0
    
   do i=isr,ier       
!      adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma_d(i)/drzm1(i)
      adot(iplanet) = adot(iplanet) + torquei(iplanet,i)*sigma_d(i)/drzm1(i)
   enddo

    ! Multiply by appropriate factors to get adot    
!   adot(iplanet) = -adot(iplanet)*4.0*pi* &
!        omega_d(iplanetrad(iplanet))*ap(iplanet)*ap(iplanet)/mp(iplanet)
    adot(iplanet) = -adot(iplanet)*(ap(iplanet)*G*mstar)**0.5*(4.0d0*pi/mp(iplanet))
   !adot(iplanet) = -adot(iplanet)*4.0*pi*ap(iplanet)/mp(iplanet)

!    print*, adot(iplanet)

    ! get numerically determined migration timescale

    if(abs(adot(iplanet))>1.0e-30) then
       tmig(iplanet) = ap(iplanet)/(-adot(iplanet))
    else
       tmig(iplanet) = 0.0
    endif
           
    ! Move planets
    ap(iplanet) = ap(iplanet) + adot(iplanet)*dt

    if (ap(iplanet) .gt. rz(nrannuli)) ap(iplanet) = rz(nrannuli)
enddo

end subroutine migrate_planets
