subroutine find_planets_in_disc
!
! Finds grid cell planets are in
!

use stardata, only: rz,isr,nrannuli
use embryodata
use planetdata
use unitdata
implicit none

integer :: icheck,iplanet

icheck = 0
do iplanet=1,nplanet

      ! Check from vicinity of last known location
      iplanetrad(iplanet) = iplanetrad(iplanet)-int(0.1*nrannuli)

      ! If at start or near inner boundary, use it as starting point
      if(iplanetrad(iplanet)<isr) iplanetrad(iplanet) = isr

   do while((rz(iplanetrad(iplanet))<ap(iplanet)).and.(iplanetrad(iplanet)<(nrannuli-1)))
      iplanetrad(iplanet) = iplanetrad(iplanet)+1
      if (iplanetrad(iplanet)>(nrannuli-1)) iplanetrad(iplanet)=nrannuli-1
   enddo

   embryo(iplanet)%icurrent = iplanetrad(iplanet)

enddo


end subroutine find_planets_in_disc
