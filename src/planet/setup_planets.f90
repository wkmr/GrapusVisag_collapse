subroutine setup_planets
!
! Subroutine sets up planets to be added to disc
! Planet data read from separate input file
!
!

use stardata
use planetdata
use unitdata
use embryodata

implicit none

integer :: iplanet

! Open planet file

nplanet = nembryo 

!read(10,*) nplanet

print*, 'There are ',nplanet, 'planets'
nactive = nplanet

if(allocated(mp)) deallocate(mp)
if(allocated(ap)) deallocate(ap)
if(allocated(alive)) deallocate(alive)
if(allocated(iplanetrad)) deallocate(iplanetrad)
if(allocated(lambdaI)) deallocate(lambdaI)
if(allocated(lambdaII)) deallocate(lambdaII)
if(allocated(fII)) deallocate(fII)
if(allocated(adot)) deallocate(adot)
if(allocated(tmig)) deallocate(tmig)
if(allocated(tmigI)) deallocate(tmigI)
if(allocated(torquei)) deallocate(torquei)
if(allocated(torque_term)) deallocate(torque_term)
if(allocated(total_planet_torque)) deallocate(total_planet_torque)

allocate(mp(nplanet),ap(nplanet),alive(nplanet), iplanetrad(nplanet))
allocate(lambdaI(nplanet,nrannuli), lambdaII(nplanet,nrannuli))
allocate(fII(nplanet))
allocate(adot(nplanet),tmig(nplanet),tmigI(nplanet))
allocate(torquei(nplanet,nrannuli), torque_term(nrannuli), total_planet_torque(nrannuli))


alive(:) = 1
mp(:) = 0.0
ap(:) = 0.0
iplanetrad(:) = 0

lambdaII(:,:) = 0.0
lambdaI(:,:) = 0.0
fII(:) = 0.0
adot(:) = 0.0
tmig(:) = 0.0
tmigI(:) = 0.0
torquei(:,:) = 0.0
total_planet_torque(:) = 0.0

do iplanet=1,nplanet
!   read(10,*) mp(iplanet), ap(iplanet)
  mp(iplanet) = embryo(iplanet)%m/mjup
  ap(iplanet) = embryo(iplanet)%a/AU
enddo
 

! Convert to correct units
mp(:) = mp(:)*mjup
ap(:) = ap(:)*AU

call find_planets_in_disc

do iplanet=1,nplanet
      print*, 'Planet ', iplanet, 'initially located at cell ', iplanetrad(iplanet)
      print*, 'Radius: ', rz(iplanetrad(iplanet))/AU
      print*, 'Mass: ' , mp(iplanet)/mjup
enddo

end subroutine setup_planets
