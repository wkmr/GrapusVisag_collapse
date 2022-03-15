subroutine nbody_deallocate_arrays
! Code deallocates nbody arrays ready for the next disc model

use embryodata
implicit none

if(allocated(pos)) deallocate(pos)
if(allocated(vel)) deallocate(vel)
if(allocated(acc)) deallocate(acc)
if(allocated(newpos)) deallocate(newpos)
if(allocated(newvel)) deallocate(newvel)
if(allocated(mass)) deallocate(mass)
if(allocated(angmom)) deallocate(angmom)
if(allocated(angmag)) deallocate(angmag)
if(allocated(ekin)) deallocate(ekin)
if(allocated(epot)) deallocate(epot)
if(allocated(etot)) deallocate(etot)

end subroutine nbody_deallocate_arrays
