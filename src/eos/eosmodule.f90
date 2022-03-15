      module eosdata

!-----------------------------------------------------------------------
! Data module for saving equation of state values
! PJC 22/05/2008
! DHF 09/09/2009
!-----------------------------------------------------------------------
use unitdata         

implicit none
save

! Units
integer :: nrhopoints,nUpoints
real(kind=8),parameter :: udist = 1.496d13
real(kind=8),parameter :: umass = 1.989d33
real(kind=8),parameter :: utime = sqrt((udist**3)/(G*umass))
real(kind=8),parameter :: uergg = udist*udist/(utime*utime)

! Arrays
real,dimension(5) :: gammamuT
real,allocatable,dimension(:,:,:) :: eostable
real,allocatable,dimension(:,:) :: cstab

!$OMP threadprivate(gammamuT)

end module eosdata
