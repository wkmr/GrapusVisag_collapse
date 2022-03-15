subroutine midplane_properties(t)
use unitdata, only:runmode
use stardata
implicit none

real :: t  

!real, dimension(nrannuli) :: tauplus

! Subroutine decides how to calculate midplane properties
! Given upper layer optical depth (tauplus)
! Depending on user selection

if(runmode=='g1') then
   call midplane_properties_grav_fixedalpha_1(t)
else
   call midplane_properties_grav_fixedalpha_2(t)
endif

end subroutine midplane_properties
