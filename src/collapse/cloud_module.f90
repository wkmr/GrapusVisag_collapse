module clouddata
!---------------------------------------------------------
! This module contains all the properties of the cloud
!-------------------------------------------------------
  
 real(kind=8) :: Mcloud_0, Mcloud_1, Rcloud_0, Rcloud_1
 real(kind=8) :: f_cloud_0, f_cloud_1, t_frag
 real(kind=8) :: rhocloud, beta_cloud, amax_cloud, tff, mdot_cloud, mdot 
 real(kind=8) :: omega_cloud, dmstar_cloud
 real(kind=8) :: dMcloud, dMtot, MRouter, Router, Mouter, Azero

 real, allocatable, dimension(:) :: Mass_enc 
 real, allocatable, dimension(:) :: fM, dsigma_cloud, E_acc

!$OMP threadprivate(rhocloud, beta_cloud, amax_cloud, tff, mdot_cloud, mdot)
!$OMP threadprivate(omega_cloud, dmstar_cloud)
!$OMP threadprivate(dMcloud, dMtot, MRouter, Router, Mouter, azero)
!$OMP threadprivate(Mass_enc, fM, dsigma_cloud, E_acc)

end module clouddata
