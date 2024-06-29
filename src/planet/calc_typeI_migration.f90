subroutine calc_typeI_migration(iplanet, tmig1)
  !
  ! Computes the Type I migration timescale of a given planet
  !

  
  use planetdata
  use embryodata, only: c_mig
  use stardata, only: mstar, omega_d,cs_d,sigma_d
  use unitdata, only: pi,yr

  implicit none

  integer, intent(in) :: iplanet
  real, intent(inout) :: tmig1
  real :: hp, aspectratio,mdiscmig

  hp = cs_d(iplanetrad(iplanet))/omega_d(iplanetrad(iplanet))
  aspectratio = hp/ap(iplanet)

  mdiscmig = pi*ap(iplanet)*ap(iplanet)*sigma_d(iplanetrad(iplanet))

  ! Baruteau et al (2013) expression
  if (mdiscmig .gt. 0.0d0) then
    tmig1 = mstar*mstar*aspectratio*aspectratio/(mdiscmig*mp(iplanet)*omega_d(iplanetrad(iplanet)))
    tmig1 = c_mig*tmig1
  else
    tmig1 = 1.0d30
  endif 

!  print*, 'Migration timescale in calc Type I = ', tmig1/yr 

  ! Old Bate et al expression
  !tmig1 = aspectratio*mstar/(mp(iplanet)*omegaK(iplanetrad(iplanet)))
  

end subroutine calc_typeI_migration
