module planetdata
! Contains all data relating to planets inside the disc

  real, parameter :: tdelay_planettorque = 1.0e4 ! Time delay for torques to activate in disc (years)

  ! Variables
  integer :: nplanet, nactive
  
  real(kind=8) :: rremove, p_create
character(100) :: planetfile
character(1) :: planetchoice

  ! Arrays
  integer, allocatable,dimension(:)  :: alive, iplanetrad
  real, allocatable, dimension(:)  :: mp,ap, total_planet_torque, torque_term
  real, allocatable, dimension(:) :: adot,tmig,tmigI, fII
real, allocatable,dimension(:,:) :: lambdaI, lambdaII
real, allocatable, dimension(:,:) :: torquei

!$omp threadprivate(nplanet,nactive,rremove,p_create,planetfile,planetchoice)
!$omp threadprivate(alive,iplanetrad,mp,ap,total_planet_torque,torque_term,adot)
!$omp threadprivate(tmig,tmigI,fII,lambdaI,lambdaII,torquei)

end module planetdata
 
