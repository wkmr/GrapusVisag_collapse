module winddata


real(kind=8), parameter :: rwind_xray = 100.0;  ! Maximum radius for X Ray winds (AU)
real(kind=8) :: Lx_0, Lx ! Xray luminosity, wind position
real(kind=8) :: windnorm, mdot_wind

real, allocatable, dimension(:) :: sigdot_wind,sigdot_accrete

!$omp threadprivate(windnorm,mdot_wind,sigdot_wind,sigdot_accrete,Lx)

end module winddata
