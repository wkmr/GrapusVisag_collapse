MODULE embryodata

real, parameter :: mearth = 5.97360d27
real, parameter :: mjup = 1.8986d30
real, parameter :: rearth = 6.3710d8 
real, parameter :: rjup = 7.1492d9
real, parameter :: zeta_1 = 3.65375  ! These parameters only work for n=1.5
real,parameter :: theta_grad = 5.5642e-2
real, parameter :: n = 1.5  ! Constants for n=1.5 polytropes
real, parameter :: mu = 2.4 

real, parameter :: rho_s = 1.0  ! Density of solids in g cm^-3
real, parameter :: scrit_min = 10.0 ! Minimum critical size for sedimentation (cm)
real, parameter :: s0 = 1.0e-4 ! Initial grain size (cm)
real, parameter :: vfrag = 100.0 ! Critical velocity for destructive collisions
real, parameter :: mfp0 = 40.0  ! Typical mean free path in H2 at 1e-10 g cm^3
real, parameter :: alpha_dust = 0.004 ! Turbulent viscosity parameter for the dust layer
real,parameter :: f_ice = 0.3333 ! Fraction of the solid component in ices
real,parameter :: gamma = 1.666 ! Specific heat ratio of the polytrope

real,parameter :: Tmelt = 160.0  ! Melting Temperature of ices
real, parameter :: Tvap = 1600.0 ! Vapourisation temperature of the grains
real, parameter :: Tdiss = 3000.0 ! Dissociation Temperature for H2

real, parameter :: tolerance = 1.0e-5
real, parameter :: G_nbody =1.0
real, parameter :: ecc_mu = 0.094
real, parameter :: ecc_sig = 0.095
real, parameter :: inc_mu = 0.00091
real, parameter :: inc_sig = 0.00005

integer,parameter :: istart = 10
integer,parameter :: ifinal = 11
integer,parameter :: ilog = 12
integer,parameter :: inbodylog = 20
integer,parameter :: isnapfile = 40

integer :: nembryo, nbodies,isnap,nsnaps,nzeros,finishcheck

real :: fg,kappa_0,kappa_star,rho_ad, m1,T1,mfp,dt, p_kap, p_grow
real :: c_mig,c_gap,c_collapse, maxerror,fragsep,tsnap,maxsnap
character(1) :: core_feedback,accr_on,accr_on_disc,nbody,multishot,initialecc,zerostring
character(6) :: zeroformat,snapchar

real, allocatable, dimension(:) :: sigdot_accr

character(100),allocatable,dimension(:) :: snapshotfile

type GE_embryo

integer :: iform, icurrent, finished
integer :: ivap,imelt,idiss,igrown,iself,ijeans,itidal ! Identifies current state of embryo
integer :: migtype ! What sort of migration is embryo undergoing (Type I/II)
real :: m, a,R,Rg,Rg0, Rcore, mcore, msolids
real :: Nsteps ! Number of steps left at an orbital radius
real :: R0, T0, cs0, t_cool0, t_0_cool, t_grow0,t_sed0
real :: T,L, t_cool, cs,rhoc, fg
real :: t_grow, t_sed, tmig,tgap,tcross
real :: t_spent ! Time spent at a particular orbital radius
real :: tself,rself,scrit
! N Body variables

real :: rmag,semimaj,ecc,inc,longascend,argper,trueanom


end type GE_embryo

type(GE_embryo), allocatable :: embryo(:)

! N Body variables

real,parameter :: small = 1.0e-20
real,parameter :: dampfac = 10.0
real,parameter :: rsoft = 1.0e-5

real :: system_ang, system_energy, initial_system_ang,initial_system_energy
real :: dE, dL, dt_nbody, totalmass

! Body data

real, allocatable, dimension(:,:) :: pos,vel,acc
real, allocatable,dimension(:,:) :: newpos,newvel
real,allocatable,dimension(:,:) :: angmom

real,dimension(3) :: system_angmom,rcom,vcom,acom
real,allocatable,dimension(:) :: mass, ekin,epot,etot,angmag
real,allocatable,dimension(:) :: r

!$OMP threadprivate(nembryo, nbodies,isnap,nsnaps,nzeros,finishcheck)
!$OMP threadprivate(fg,kappa_0,kappa_star,rho_ad, m1,T1,mfp,dt, p_grow)
!$OMP threadprivate(maxerror)
!$OMP threadprivate(zerostring)
!$OMP threadprivate(zeroformat,snapchar,snapshotfile)

!$OMP threadprivate(embryo,sigdot_accr)

!$OMP threadprivate(system_ang, system_energy, initial_system_ang,initial_system_energy)
!$OMP threadprivate(dE, dL, dt_nbody, totalmass,pos,vel,acc,newpos,newvel,angmom)
!$OMP threadprivate(system_angmom,rcom,vcom,acom,mass, ekin,epot,etot,angmag,r)

contains

subroutine get_magnitude(vector,magnitude)

real,dimension(3),intent(in) :: vector
real,intent(out):: magnitude

magnitude = sqrt(vector(1)*vector(1)+vector(2)*vector(2)+vector(3)*vector(3))

end subroutine get_magnitude


END MODULE embryodata
