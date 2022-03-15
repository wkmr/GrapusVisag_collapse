module unitdata
!-------------------------------------------------
! Module stores all physical constants, EOS table,
! and other parameters relating to file storage
!--------------------------------------------------
  
  integer :: itime = 80
  integer :: iprof = 83
  integer :: ispec = 178
  integer :: itorque = 275
  integer :: snapshotcounter
  real :: tdump, t_disc_dump, trun,nfiles

  real(kind=8), parameter :: G = 6.67d-8           ! cgs
  real(kind=8), parameter :: AU = 1.496e13	   ! AU in cm
  real(kind=8), parameter :: yr = 3.15e7	   ! year in s
  real(kind=8), parameter :: solarmass = 1.99e33   ! solar mass
  real(kind=8), parameter :: msolyr = solarmass/yr

  real(kind=8), parameter :: Qcrit = 1.5	! Critical Toomre Parameter

  real(kind=8), parameter :: pi =3.1415926535
  real(kind=8), parameter :: twopi = 2.0*pi
  real(kind=8), parameter :: pibytwo = 0.5*pi
  real(kind=8), parameter :: stefan = 5.67d-5    ! Stefan-Boltzmann constant
  real(kind=8), parameter :: planck = 6.626d-27  ! Planck's constant
  real(kind=8), parameter :: pc = 3.08d18        ! parsec			
  real(kind=8), parameter :: k_B = 1.38d-16  	 ! Boltzmann constant
  real(kind=8), parameter :: m_H = 1.67d-24	 ! mass of hydrogen atom

  character(len=2) :: runmode, layerchoice, tempchoice, zerostring
  character(len=6) :: snapshotformat
  character(len=100) :: fileno,prefix_orig,prefix
  character(len=25), parameter :: paramfile = 'grapusvisag.params'

  !$omp threadprivate(prefix,fileno,snapshotcounter,tdump,nfiles)

end module unitdata
