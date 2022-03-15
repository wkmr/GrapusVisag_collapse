MODULE stardata
  integer, parameter :: nrannuli = 300
  real,parameter :: tolerance = 1.0e-4

  integer,parameter :: idisc = 5
  integer,parameter :: nvar = 7
!  integer,parameter:: ntime_mod = 95
!  integer,parameter :: ntime_mod = 450
!  integer,parameter :: nrad_mod = 200

  integer :: nmodels,ntime_mod,nrad_mod
  integer :: Nstar, istar,ntries,irout,irfrag, iseed, imodel
  integer :: isr, ier
  real :: mstar,mdisc,q_disc,dr,rin,rout, rtrunc, rtruncmax
  real :: mdotvisc,rhomid,rfrag,rmax,tmax
  real :: maxstep
  character(100) :: discfile,datafilepath
  character(1) :: truncate_disc,debug

  !	Zone centred arrays
  real, allocatable, dimension(:) :: rz,rzm1,rz1_2,drzm1
  !	Face centred arrays
  real, allocatable, dimension(:) :: rf,rf1_2,drfm1

  real,allocatable,dimension(:) :: sigma_d, cs_d, omega_d, gamma_J
  real,allocatable,dimension(:) :: betac_d, mjeans, ljeans, r_d, H_d,alpha_d, alpha_g
  real,allocatable, dimension(:) :: T_d, kappa_d, gamma_d, tau_d,nu_d,nu_tc
  real,allocatable, dimension(:) :: T_source, sigdot, coolfunc, Q, tcool
  real,allocatable, dimension(:) :: snew, Tnew
  real,allocatable, dimension(:) :: timemod,mstar_mod
  real,allocatable,dimension(:,:) :: r_mod,sigma_mod,nu_mod,Tc_mod
  real,allocatable,dimension(:,:) :: alpha_mod,tau_mod  


  real,dimension(nvar) :: holder




CONTAINS

  subroutine get_zero_padding_format(nfiles,zeroformat)
    implicit none

    integer,intent(in) :: nfiles
    integer :: nzeros
    character(1) :: zerostring
    character(6),intent(inout) :: zeroformat
    
    nzeros = int(log10(real(nfiles+1)))+2
    write(zerostring,'(I1)') nzeros
    zeroformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"
    
    return 
  end subroutine get_zero_padding_format


  FUNCTION ran2(idum)
    INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
    REAL AM,EPS,RNMX
    PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1, &
         IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791, &
         NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
    INTEGER idum2,j,k,iv(NTAB),iy
    SAVE iv,iy,idum2
    DATA idum2/123456789/, iv/NTAB*0/, iy/0/
    if (idum.le.0) then
       idum=max(-idum,1)
       idum2=idum
       do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11        continue
          iy=iv(1)
       endif
       k=idum/IQ1
       idum=IA1*(idum-k*IQ1)-k*IR1
       if (idum.lt.0) idum=idum+IM1
       k=idum2/IQ2
       idum2=IA2*(idum2-k*IQ2)-k*IR2
       if (idum2.lt.0) idum2=idum2+IM2
       j=1+iy/NDIV
       iy=iv(j)-idum2
       iv(j)=idum
       if(iy.lt.1)iy=iy+IMM1
       ran2=min(AM*iy,RNMX)
       return

     END FUNCTION

   END MODULE stardata
