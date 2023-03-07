SUBROUTINE generate_embryos(t)
  ! Subroutine finds fragment radius, and then calculates fragment masses at equal log a spacing

  use stardata
  use embryodata
  use eosdata
  use unitdata

  implicit none


  integer :: i,j,ibody
  real :: kappa,r_hill,rtest,exp1,exp2
  real :: t, rnum
  real, dimension(100) :: cspace


  ! Debug variables - for nbody output files
  character(100) :: outputfile
  character(6) :: filenumformat
  character(10) :: file_no,runno

  IF(allocated(sigdot_accr)) deallocate(sigdot_accr)
  allocate(sigdot_accr(nrannuli))

  i=0

  DO WHILE(i<nrannuli) 
     i=i+1

     IF ((alpha_d(i) .gt. alpha_frag).and.(rz(i)/AU .lt. 400.0d0)) THEN
       call random_number(rnum)
       If (MJeansdot == 'y') Then
        IF(gamma_j(i)>-10.0 .and. gamma_j(i)<0.0) THEN
          If (rnum .gt. fragprob) Then
            rfrag = r_d(i)
            irfrag = i
            exit
          EndIf
        EndIf
       EndIf
       If (MJeansdot == 'n') Then
         If (rnum .gt. fragprob) Then
           rfrag = r_d(i)
           irfrag = i
           exit
         EndIf
       EndIf
     EndIf
  ENDDO

  IF(i==nrannuli) THEN
     IF(allocated(embryo)) deallocate(embryo)
     nembryo=0
     return
  ELSE
     write(*,'(A,X,F10.2,A,X,I5)') 'Fragmentation radius is ',rfrag/udist, ' AU: annulus ',irfrag
  ENDIF

  ! Generate number of possible embryos
  ! Embryo spacing is C Hill Radii, where C uniformly distributed
  ! between 1.5 and 'fragsep' (stored for later use)
  
  rtest = rfrag
  i=irfrag
  nembryo=1

  j = 0
  DO WHILE ((i < irout) .and. (alpha_d(i) .gt. alpha_frag)) 
     j =j+1
     cspace(j) = 1.5 + ran2(iseed)*(fragsep-1.5)

     IF(r_d(i)>rmax) exit

     r_hill = rtest*(mjeans(i)/(3.0*mstar))**0.333

     DO WHILE (r_d(i) < rtest+cspace(j)*r_hill .and. i<irout)
        i=i+1
     END DO

     if (i .lt. irout) Then
       rtest = r_d(i)
       nembryo = nembryo+1
     endif
 
  ENDDO
  
  write(*,'(A,I1,A)') 'There are ',nembryo, ' embryos'

  ! Allocate embryo type array

  IF(allocated(embryo)) deallocate(embryo)
  allocate(embryo(nembryo))

  ! Calculate initial embryo properties

  i=irfrag
  rtest =rfrag

  totalmass = mstar

  DO j=1,nembryo

     embryo(j)%m = mjeans(i)
!     embryo(j)%m = 10.0d0*mjup
     embryo(j)%a = rtest 
     embryo(j)%semimaj = embryo(j)%a/udist
     embryo(j)%iform = i
     embryo(j)%icurrent = i

     embryo(j)%t_spent = 0.0
     embryo(j)%finished = 0

     embryo(j)%t_form = t

     totalmass = totalmass + embryo(j)%m

    ! If this is an n-body run, then set up orbital data here

     ! If using a non-zero initial eccentricity/inclination, sample from
     ! Gaussians here
     if(initialecc=='y') then

        call sample_gaussian(embryo(j)%ecc, ecc_mu, ecc_sig)
        if(embryo(j)%ecc <0.0) embryo(j)%ecc = 0.0

        call sample_gaussian(embryo(j)%inc, inc_mu,inc_sig)
        if(embryo(j)%inc <0.0) embryo(j)%inc = twopi - embryo(j)%inc
     else

        embryo(j)%ecc = 0.0
        embryo(j)%inc = 0.0
     endif

    embryo(j)%longascend = 0.0
    embryo(j)%argper = 0.0
    embryo(j)%trueanom = twopi*ran2(iseed)

     exp1 = (1.0-n)/n
     exp2 = (3.0-n)/n

     m1 = embryo(j)%m/(0.01*umass)
     T1 = T_d(i)/10.0

!     kappa_0 = kappa_d(i)/(T1)**p_kap
     kappa_0 = 0.1d0*T1**p_kap         ! 0.1 is the assumed cgs opacity at 10K
     kappa_star = kappa_0/0.01

     !T1 = 1.0
     !kappa_star = 1.0

     rho_ad = 5.0e-13*kappa_star**(-0.666)*(T1)**((4.0-2.0*p_kap)/3.0)

!     embryo(j)%R0 = 17.5*udist*m1**(-0.333)*T1*(1.0e13*rho_ad)**(-0.666)
!     embryo(j)%R0 = 6.0*udist*m1**(-0.33333)*T1**((1.0+4.0*p_kap)/9.0)*kappa_star**(4.0/9.0)
!     embryo(j)%R0 = 6.0d0*udist*m1**(-0.33333d0)*kappa_star**(4.0d0/9.0d0)
     embryo(j)%R0 = ljeans(i)
!     embryo(j)%R0 = rjup 
!     embryo(j)%T0 = 146.0*m1**1.3333*T1**(-(1.0+4.0*p_kap)/9)*kappa_star**(-4.0/9.0)
     embryo(j)%T0 = G*mjeans(i)*2.3d0*m_H/k_B/embryo(j)%R0/3.0d0
     embryo(j)%cs0 = gamma_d(i)*k_B*embryo(j)%T0/(mu*m_H)
     embryo(j)%cs0 = sqrt(embryo(j)%cs0)

     print*, embryo(j)%R0/au, ljeans(i)/au, alpha_d(i), embryo(j)%t_form

     embryo(j)%R = embryo(j)%R0
     embryo(j)%rhoc = embryo(j)%M/(4.0*pi*theta_grad*embryo(j)%R**3)

     embryo(j)%RG = embryo(j)%R0
     embryo(j)%fg = fg
     embryo(j)%msolids = embryo(j)%fg*embryo(j)%m
     embryo(j)%rcore = 0.0  ! No solid core formed yet
     embryo(j)%tself = 0.0
     embryo(j)%rself = 0.0

    ! Initialise state of embryo in terms of evolutionary phases

     embryo(j)%iself = 0
     embryo(j)%ivap = 0
     embryo(j)%imelt = 0
     embryo(j)%idiss = 0
     embryo(j)%igrown = 0
     embryo(j)%ijeans = 0
     embryo(j)%itidal = 0

     ! If embryo hot enough, ice can melt, grains will evaporate or H2 dissociates
     IF(embryo(j)%T0> Tmelt) embryo(j)%imelt = 1
     IF(embryo(j)%T0> Tvap) embryo(j)%ivap = 1
     IF(embryo(j)%T0> Tdiss) embryo(j)%idiss = 1

     kappa = kappa_d(i)

     IF(p_kap==1.0) THEN
        embryo(j)%t_cool0 = 380.0*yr*m1**0.666*T1**(-1.333)*kappa_star**(1.0/9.0)
     ELSE IF(p_kap==2.0) THEN
        embryo(j)%t_cool0 = 5700.0*yr*m1**2*T1**(-2.1111)*kappa_star**(-0.333)
     ENDIF

!     embryo(j)%t_0_cool = 0.0 

     ! Multiply by tunable parameter c_collapse
     embryo(j)%t_cool0 = c_collapse*embryo(j)%t_cool0

     mfp = 1.0e10*mfp0*embryo(j)%m/(embryo(j)%R0)**3

     ! Calculate critical grain size at which turbulence loses to sedimentation

     embryo(j)%scrit = -mfp + sqrt(mfp*mfp +6.0*mfp*alpha_dust*embryo(j)%M/(pi*embryo(j)%R0*embryo(j)%R0*rho_s))
     embryo(j)%scrit = embryo(j)%scrit/2.0

     embryo(j)%scrit = max(scrit_min, embryo(j)%scrit)
    ! embryo(j)%scrit = 10.0
     ! Calculate grain growth timescales for this critical size

     embryo(j)%t_grow0 = 3.0*embryo(j)%cs0*embryo(j)%R0*embryo(j)%R0/(embryo(j)%m*pi*fg*G)*log(embryo(j)%scrit/s0)

     ! Modify this growth timescale due to the core's contraction

     embryo(j)%t_grow = 1.0+embryo(j)%t_grow0*(2.5+p_kap)/embryo(j)%t_cool0
     embryo(j)%t_grow = embryo(j)%t_grow**(p_grow) -1
     embryo(j)%t_grow = embryo(j)%t_cool0*embryo(j)%t_grow/(1.0+p_kap)

     ! Calculate sedimentation timescale for this critical grain size

     embryo(j)%t_sed0 = 3.0*embryo(j)%cs0*mfp/(4.0*pi*G*rho_s*embryo(j)%scrit*(embryo(j)%scrit+mfp))
     embryo(j)%t_sed = embryo(j)%t_sed0

     ! Core mass, radius

     embryo(j)%mcore = 0.0
     embryo(j)%rcore = 0.0
	
     ! Now find location of next embryo

     r_hill = rtest*(mjeans(i)/(3.0*mstar))**0.333

     DO WHILE (r_d(i) < rtest+cspace(j)*r_hill .and. i<irout)
        i=i+1
     END DO

     rtest = r_d(i)

     write(istart,'(1P,8E18.10)') embryo(j)%a/udist, embryo(j)%m/mjup, &
          embryo(j)%R0/rjup, embryo(j)%T0, embryo(j)%scrit, &
          embryo(j)%t_cool0/yr, embryo(j)%t_grow/yr, embryo(j)%t_sed0/yr

     r_hill = embryo(j)%a*(embryo(j)%m/(3.0*mstar))**0.333

     if(debug=='y') then
        write(*, '(I2,1X,I4,1P,6E18.4)') j, embryo(j)%iform, embryo(j)%m/mjup,&
             embryo(j)%a/udist, embryo(j)%R/rjup,&
             embryo(j)%t_cool0/yr, embryo(j)%t_grow/yr, embryo(j)%t_sed0/yr
     endif

  ENDDO

! Write data referring to this star-disc-planets system to log file

  WRITE(ilog,'(I6,5E18.10, I3)') istar, mstar/umass, mdisc/umass, q_disc, rout/udist, rfrag/udist, nembryo
  call flush(istart)
  call flush(ilog)


! If this is an N Body run, then create arrays for N body calculation
! Easier to do the N body calculation in separate arrays (which include star)
! N Body calculation units: (M=Msol, r=AU, t=1 yr/(2pi))

! Remember iembryo and ibody exclude/include star respectively

if(nbody=='y') then
    nbodies = nembryo+1 ! Must include the star

    allocate(pos(3,nbodies),vel(3,nbodies),acc(3,nbodies), mass(nbodies))
    allocate(angmom(3,nbodies),angmag(nbodies))
    allocate(ekin(nbodies),epot(nbodies),etot(nbodies))
    allocate(newpos(3,nbodies),newvel(3,nbodies))

    pos(:,:) = 0.0
    vel(:,:) = 0.0
    acc(:,:) = 0.0

    angmom(:,:) = 0.0
    angmag(:) = 0.0
    ekin(:) = 0.0
    epot(:) = 0.0
    etot(:) = 0.0

    newpos(:,:) = 0.0
    newvel(:,:) = 0.0

    totalmass = totalmass/umass
    dt_nbody = 1.0e-3 ! Set arbitrary small timestep initially
    mass(1) = mstar/umass

    do ibody=2,nbodies
    	mass(ibody) = embryo(ibody-1)%m/umass 
    enddo

    ! Debug lines - open N Body files to check orbits

    if(debug=='y') then
       ! First set up character for run ID

       call get_zero_padding_format(Nstar,filenumformat)
       write(runno,filenumformat) istar
       
       ! Now character for embryo IDs

       call get_zero_padding_format(nbodies,filenumformat)
       
       ! Open output files
       do ibody=2,nbodies          
          write(file_no, filenumformat) ibody          
          outputfile = TRIM(prefix)//".nbodyrun."//TRIM(runno)//".body."//TRIM(file_no)          
          open(ibody+inbodylog,file=outputfile, form="formatted")
       enddo
       
       ! Now set up log file
       outputfile = TRIM(prefix)//".nbodyrun."//TRIM(runno)//".log"

       open(inbodylog,file=outputfile,form="formatted")
    endif

    call calc_vector_from_orbit

endif


END SUBROUTINE generate_embryos
