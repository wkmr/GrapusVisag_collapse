subroutine initial

  ! ******************************************
  ! Reads in parameter file and sets up run
  ! *******************************************

  ! Program carries out population synthesis of gravitational instability 
  ! using the tidal downsizing paradigm of Nayakshin et al 
  ! see:
  ! Forgan & Rice 2013, MNRAS, 432, pp 3168-3185 (v1.0)
  ! Forgan et al, 2018, MNRAS, 474, pp 5036-5048  (v2.0: adding N body physics)
  !
  
  use stardata
  use embryodata
  use eosdata
  use unitdata
  use winddata
  use clouddata

  IMPLICIT NONE


  ! Print introductory text

  print*, '*****************************************************************'
  print*, '*   GRAPUS: GRavitational instability PopUlation Synthesis'
  print*, '*'
  print*, '*   Population Synthesis of Self-Gravitating  '
  print*, '*   Disc Fragmentation and Tidal Downsizing   '
  print*, '*'
  print*, '*   Written by Duncan Forgan (dh4gan)'
  print*, '*'
  print*, '*   v1.0: 2013 - see Forgan & Rice 2013, MNRAS, 432, pp 3168-3185 '
  print*, '*   v2.0: 2017 - addition of C parameters and N Body engine'
  print*, '*                (Forgan et al 2018, MNRAS, 474, pp 5036-5048)'
  print*, '*'   
  print*, '*   For best results, compile in gfortran '
  print*, '*   Reads inputs from grapusvisag.params '
  print*, '*   File path to EOS and disc files contained therein '
  print*, '*'
  print*, '*****************************************************************'


  print*, 'Reading parameters'
  !		Read in input parameters

  open(10, file='grapusvisag.params', status='unknown')
  read(10,*) prefix_orig           ! Output file prefix
  read(10,'(a2)') runmode          ! g1=self-gravitating 1, g2=self-gravitating 2
  read(10,*) debug                 ! Run in debug mode? (y/n)
  read(10,*) multishot             ! Multiple time snapshots of popn? (y/n)
  read(10,*) tsnap                 ! Time between population snapshots
  read(10,*) maxsnap               ! Popn time snapshots not made after this time (except final snapshot!)
  read(10,*) alpha_visc0           ! Value of constant alpha in no GI regime
  read(10,*) alpha_frag            ! Fragmentation boundary represented in terms of alpha
  read(10,*) MJeansdot             ! Include rate of change of Jeans mass (-10 > MdotJ > 0) (y/n)
  read(10,*) fragprob              ! One minus probability of fragmentation in fragmentation region (i.e., 0.95 means 1 in 20).
  read(10,*) Lx_0                  ! X Ray luminosity (in ergs s-1)   
  read(10,*) nbody                 ! Use N Body integrator? (y/n)
  read(10,*) Nstar                 ! Number of star systems to simulate
  read(10,*) mstar0                ! Mid-point of range of star masses (can be up to 1/8xmstar0 higher) if in g mode.
  read(10,*) q_disc0               ! Lower limit to disc-to-star mass ratio (range from q_disc0 to q_disc0 + 0.2) if in g mode.
  read(10,*) mdotvisc0             ! Lower limit to log10 of mass accretion rate (range from mdotvisc0 to mdotvisc0 + 1.5) if in g mode. 
  read(10,*) stell_irr             ! Assume stellar irradiated disc (n = constant background Tirr)
  read(10,*) Lstar                 ! Luminosity of central star in solar luminosities
  read(10,*) Tirr0                 ! Constant background irradiation temperature (in K) if not using stellar irradiation
  read(10,*) datafilepath          ! File path to location of disc file
  read(10,*) Mcloud_0              ! Mass of cloud in collapsing cloud runs (solar masses) - C mode.
  read(10,*) Rcloud_0              ! Radius of cloud (pc) - C mode.
  read(10,*) f_cloud_0             ! f, for rotation of cloud - C mode.    
  read(10,*) t_frag                ! time after which fragmentation can occur - C mode.
  read(10,*) tmax                  ! Maximum runtime of the model (in years)
  read(10,*) t_disc_dump           ! Time between dumps for disc model
  read(10,*) iseed                 ! Random number seed
  read(10,*) rin                   ! Disc inner boundary
  read(10,*) dr                    ! Radial separation
  read(10,*) p_kap                 ! opacity index (kappa = kappa_0 T^p_kap)
  read(10,*) fragsep               ! Maximum frag separation (Hill Radii)
  read(10,*) initialecc            ! Give fragments non-zero initial e/i? (y/n)
  read(10,*) c_mig                 ! Migration Factor (tmig = c_mig*tmig)
  read(10,*) c_gap                 ! Gap opening factor (tgap = c_gap*tgap)
  read(10,*) c_collapse            ! Collapse timescale factor
  read(10,*) core_feedback         ! Radiative Feedback of Core Formation? (y/n)
  read(10,*) accr_on               ! Mass accretion onto protoplanet (y/n)
  read(10,*) accr_on_disc          ! Include fraction of mass accretion through disc (y/n)
  close(10)

!  q_disc0 = q_disc0 + 0.1

  print*, accr_on

  print*, runmode

  print*, 'Parameter file read complete', tdump

  rin = rin*udist
  dr = dr*udist

  tsnap = tsnap*yr
  maxsnap = maxsnap*yr

  p_grow = (1.0+p_kap)/(2.5+p_kap)

  iseed = -abs(iseed)

  print*, 'iseed 1 = ', iseed 

  ! Disc model counter (only used when interpolating from file)
  imodel = 0

  !
  ! Set up output files for initial and final parameters
  !

  OPEN(istart, file=TRIM(prefix_orig)//'.initial', status='unknown')
  OPEN(ifinal, file=TRIM(prefix_orig)//'.final',status='unknown')
  OPEN(ilog, file=TRIM(prefix_orig)//'.log', status='unknown')

  !
  ! Set up files for snapshots
  !
  nsnaps = INT(maxsnap/tsnap)+1
  call get_zero_padding_format(nsnaps,zeroformat)

  allocate(snapshotfile(nsnaps))

  do isnap=1,nsnaps

     write(snapchar,zeroformat) isnap

     snapshotfile(isnap) = trim(prefix_orig)//'.snap.'//trim(snapchar)
     open(isnapfile+isnap,file=snapshotfile(isnap),status='unknown',form='formatted')

     ! Write the time at the beginning of every snapshot file
     write(isnapfile+isnap,*) tsnap*isnap
  enddo

  ! Read in Equation of State for self-gravitating disc
  CALL eosread


end subroutine initial
