PROGRAM grapus

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
use planetdata
use clouddata

IMPLICIT NONE

character(len=8) :: fnum
integer is

INTERFACE

   subroutine initial
   end subroutine initial

   SUBROUTINE generate_star
   END SUBROUTINE generate_star

   SUBROUTINE generate_disc
   END SUBROUTINE generate_disc

   SUBROUTINE generate_embryos
   END SUBROUTINE generate_embryos

   SUBROUTINE evolve
   END SUBROUTINE evolve

   SUBROUTINE evolve_disc
   END SUBROUTINE evolve_disc

END INTERFACE


! Initialise run
call initial

! Loop over total number of stars
!$OMP PARALLEL default(private) &
!$OMP shared(Nstar,mstar0,prefix_orig,runmode,debug,multishot,tsnap,maxsnap,alpha_visc,Lx_0) &
!$OMP shared(q_disc0,mdotvisc0) &
!$OMP shared(nbody,datafilepath,tmax,t_disc_dump,iseed,rin) &
!$OMP shared(stell_irr, Lstar, Tirr) &
!$OMP shared(Mcloud_0, Rcloud_0, f_cloud_0) &
!$OMP shared(p_kap,fragsep,initialecc,c_mig,c_gap,c_collapse,core_feedback) &
!$OMP shared(accr_on, accr_on_disc)
!$OMP do schedule(dynamic)
  DO is=1,Nstar
     istar = is

     write (fnum,'(I5.5)') istar   

     prefix = trim(prefix_orig)//trim(fnum)

     print*, prefix

     ! Create a T_Tauri star
     CALL generate_star

     ! Generate a self-consistent self-gravitating disc
     CALL generate_disc

     print*, 'generate disc finished'

!     if (rout/AU.gt.50 .and. rout/AU.lt.500) then
!      HERE
     call setup_wind

!       call generate_embryos

!       call setup_planets

       ! Where possible, generate embryos from fragmenting disc
       ! Subroutine outputs initial embryo M, R, a,  to file

       ! Evolve these embryos towards the inner disc
       ! Subroutine outputs final M,R,a to file
!     IF(nembryo>0) CALL evolve
      snapshotcounter = 1

      print*, nembryo

      If (nembryo.ge.0) then
!        call write_dump(0.0) 
        CALL evolve
      else
        call nbody_deallocate_arrays 
      endif

      call nbody_deallocate_arrays
!    else
!      print*, 'ERROR: Disc radius too small/large, trying again...'
!      call nbody_deallocate_arrays 
!    endif

     print*, 'closing logfile'

     close(itime)

    !If (nplanet .gt. 0) CALL evolve
    !If ((nplanet .eq. 0) .and. (nbody=='y')) call nbody_deallocate_arrays 
     print*, ''
     print*, ''
     print*, ''
     print*, ''

  ENDDO

!$omp end do
!$omp end parallel

  close(istart)
  close(ifinal)
  close(ilog)

END PROGRAM grapus