SUBROUTINE evolve
  !****************************************************************************
  ! Subroutine evolves both the embryo and disc population simultaneously 
  ! This routine evolves embryos simultaneously on a shared timestep
  !****************************************************************************

  use stardata
  use embryodata
  use eosdata
  use winddata, only: sigdot_wind
  use planetdata
  use unitdata

  implicit none

  integer :: i,j, iplanet, timeup,nsurvive,ifile
  real :: t, tout, tnext
  real (kind=8) :: dtmin,CO,C1,dr2
  real (kind=8) :: term1, term2
  real :: dtorque, dTcdr, vr, dt_torque
  real :: Mdotwind
  logical :: timestepOK

  tdump = 0.0
  isnap = 0

  ! Initialise all embryos

  do j=1,nembryo

     embryo(j)%icurrent = embryo(j)%iform
     embryo(j)%R = embryo(j)%R0
     embryo(j)%t_spent = 0.0
     embryo(j)%finished = 0
 
  enddo

  t = 0.0d0     
  tout = 0.0d0
  dt = 0.0d0

  tnext = t_disc_dump

  Do While (t .lt. tmax*3.15d7)

    timestepOK = .false.

    ! Do timestep while checking that timestep is okay

!    call evolve_embryos(t)
    
    DO WHILE (.not.timestepOK)

     ! Check to see if all embryos have finished, and update global timestep

       snew(:) = 0.0d0
       Tnew(:) = 0.0d0

       call disc_properties 

       call compute_planet_torques(t)

       timestepOK = .true.

       CALL timestep(t)

       nu_tc(isr-1) = nu_tc(isr)
       sigma_d(isr-1) = sigma_d(isr)
       nu_tc(ier+1) = nu_tc(ier)
       sigma_d(ier+1) = sigma_d(ier)*sqrt(rz(ier)/rz(ier+1))

       call compute_wind

       Mdotwind = 0.0d0

       do i = isr, ier
       
          term1=rf1_2(i+1)*(nu_tc(i+1)*sigma_d(i+1)*rz1_2(i+1)-nu_tc(i)*sigma_d(i)*rz1_2(i))*drfm1(i+1)
          term2=rf1_2(i)*(nu_tc(i)*sigma_d(i)*rz1_2(i)-nu_tc(i-1)*sigma_d(i-1)*rz1_2(i-1))*drfm1(i)

          dtorque = 0.0

         ! Symmetrised planet torque
         ! 1/2[ ((i+1) - (i)) + (i)-(i-1)) ]
          dtorque = 0.5*(torque_term(i+1) - torque_term(i-1))

          if(3.0*(term1-term2)-dtorque>0.0) then
              dtorque = 0.0
          endif
           
          do iplanet=1,nplanet
             if(i==iplanetrad(iplanet) .or. i-1==iplanetrad(iplanet))then
                dtorque = 0.0
             endif
          enddo

          vr = -3.0*(term2)/(rf(i)*sigma_d(i))

!          if (t .gt. 3.0d9) print*, 1.0d0/rzm1(i)/AU, dtorque, ap(1)/AU, ap(2)/AU, ap(3)/AU
 
!           sigdot_wind(i) = 0.0

!          print*, i, term1, term2, sigdot_wind(i), dtorque, dt

          snew(i) = sigma_d(i) + rzm1(i)*drzm1(i)*(3.0*(term1-term2) - dtorque)*dt - sigdot_wind(i)*dt

          Mdotwind = Mdotwind + twopi*sigdot_wind(i)*rz(i)*(rz(i+1)-rz(i)) 

          if (snew(i) .lt. 0.0d0) snew(i) = 0.0d0

          Tnew(i) = T_d(i)
       enddo    

       timestepOK = .true.
       do i = isr, ier
         if (snew(i) .lt. 0.0d0) timestepOK = .false.
       enddo

       if (.not.timestepOK) then
         print*, 'Reducing maximum allowed timestep ', maxstep, maxstep/10.0d0
         maxstep = maxstep/10.0d0

         if(maxstep<1.0e-10) then
            print*, 'MAXIMUM TIMESTEP TOO LOW: ending run'
            stop
         endif
       endif 

    ENDDO
! Copy back surface density

    mdisc = 0.0d0
    do i = isr, ier
       sigma_d(i) = snew(i)

       mdisc = mdisc + twopi*sigma_d(i)*rz(i)*(rz(i+1)-rz(i)) 

       T_d(i) = Tnew(i)
    enddo

    call evolve_embryos(t)

    If (t/3.15d7 .gt. 1.0d3) call migrate_planets 

    do iplanet = 1, nplanet
      embryo(iplanet)%a = ap(iplanet)
      mp(iplanet) = embryo(iplanet)%m

!      print*, iplanet, embryo(iplanet)%a/AU, embryo(iplanet)%m
    enddo 

    t = t + dt
    tout = tout + dt

    If (tout/3.15d7 .gt. 10000.0d0) Then
        print*, t/3.15d7, dt/3.15d7, sigma_d(50), mstar, mdisc/mstar
        tout = 0.0d0
        write(27,*) t/3.15d7, mdisc/mstar
    EndIf

    If (t/3.15d7 .gt. tnext) Then
      snapshotcounter = snapshotcounter + 1
      call write_dump(t)
      tnext = t_disc_dump*snapshotcounter     
      print*, tnext
    endif

    tdump = tdump + dt

   !********************************************************************
   ! If necessary, write a snapshot of the population at this timestep
   !********************************************************************

    if(tdump>tsnap) then
       isnap = isnap +1

       if(debug=='y') write(*,'(A,1P,3e18.4)'), 't, dt, mdisc/mstar: ', t/yr, dt/yr,q_disc
         
       if(isnap.le.nsnaps) then
          ifile = isnapfile+isnap
          nsurvive = 0
               
          OPEN(ifile,file=snapshotfile(isnap),form='formatted',status='old',position='append')

          call write_population_snapshot(ifile,nsurvive)
       endif

       tdump = 0.0
    endif

  enddo   

     ! If all embryos have finished, exit the loop
!     IF(finishcheck==1) exit

!     t = t+dt

     ! Now evolve the disc for this timestep
 
!     IF(mstar<mdisc.or.t<dt.or.t/yr>tmax) exit 

     ! Evolve the disc
     ! Check for the end of the disc simulation using timeup

!     timeup = 0
!     CALL evolve_disc(t,dt,timeup)

!     tdump = tdump + dt

     ! If debugging, write nbody data to separate output files
!     if(nbody=='y' .and.debug=='y') call nbody_output(t) ! Debug line - check nbody outputs

!     IF(timeup==1) exit

!     i=0
     !
     ! Recalculate ancillary variables
     ! (omega, cs, H)
     !
!     q_disc = mdisc/mstar
     !print*, mstar/umass, mdisc/umass, q_disc

!     DO WHILE(i < irout)
!        i=i+1
!        omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

!        IF(alpha_d(i)/=0.0.and.sigma_d(i)/=0.0) THEN
!           cs_d(i) = sqrt(ABS(nu_d(i)*omega_d(i)/alpha_d(i)))
!        ELSE
!           cs_d(i) = 0.0
!        ENDIF
!        H_d(i) = cs_d(i)/omega_d(i)
!     ENDDO

!ENDDO

!********************************
! Evolution of system complete
!******************************** 

!
! Print out data pertaining to all surviving embryos
!

if(debug=='y') then
   print*, 'Resulting ', nembryo, ' Objects:'
   print*, '(istar, iproperties, a,e,i,m,r,rcore,mcore)'

   do j=1,nembryo
      WRITE(*,'("Run ",I5, ", Embryo ", I2,": ",7I2,7F12.4)') istar, j, &
           embryo(j)%imelt, &
           embryo(j)%ivap,embryo(j)%idiss, embryo(j)%igrown, &
           embryo(j)%iself, embryo(j)%ijeans, embryo(j)%itidal, &
           embryo(j)%a/udist, embryo(j)%ecc, embryo(j)%inc, &
           embryo(j)%M/Mjup, embryo(j)%R/Rjup, &
           embryo(j)%mcore/mearth, embryo(j)%rcore/rearth
   enddo
endif


nsurvive = 0

! Write final population
call write_population_snapshot(ifinal,nsurvive)

write(*,'(A,I1)') 'Number of survivors: ',nsurvive
call flush(ifinal)

! Deallocate nbody arrays ready for next run
if(nbody=='y') call nbody_deallocate_arrays

END SUBROUTINE evolve
