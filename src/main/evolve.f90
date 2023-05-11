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
  use clouddata

  implicit none

  integer :: i,j, iplanet, timeup,nsurvive,ifile
  real :: t, tout, tnext
  real (kind=8) :: dtmin,CO,C1,dr2
  real (kind=8) :: term1, term2
  real (kind=8) :: diskmass, diskmass_added
  real (kind=8) :: massindisc1,massindisc2,mass_added
  real (kind=8) :: angmom_disc, angmom_added, angmom_tot
  real (kind=8) :: ro, M_ro, vphi
  real (kind=8) :: masstot, radnew 
  real (kind=8) :: sigma_d_max
 
  real,allocatable, dimension(:) :: massenc, massenc_collapse

  real :: dtorque, dTcdr, vr, dt_torque
  real :: cp
  real :: mass_accr, Mdotwind
  logical :: timestepOK

  If(allocated(massenc)) deallocate(massenc)
  If(allocated(massenc_collapse)) deallocate(massenc_collapse)
  allocate(massenc(nrannuli))
  allocate(massenc_collapse(nrannuli))

  tdump = 0.0
  isnap = 0

  ! Initialise all embryos

!  do j=1,nembryo
!
!     embryo(j)%icurrent = embryo(j)%iform
!     embryo(j)%R = embryo(j)%R0
!     embryo(j)%t_spent = 0.0
!     embryo(j)%finished = 0
! 
!  enddo

  t = 1.0d-3     
  tout = 0.0d0
  dt = 0.0d0

  tnext = t_disc_dump

  Do While (t .lt. tmax*3.15d7)

    timestepOK = .false.

    ! Do timestep while checking that timestep is okay

    If (nembryo .gt. 0) Then
      call find_planets_in_disc    

      call evolve_embryos(t)
    EndIf   

    massenc(:) = 0.0d0
    massenc_collapse(:) = 0.0d0

    massenc(1) = mstar 
    massenc_collapse(1) = mstar_collapse 
    do i = isr, ier
      massenc(i) = sigma_d(i-1)*pi*(rf(i)**2.0d0 - rf(i-1)**2.0d0) + massenc(i-1)
      massenc_collapse(i) = sigma_d_collapse(i-1)*pi*(rf(i)**2.0d0 - rf(i-1)**2.0d0)
      massenc_collapse(i) = massenc_collapse(i) + massenc_collapse(i-1)
    enddo
     

    DO WHILE (.not.timestepOK)

     ! Check to see if all embryos have finished, and update global timestep

       snew(:) = 0.0d0
       Tnew(:) = 0.0d0

       call disc_properties(t) 

       If (nembryo .gt. 0) Then
         torque_term(:) = 0.0d0

         call compute_planet_torques(t)
       EndIf

       timestepOK = .true.

       CALL timestep(t)

       T_d(isr-1) = T_d(isr)
       nu_tc(isr-1) = 0.0d0
       sigma_d(isr-1) = 0.0d0
       nu_tc(ier+1) = nu_tc(ier)
       sigma_d(ier+1) = sigma_d(ier)*sqrt(rz(ier)/rz(ier+1))
       T_d(ier+1) = T_d(ier)

       heatfunc(:) = 9.0d0*nu_tc(:)*sigma_d(:)*omega_d(:)*omega_d(:)/8.0d0 
   
       If ((runmode == 'C1') .and. (t .lt. 2.0d0*tff)) Then
         call setup_wind
       EndIf
       call compute_wind

       Mdotwind = 0.0d0
       mass_accr = 0.0d0

       If (runmode == 'C1') then
         call sigma_mdot(t)
         call Eacc_calc(t)
       Else
         dsigma_cloud(:) = 0.0d0
         E_acc(:) = 0.0d0
       endif

       do i = isr, ier

          term1=rf1_2(i+1)*(nu_tc(i+1)*sigma_d(i+1)*rz1_2(i+1)-nu_tc(i)*sigma_d(i)*rz1_2(i))*drfm1(i+1)
          term2=rf1_2(i)*(nu_tc(i)*sigma_d(i)*rz1_2(i)-nu_tc(i-1)*sigma_d(i-1)*rz1_2(i-1))*drfm1(i)

!          dtorque = 0.0

         ! Symmetrised planet torque
         ! 1/2[ ((i+1) + (i)) - ((i)+(i-1)) ]
!          dtorque = 0.5*(torque_term(i+1) - torque_term(i-1))

          if (nembryo .gt. 0) Then
            if (i .le. iplanetrad(1)) then
              dtorque = (torque_term(i+1) - torque_term(i))
            endif
            if (i .gt. iplanetrad(1)) then
              dtorque = (torque_term(i) - torque_term(i-1))
            endif

            if (i .eq. isr) then
              do iplanet = 1, nembryo
                if (i .lt. iplanetrad(iplanet)) then
                  dtorque = torque_term(i+1)
                endif
                if (i .gt. iplanetrad(iplanet)) then
                  dtorque = torque_term(i)
                endif
              enddo
            endif
          endif

!          if (3.0*(term1-term2)-dtorque>0.0) then
!              dtorque = 0.0
!          endif

          do iplanet=1, nembryo
             if(i==iplanetrad(iplanet) .or. i-1==iplanetrad(iplanet))then
                dtorque = 0.0
             endif
          enddo

          if (sigma_d(i) .gt. 0.0d0) Then
            vr = -3.0*(term2)/(rf(i)*sigma_d(i))
          else
            vr = 0.0d0
          endif

          If (i .lt. ier) Then  
            dTcdr = (T_d(i+1) - T_d(i))*drzm1(i) 
          Else
            dTcdr = 0.0d0
          EndIf

          If (nembryo .eq. 0) dtorque = 0.0d0

!          print*, t/yr, rz(i)/AU, dsigma_cloud(i)*dt

          snew(i) = sigma_d(i) + rzm1(i)*drzm1(i)*(3.0*(term1-term2) - dtorque)*dt - sigdot_wind(i)*dt 
!          snew(i) = sigma_d(i) 

!          If (t/yr .gt. 1.0d5) Then
!            If ((i .gt. iplanetrad(1) - 10).and.(i.lt.iplanetrad(1)+10)) Then
!              print*, i, rz(i)/au, sigma_d(i), snew(i),nu_tc(i),alpha_d(i),term1, term2, dtorque
!            EndIf
!          EndIf

          If (runmode == 'C1') Then
            angmom_disc = snew(i)*pi*(rf(i+1)**2.0d0-rf(i)**2.0d0)*DSqrt(G*massenc(i)*rz(i))

            snew(i) = snew(i) + dsigma_cloud(i)*dt

            masstot = snew(i)*pi*(rf(i+1)**2.0d0 - rf(i)**2.0d0)
   
            mass_added = dsigma_cloud(i)*dt*pi*(rf(i+1)**2.0d0-rf(i)**2.0d0)

            ro = 4.0d0*pi*azero*rhocloud*G/3.0d0/omega_cloud**2.0d0

            If (azero .gt. 0.0d0) Then
              M_ro = omega_cloud**2.0d0*ro**2.0d0/G/azero
              M_ro = M_ro*ro**2.0d0

              vphi = DSqrt(G*M_ro/rz(i)) 
            Else
              M_ro = 0.0d0

              vphi = 0.0 
            EndIf
 
!            print*, i, DSqrt(G*massenc(i)*rf(i)), DSqrt(G*massenc(i+1)*rf(i+1)), vphi*rz(i),'  ang mom'
!            print*, i, DSqrt(G*massenc(i)*rz(i)), vphi*rf(i)
!            print*, i, DSqrt(G*massenc_collapse(i)*rz(i))
!            print*, rz(i)/AU, rf(i)/AU, rf(i+1)/AU, azero/AU
!            print*, massenc(i)/solarmass, M_ro/solarmass, dMtot/solarmass  
!            print*, massenc_collapse(i)/solarmass 
!            print*, sigma_d_collapse(i), snew(i) 
 
            angmom_added = mass_added*DSqrt(G*massenc_collapse(i)*rz(i))

            angmom_tot = angmom_disc + angmom_added

            If ((masstot .gt. 0.0d0).and.(massenc(i) .gt. 0.0d0)) Then
              radnew = (angmom_tot/masstot)**2.0d0/(G*massenc(i))

!              If (radnew .lt. rf(i)) Then
!                print*, i, rf(i)/AU, rz(i)/AU, rf(i+1)/AU, radnew/AU
!                print*, massenc(i)/solarmass, M_ro/solarmass
!              EndIf

!              print*, i, rz(i)/AU, rf(i)/AU, rf(i+1)/AU, radnew/AU, azero/AU
!              print*, ro/AU, massenc(i)/solarmass, M_ro/solarmass
              
!              j = i
!              do while ((radnew .lt. rz(j)) .and. (j .gt. isr)) 
!                print*, i,j,  radnew/AU, rz(j)/AU
!                j = j - 1
!                angmom_tot = angmom_tot + sigma_d(j)*pi*(rf(j+1)**2.0d0-rf(j)**2.0d0)*Sqrt(G*massenc(j)*rz(j))
!                masstot = masstot + sigma_d(j)*pi*(rf(j+1)**2.0d0-rf(j)**2.0d0)

!                radnew = (angmom_tot/masstot)**2.0d0/(G*massenc(i))
!              enddo 
              
!              print*, t/dt, i, radnew/AU, rz(i)/AU, rf(i)/AU
            Else
              radnew = rz(i)
            EndIf

            masstot = snew(i)*pi*(rf(i+1)**2.0d0 - rf(i)**2.0d0)

!            If (snew(i) .gt. 0.0d0) Then           
!              print*, snew(i-1), snew(i), frac1, frac2, '  1'

!              If (i .gt. isr) Then 
!                snew(i) = masstot*0.9/(pi*(rf(i+1)**2.0d0 - rf(i)**2.0d0))
!                snew(i-1) = snew(i-1) + masstot*0.1/(pi*(rf(i)**2.0d0 - rf(i-1)**2.0d0))
!              EndIf
             
!              print*, snew(i-1), snew(i)
!            EndIf
          EndIf

          if ((accr_on == 'y') .and. (nembryo .gt. 0)) then
            snew(i) = snew(i) - sigdot_accr(i)*dt
            mass_accr = mass_accr + twopi*sigdot_accr(i)*rz(i)*(rz(i+1)-rz(i))*dt
          endif
          
          Mdotwind = Mdotwind + twopi*sigdot_wind(i)*rz(i)*(rz(i+1)-rz(i)) 

          Tnew(i) = T_d(i)

          timestepOK = .true.

          if (T_d(i) .gt. 0.0d0) then  
            cp = cs_d(i)**2.0d0/(T_d(i)*1.667d0*(1.667d0-1.0d0))
          else
            cp = cs_d(i)**2.0d0/(T_source(i)*1.667d0*(1.667d0-1.0d0))
          endif
          if (((runmode=='g2').or.(runmode=='C1')).and.(alpha_d(i) .le. (alpha_visc + 1.0e-8))) then
            if ((sigma_d(i) .gt. 0.0d0).and.(cp .gt. 0.0d0)) Then
              Tnew(i) = T_d(i) + 2.0*dt*(heatfunc(i)-coolfunc(i))/(cp*sigma_d(i)) - vr*dTcdr*dt
              If (runmode == 'C1') Then 
                Tnew(i) = Tnew(i) + E_acc(i)*dsigma_cloud(i)*dt/(cp*sigma_d(i)) 
              EndIf
              If (isnan(Tnew(i))) Then
                print*, 'Nan', T_d(i), E_acc(i), dsigma_cloud(i), dt, cp, sigma_d(i)
              EndIf  
              If (Tnew(i) .gt. T_d(i)*5.0d0) Then
                timestepOK = .false.       
!                Tnew(i) = T_d(i)*1.5d0
              EndIf
              If (Tnew(i) .lt. T_d(i)*0.2d0) Then
                timestepOK = .false.      
!                Tnew(i) = T_d(i)*0.75d0
              EndIf
            else
              Tnew(i) = T_d(i)
            endif
          endif          
 
          if (Tnew(i).lt.T_source(i)) Then
            Tnew(i) = T_source(i)
          endif
          if (Tnew(i).gt.1000.0d0) Then
            Tnew(i) = 1000.0d0
          endif
          if (snew(i) .lt. 0.1d0) Then
            Tnew(i) = T_source(i)
          endif

          if (snew(i) .lt. 1.0d-12) snew(i) = 0.0d0
!          if (snew(i) .lt. 1.0d-15) snew(i) = 1.0d-15
       enddo    

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

!    call evolve_embryos(t)

    diskmass = 0.0d0
    massindisc1 = 0.0d0
    massindisc2 = 0.0d0
    mass_added = 0.0d0
    do i = isr, ier
      massindisc1 = massindisc1 + 2.0d0*pi*rz(i)*sigma_d(i)*(rf(i+1)-rf(i))
      massindisc2 = massindisc2 + 2.0d0*pi*rz(i)*snew(i)*(rf(i+1)-rf(i))

      diskmass = massindisc2

      If (runmode == 'C1') Then
        mass_added = mass_added + 2.0d0*pi*rz(i)*dsigma_cloud(i)*dt*(rf(i+1)-rf(i))
      Else
        mass_added = 0.0d0
      EndIf
    enddo

    If (runmode == 'C1') Then
      mstar = mstar + dmstar_cloud*dt
      mstar_collapse = mstar_collapse + dmstar_cloud*dt
    EndIf

    if (massindisc1-(massindisc2-mass_added) .ge. 0.0d0) then
      mstar = mstar + (massindisc1-(massindisc2-mass_added))
    endif

    dMtot = dMtot + mdot*dt

! Copy back surface density

    mdisc = 0.0d0
    sigma_d_max = 0.0d0
    do i = isr, ier
       sigma_d(i) = snew(i)
       If (sigma_d(i) .gt. sigma_d_max) sigma_d_max = sigma_d(i)
       T_d(i) = Tnew(i)

       mdisc = mdisc + twopi*sigma_d(i)*rz(i)*(rf(i+1)-rf(i))

       sigma_d_collapse(i) = sigma_d_collapse(i) + dsigma_cloud(i)*dt

!       T_d(i) = Tnew(i)
    enddo

    if ((sigma_d_max .lt. 0.1d0) .and. (t/yr .gt. 1.0d6)) Then
      exit
    endif

    T_d(isr-1) = T_d(isr)
    nu_tc(isr-1) = 0.0d0
    sigma_d(isr-1) = 0.0d0
    nu_tc(ier+1) = nu_tc(ier)
    sigma_d(ier+1) = sigma_d(ier)*sqrt(rz(ier)/rz(ier+1))
    T_d(ier+1) = T_d(ier)

    do iplanet = 1, nembryo
      If (runmode .eq. 'C1') Then
        If (t/yr .gt. (embryo(iplanet)%t_form/yr + tdelay_planettorque)) Then
          call migrate_planets
        EndIf
      Else If (t/yr .gt. (tdelay_planettorque)) Then
        call migrate_planets
      End If       
      embryo(iplanet)%a = ap(iplanet)
      mp(iplanet) = embryo(iplanet)%m
    enddo 

    t = t + dt
    tout = tout + dt

    If (tout/3.15d7 .gt. t_disc_dump/10.0d0) Then
        print*, t/3.15d7, dt/3.15d7, sigma_d(isr), mstar/solarmass, mdisc/solarmass, mass_accr/solarmass
!        print*, embryo(1)%icurrent, embryo(1)%m/mjup, embryo(1)%r/rjup, embryo(1)%a/au, embryo(1)%t_cool0/yr
!        print*, embryo(1)%T, embryo(1)%mcore/mearth, embryo(1)%rcore/rearth
!        print*, 0.8d0*au*kappa_star**0.5d0*(10.0d4*yr/t)**0.5d0/au
        tout = 0.0d0

!        write(16,*) t/3.15d7, mp(1)/mjup, ap(1)/au
 
        diskmass = 0.0d0
        do i = isr, ier
          diskmass = diskmass + 2.d0*pi*rz(i)*sigma_d(i)*(rz(i+1)-rz(i))
        enddo

        !write(27,*) t/3.15d7, mdisc/mstar
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

    If (t/yr .gt. 2.0d7) Then
      If (nembryo .eq. 0) Then
        exit
      EndIf
      If ((nembryo .eq. 1) .and. (embryo(1)%a/AU .lt. 1.0d0)) Then
        exit
      EndIf       
    EndIf  
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
