subroutine compute_planet_torques(t)
  !
  ! For each planet, compute the torque it exerts on the disc
  ! cf Nayakshin (2014)
  !

  use planetdata
  use stardata
  use unitdata

  implicit none


  integer :: iplanet,i
  real :: rhill,mratio,deltap, Pcrit, typeInorm
  real :: tmig1, lambda_dash, deltamax,aspectratio, softenfactor
  real :: gap_crit
  real :: t
  logical :: soften

  call find_planets_in_disc

  torquei(:,:) = 0.0
  total_planet_torque(:) = 0.0
  torque_term(:) = 0.0

  ! Loop over each planet

  do iplanet =1,nplanet

     mratio = mp(iplanet)/mstar
     rhill = ap(iplanet)*(mratio/3.0)**0.333

     typeInorm = 0.0    

     !**************************************************************
     ! Compute Type II and Type I torques (appropriately normalised)
     !**************************************************************

     do i = isr,ier

        !*****************************************************
        ! Compute the Type II specific torque at this radius
        !*****************************************************

        deltap = abs(rz(i)-ap(iplanet))  

        deltamax = H_d(i)

        soften = .true.
        if(1.44*rhill>deltamax) deltamax = 1.44*rhill
!        if (rhill > deltamax) deltamax = rhill 
       
        if(deltap<deltamax) then
           deltap =deltamax
           soften = .true.
!           soften = .false.
        endif

        lambdaII(iplanet,i) = 0.0d0
        If (deltap .ne. 0.0d0) Then
          lambdaII(iplanet,i) = 0.5*mratio*mratio/(deltap)**4
!          lambdaII(iplanet,i) = 0.1*mratio*mratio/(deltap)**4.0d0
!          lambdaII(iplanet,i) = 0.1*mratio*G*mp(iplanet)/(deltap)**4.0d0
          if(rz(i) < ap(iplanet)) then
             lambdaII(iplanet,i) = -lambdaII(iplanet,i)*(rz(i))**4
!             lambdaII(iplanet,i) = -lambdaII(iplanet,i)*(ap(iplanet))**3.0d0 
          else
             lambdaII(iplanet,i) = lambdaII(iplanet,i)*(ap(iplanet))**4
!             lambdaII(iplanet,i) = lambdaII(iplanet,i)*(ap(iplanet))**3.0d0 
          endif
        endif

        if(soften) then
           softenfactor = 0.0001
           If (deltamax .gt. 0.0d0) Then
             softenfactor = abs(rz(i)-ap(iplanet))/deltamax
           EndIf
           if(softenfactor < 0.0001) softenfactor = 0.0001           
           lambdaII(iplanet,i)=lambdaII(iplanet,i)*softenfactor           
        endif

        !*************************************
        ! Compute the Type I specific torque (cf Nayakshin)
        !*************************************

        ! Must first compute the total torque and then normalise
        ! to obtain the correct Type I migration timescale
           
        ! Integrate the torque over all radii

        if ((H_d(i)+rhill) .gt. 0.0d0) Then 
          if(rz(i)<ap(iplanet)) then
!             typeInorm = typeInorm + exp(-deltap/(H_d(i)+rhill))*sigma_d(i)/drzm1(i)    
             typeInorm = typeInorm + 2.0d0*pi*rz(i)*exp(-deltap/(H_d(i)+rhill))*sigma_d(i)/drzm1(i)        
          else
!             typeInorm = typeInorm + exp(-deltap/(H_d(i)+rhill))*sigma_d(i)/drzm1(i)    
             typeInorm = typeInorm + 2.0d0*pi*rz(i)*exp(-deltap/(H_d(i)+rhill))*sigma_d(i)/drzm1(i)        
          endif
        endif  
       
     enddo

     ! Compute Type I migration timescale
     call calc_typeI_migration(iplanet,tmig1)

     ! Find normalisation constant to ensure correct Type I timescale

     if(tmig1 > 1.0e-30) then
!        lambda_dash = mp(iplanet)/(omega_d(iplanetrad(iplanet))*ap(iplanet)*4.0*pi*tmig1*typeInorm)
        lambda_dash = mp(iplanet)*dsqrt(G*mstar*ap(iplanet))/2.0d0/tmig1
      else
        lambda_dash = 0.0
      endif

     ! Now compute functional form of lambda 
     ! (assuming concentration of torque in planet local vicinity)
     do i=isr,ier

        deltap = abs(rz(i)-ap(iplanet))

        !if(deltap<H_d(i)) deltap =H_d(i)
        !if(deltap<rhill) deltap=rhill

        lambdaI(iplanet,i) = lambda_dash*exp(-deltap/(H_d(i)+rhill))/typeInorm
     
        lambdaI(iplanet,i) = lambdaI(iplanet,i)*rz(i)/(G*mstar)

        !if(rz(i)<ap(iplanet)) lambdaI(iplanet,i) = -lambdaI(iplanet,i)
       
     enddo


     !**************************************************
     ! Now compute the relative dominance of each torque 
     !**************************************************

     Pcrit = 0.75 * H_d(iplanetrad(iplanet))/rhill 
     Pcrit = Pcrit + 50.0*alpha_d(iplanetrad(iplanet))*(H_d(iplanetrad(iplanet))/ap(iplanet))**2/mratio

     gap_crit = 40.0*mstar*alpha_d(iplanetrad(iplanet))*(H_d(iplanetrad(iplanet))/ap(iplanet))**2

     fII(iplanet) = exp(-(Pcrit-1.0)**2.0d0)
     if (Pcrit .lt. 1.0d0) fII(iplanet) = 1.0d0
     if(fII(iplanet) > 1.0) fII(iplanet) = 1.0d0

!     if (t/yr .gt. 5.0d5) Then
!        print*, Pcrit, fII(iplanet), H_d(iplanet)/au, mp(iplanet)/1.898d30 
!        print*, ap(iplanet)/au, rhill/au
!        print*, alpha_d(iplanetrad(iplanet)), gap_crit/1.898d30
!     endif

!     fII(iplanet) = 1.0d0

      !********************************************************
      ! Compute the total effective planet torque at this radius
      !*********************************************************

     torquei(iplanet,:) = lambdaI(iplanet,:)*(1.0-fII(iplanet)) + lambdaII(iplanet,:)*fII(iplanet)

     ! Add planet contribution to the total torque exerted on the disc

     total_planet_torque(:) = total_planet_torque(:) + torquei(iplanet,:)
!     total_planet_torque(:) = total_planet_torque(:) + lambdaII(iplanet,:)*fII(iplanet)  

  enddo

  ! Ensure torque magnitude does not exceed maximum permitted value
  ! (Alexander & Armitage 2009, ApJ, 704, 989)

  do i=isr,ier
     aspectratio = H_d(i)/rz(i)

     if(abs(total_planet_torque(i)) > 0.1*aspectratio) then
        total_planet_torque(i) = 0.1*aspectratio*total_planet_torque(i)/abs(total_planet_torque(i))
     endif
   
  enddo

  ! Set brief time delay for planet torque activation

  if (t .lt. tdelay_planettorque*yr) then
     total_planet_torque(:) = total_planet_torque(:)*1.0e-30
     total_planet_torque(:) = 0.0d0
  endif

   torque_term(:) = 2.0*omega_d(:)*rf(:)*rf(:)*sigma_d(:)*total_planet_torque(:)
!   torque_term(:) = 2.0d0*omega_d(:)*dsqrt(rf(:))*sigma_d(:)*total_planet_torque(:) 

   ! Zero torque_term at planet locations

   do iplanet=1,nplanet
      torque_term(iplanetrad(iplanet)) = 0.0
      if (iplanetrad(iplanet) .lt. nrannuli) torque_term(iplanetrad(iplanet)+1) = 0.0
   enddo

end subroutine compute_planet_torques
