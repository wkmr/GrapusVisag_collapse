  SUBROUTINE midplane_properties_grav_fixedQ
! Routine calculates state properties of the self-gravitating disc midplane
! tauplus is passed in case there is a layer above

  use stardata
  use unitdata
  use eosdata

  implicit none
    
!    real,dimension(nrannuli) :: tauplus
  real(kind=8) :: twoDint,fine,oldtry
  real(kind=8) :: rho, Teff
  real :: heatfunc, T_dold, dT_d, gamma_d_old, kappa_d_old
  integer :: i

  IF(allocated(coolfunc)) deallocate(coolfunc)
  IF(allocated(Q)) deallocate(Q)
  IF(allocated(omega_d)) deallocate(omega_d)
  IF(allocated(cs_d)) deallocate(cs_d)
  IF(allocated(T_d)) deallocate(T_d)
  IF(allocated(H_d)) deallocate(H_d)
  IF(allocated(kappa_d)) deallocate(kappa_d)
  IF(allocated(tcool)) deallocate(tcool)
  IF(allocated(alpha_g)) deallocate(alpha_g)
  IF(allocated(nu_tc)) deallocate(nu_tc)

  allocate(coolfunc(nrannuli))
  allocate(Q(nrannuli))
  allocate(omega_d(nrannuli))
  allocate(cs_d(nrannuli))
  allocate(T_d(nrannuli))
  allocate(H_d(nrannuli))
  allocate(kappa_d(nrannuli))
  allocate(tcool(nrannuli))
  allocate(alpha_g(nrannuli))
  allocate(nu_tc(nrannuli))

  coolfunc(:) = 0.0

  do i = ier, isr, -1

     omega_d(i) = Sqrt(G*mstar/rz(i)**3.0d0)

     Q(i) = Qcrit
     cs_d(i) = pi*G*sigma_d(i)*Q(i)/omega_d(i)

     H_d(i) = cs_d(i)/omega_d(i)
    
!     print*, i, H_d(i), cs_d(i), omega_d(i)

     IF(H_d(i)>1.0e-40) THEN
        rho = 0.5d0*sigma_d(i)/H_d(i)
     ELSE
        rho = 0.0
     ENDIF
  
     !	Interpolate over rho,T to get cs,kappa, mu, gamma

!     print*, i, rz(i), cs_d(i), sigma_d(i), rho

     IF((rho>=1.0e-25).and.(cs_d(i).gt.50.0d0)) THEN
 
        call eos_cs(rho,cs_d(i))

        gamma_d(i) = gammamuT(1)
        T_d(i) = gammamuT(3)
        kappa_d(i) = gammamuT(4)

!        print*, i, rho, T_d(i), T_source(i), kappa_d(i), '  here' 

        if(T_d(i)<T_source(i)) then
           T_d(i) = T_source(i)
           call eos_T(rho,T_d(i))
!           Q(i) = cs_d(i)*omegaK(i)/(pi*G*sigma_d(i))
!           H_d = cs_d(i)/omegaK(i)
        endif
        
        ! Now calculate cooling time
        tau_d(i) = kappa_d(i)*sigma_d(i)
! + tauplus(i)

        ! Cooling function for disc
        coolfunc(i) = 16.0d0/3.0d0*stefan*(T_d(i)**4.0d0-T_source(i)**4.0d0)
!        coolfunc(i) = 16.0d0/3.0d0*stefan*(T_d(i)**4.0d0)
        coolfunc(i) = coolfunc(i)*tau_d(i)/(1.0d0+tau_d(i)**2.0d0)

        Teff = T_d(i)**4.0d0*tau_d(i)/(1.0d0+tau_d(i)**2.0d0)
        Teff = Teff**0.25d0

        twoDint = cs_d(i)**2.0*sigma_d(i)/gamma_d(i)/(gamma_d(i)-1.0d0)

        If (coolfunc(i).ne.0.0d0) Then
           tcool(i) = twoDint/coolfunc(i)
        Else
           tcool(i) = 1.0d35
        EndIf

        alpha_g(i) =  9.0D0/4.0D0*gamma_d(i)*(gamma_d(i)-1.0d0)*tcool(i)*omega_d(i)
        alpha_g(i) = 1.0d0/alpha_g(i)

!        print*, rz(i)/AU, alpha_g(i) 

        dT_d = 10000.0d0
        If (alpha_g(i) .lt. 5.0d-4) Then
          alpha_g(i) = 5.0d-4

          coolfunc(i) = 0.0d0
          heatfunc = 1.0d0

          T_dold = T_d(i)

          cs_d(i) = Sqrt(gamma_d(i)*k_B*T_d(i)/2.4d0/m_H)
          H_d(i) = cs_d(i)/omega_d(i)

          rho = sigma_d(i)/2.0d0/H_d(i)

          call eos_cs(rho,cs_d(i))

          gamma_d(i) = gammamuT(1)
          T_d(i) = gammamuT(3)
          kappa_d(i) = gammamuT(4)

          tau_d(i) = kappa_d(i)*sigma_d(i)

          heatfunc = 9.0d0*sigma_d(i)*alpha_g(i)*cs_d(i)**2.0d0*omega_d(i)
          T_d(i) = heatfunc*(1.0d0+tau_d(i)**2.0d0)/tau_d(i)
          T_d(i) = T_d(i) +16.0d0*stefan/3.0d0*T_source(i)**4.0d0
          T_d(i) = (T_d(i)*3.0d0/stefan/16.0d0)**0.25d0

          If (T_d(i) .lt. T_source(i)) Then
            T_d(i) = T_source(i)
          EndIf
          If (T_d(i) .gt. 1000.0d0) Then
            T_d(i) = 1000.0d0
          EndIf       
          cs_d(i) = Sqrt(gamma_d(i)*k_B*T_d(i)/2.4d0/m_H)
          H_d(i) = cs_d(i)/omega_d(i)

          rho = sigma_d(i)/2.0d0/H_d(i)

          do while (ABS((T_dold-T_d(i))/T_d(i)) .gt. 0.01d0) 

            T_dold = T_d(i)

            cs_d(i) = Sqrt(gamma_d(i)*k_B*T_d(i)/2.4d0/m_H)
            H_d(i) = cs_d(i)/omega_d(i)

            rho = sigma_d(i)/2.0d0/H_d(i)
          
            call eos_cs(rho,cs_d(i))
          
            gamma_d_old = gamma_d(i) 
            kappa_d_old = kappa_d(i)

            gamma_d(i) = gammamuT(1)
            T_d(i) = gammamuT(3)
            kappa_d(i) = gammamuT(4)

            kappa_d(i) = 0.01d0 
            gamma_d(i) = 1.66d0

!            If (kappa_d(i) .lt. 1.0d0) Then
!              kappa_d(i) = kappa_d_old
!              gamma_d(i) = gamma_d_old
!            EndIf

            If (ABS((kappa_d_old - kappa_d(i))/kappa_d(i)) .gt. 0.01) Then
              kappa_d(i) = kappa_d_old + 0.01*(kappa_d(i)-kappa_d_old)**2.0d0/ABS(kappa_d(i)-kappa_d_old)
            EndIf

            tau_d(i) = kappa_d(i)*sigma_d(i)

            heatfunc = 9.0d0*sigma_d(i)*alpha_g(i)*cs_d(i)**2.0d0*omega_d(i)
            T_d(i) = heatfunc*(1.0d0+tau_d(i)**2.0d0)/tau_d(i)
            T_d(i) = T_d(i) +16.0d0*stefan/3.0d0*T_source(i)**4.0d0
            T_d(i) = (T_d(i)*3.0d0/stefan/16.0d0)**0.25d0

            If (T_d(i) .lt. T_source(i)) Then
              T_d(i) = T_source(i)
            EndIf
!            If (T_d(i) .gt. 1000.0d0) Then
!              T_d(i) = 1000.0d0
!            EndIf
            cs_d(i) = Sqrt(gamma_d(i)*k_B*T_d(i)/2.4d0/m_H)
            H_d(i) = cs_d(i)/omega_d(i)

            rho = sigma_d(i)/2.0d0/H_d(i)

!            print*, i, heatfunc, sigma_d(i), T_d(i), T_dold, alpha_g(i)
!            print*, rz(i)/AU, gamma_d(i), kappa_d(i)
!            print*, ABS((T_dold - T_d(i))/T_d(i))*100.0, ABS((kappa_d_old - kappa_d(i))/kappa_d(i))*100.0d0
          enddo
        EndIf
     ELSE
        T_d(i) = 0.0
        cs_d(i) = 0.0
        kappa_d(i) = 0.0
!        mu(i) = 0.0
        gamma_d(i) = 0.0
        tcool(i) = 1.0d35
        coolfunc(i) = 0.0
     ENDIF


!    if(gamma_d(i)>1.000001) then
!        alpha_g(i) =  9.0d0/4.0d0*gamma_d(i)*(gamma_d(i)-1.0)*tcool(i)*omega_d(i)
!        alpha_g(i) =  1.0d0/alpha_g(i)
     
        !IF(t/yr> 999.0) print*, i, ' alpha is ', alpha_visc
!     else if (alpha_g(i) .lt. 1.0d-12.or.gamma_d(i)<1.00001) THEN
!        alpha_g(i) = 1.0d-12
!     endif

!    if (alpha_g(i) .lt. 5.0d-4) alpha_g(i) = 5.0d-4

!    heatfunc = 9.0d0*sigma_d(i)*alpha_g(i)*cs_d(i)**2.0d0*omega_d(i)
!    T_d(i) = heatfunc*(tau_d(i) + 1.0d0/tau_d(i))
!    T = T + 16.0d0*sigma_SB/3.0d0*T_irr1**4.0d0
!    T_d(i) = (T_d(i)*3.0d0/stefan/16.0d0)**0.25d0

!    if (alpha_g(i) .eq. 5.0d-4) cs_d(i) = SQRT(gamma_d(i)*k_B*T_d(i)/(2.4*m_H))

!     print*, i, rz(i)/udist, sigma_d(i), rho, alpha_g(i)
!     tcool(i), alpha_d(i), alpha_g(i)

     nu_tc(i) = alpha_g(i)*cs_d(i)*cs_d(i)/omega_d(i)               
     alpha_d(i) = alpha_g(i)

!     print*, i, T_d(i), alpha_g(i), alpha_d(i)
     
  enddo

return

END SUBROUTINE midplane_properties_grav_fixedQ
