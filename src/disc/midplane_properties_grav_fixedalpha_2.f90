  SUBROUTINE midplane_properties_grav_fixedalpha_2(t)
! Routine calculates state properties of the self-gravitating disc midplane
! tauplus is passed in case there is a layer above

  use stardata
  use unitdata
  use eosdata
  use clouddata
  use embryodata 

  implicit none
    
  real(kind=8) :: twoDint,fine,oldtry
  real(kind=8) :: rho, Teff
  real :: diskmass, collapse_term
  real :: gamma_d_old, kappa_d_old
  real :: t
  integer :: i
  character(1) :: gen_embs

  IF(allocated(coolfunc)) deallocate(coolfunc)
  IF(allocated(heatfunc)) deallocate(heatfunc)
  IF(allocated(Q)) deallocate(Q)
  IF(allocated(omega_d)) deallocate(omega_d)
  IF(allocated(cs_d)) deallocate(cs_d)
  IF(allocated(H_d)) deallocate(H_d)
  IF(allocated(kappa_d)) deallocate(kappa_d)
  IF(allocated(tcool)) deallocate(tcool)
  IF(allocated(alpha_g)) deallocate(alpha_g)
  IF(allocated(alpha_d)) deallocate(alpha_d)
  IF(allocated(nu_tc)) deallocate(nu_tc)
  IF(allocated(mjeans)) deallocate(mjeans)
  IF(allocated(ljeans)) deallocate(ljeans)

  allocate(coolfunc(nrannuli))
  allocate(heatfunc(nrannuli))
  allocate(Q(nrannuli))
  allocate(omega_d(nrannuli))
  allocate(cs_d(nrannuli))
  allocate(H_d(nrannuli))
  allocate(kappa_d(nrannuli))
  allocate(tcool(nrannuli))
  allocate(alpha_g(nrannuli))
  allocate(alpha_d(nrannuli))
  allocate(nu_tc(nrannuli))
  allocate(mjeans(nrannuli)) 
  allocate(ljeans(nrannuli))

  coolfunc(:) = 0.0

  diskmass = 0.0d0

  call sigma_mdot(t)
!  call Eacc_calc(t)  

  omega_d(:) = Sqrt(G*mstar/rz(:)**3.0d0)
  alpha_d(:) = alpha_visc

  gen_embs = 'N'

  do i = isr, ier

     diskmass = diskmass + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_d(i)

     omega_d(i) = Sqrt(G*(mstar+diskmass)/rz(i)**3.0d0)

     Q(i) = Qcrit
     if (omega_d(i) .gt. 0.0d0) then
       cs_d(i) = pi*G*sigma_d(i)*Q(i)/omega_d(i)
       H_d(i) = cs_d(i)/omega_d(i)
     else
       cs_d(i) = 0.0d0
       H_d(i) = 0.0d0
     endif
  
     IF(H_d(i)>1.0e-40) THEN
        rho = 0.5d0*sigma_d(i)/H_d(i)
     ELSE
        rho = 0.0
     ENDIF
  
     !	Interpolate over rho,T to get cs,kappa, mu, gamma

     IF((rho>=1.0e-25).and.(cs_d(i).gt.50.0d0)) THEN
 
        call eos_cs(rho,cs_d(i))

        gamma_d(i) = gammamuT(1)
        T_d(i) = gammamuT(3)
        kappa_d(i) = gammamuT(4)

        if(T_d(i)<T_source(i)) then
           T_d(i) = T_source(i)
           call eos_T(rho,T_d(i))
           Q(i) = cs_d(i)*omega_d(i)/(pi*G*sigma_d(i))
           H_d = cs_d(i)/omega_d(i)
        endif
        
        ! Now calculate cooling time
        tau_d(i) = kappa_d(i)*sigma_d(i)

        ! Cooling function for disc
        coolfunc(i) = 16.0d0/3.0d0*stefan*(T_d(i)**4.0d0-T_source(i)**4.0d0)
        coolfunc(i) = coolfunc(i)*tau_d(i)/(1.0d0+tau_d(i)**2.0d0)

        Teff = T_d(i)**4.0d0*tau_d(i)/(1.0d0+tau_d(i)**2.0d0)
        Teff = Teff**0.25d0

        twoDint = cs_d(i)**2.0*sigma_d(i)/gamma_d(i)/(gamma_d(i)-1.0d0)

        If (coolfunc(i).ne.0.0d0) Then
           tcool(i) = twoDint/coolfunc(i)
        Else
           tcool(i) = 1.0d35
        EndIf

     ELSE
        T_d(i) = T_source(i) 
        cs_d(i) = 0.0
        kappa_d(i) = 0.0
        gamma_d(i) = 0.0
        tcool(i) = 1.0d35
        coolfunc(i) = 0.0
     ENDIF

     alpha_g(i) = 0.0d0
     if(gamma_d(i)>1.000001) then
        alpha_g(i) = 9.0d0/4.0d0*gamma_d(i)*(gamma_d(i)-1.0)*tcool(i)*omega_d(i)
        If (alpha_g(i) .gt. 0.0d0) Then
          alpha_g(i) =  1.0d0/alpha_g(i)
        Else
          alpha_g(i) = 0.0d0
        EndIf

!        collapse_term = 4.0d0*E_acc(i)*dsigma_cloud(i)/9.0d0/cs_d(i)**2.0d0/sigma_d(i)/omega_d(i)

!        alpha_g(i) = alpha_g(i) - collapse_term
     else if (alpha_g(i) .lt. 1.0d-12 .or.gamma_d(i)<1.00001) THEN
        alpha_g(i) = 1.0d-12
     endif

     if (alpha_g(i) .lt. alpha_visc) then
       alpha_g(i) = alpha_visc
     endif

  enddo

  do i =isr, ier
    if (i .lt. 3) then
      alpha_g(i) = alpha_g(i) + alpha_g(i+1) + alpha_g(i+2) + alpha_g(i+3)
      alpha_g(i) = (alpha_g(i) + alpha_g(i+4))/5.0d0
    endIf
    if (i .gt. ier - 3) then
      alpha_g(i) = alpha_g(i) + alpha_g(i-4) + alpha_g(i-3) + alpha_g(i-2)
      alpha_g(i) = (alpha_g(i) - alpha_g(i-1))/5.0d0
    endIf
    if ((i .gt. 2)  .and. (i .lt. ier-2)) then
      alpha_g(i) = alpha_g(i-2) + alpha_g(i-1) + alpha_g(i)
      alpha_g(i) = (alpha_g(i) + alpha_g(i+1) + alpha_g(i+2))/5.0
    endif

    If (alpha_g(i) .gt. alpha_frag) Then
      If ((nembryo .eq. 0) .and. (rz(i) .gt. 5.0d0*au)) Then
        gen_embs = 'Y'
      EndIf
    EndIf
  enddo

  where (omega_d > 0.0d0)
     nu_tc(:) = alpha_g(:)*cs_d(:)*cs_d(:)/omega_d(:)
  elsewhere
     nu_tc(:) = 0.0d0
  end where

  alpha_g(1) = alpha_g(2)
  alpha_d(:) = alpha_g(:)

  cs_d(1) = cs_d(2)
  if ((runmode == 'C') .and. (nembryo == 0) .and. (mstar .gt. 0.0d0)) Then
 
    where (omega_d > 0.0d0)
      mjeans(:) = Sqrt(3.0)*pi*pi*pi*Sqrt(Qcrit)/(32.0d0*G)
      mjeans(:) = mjeans(:)*cs_d(:)*cs_d(:)*cs_d(:)
      mjeans(:) = mjeans(:)/(omega_d(:)*Sqrt(1.0+4.47*Sqrt(alpha_d(:))))
    elsewhere
      mjeans(:) = 0.0d0
    end where

    where ((sigma_d > 0.0d0) .and. (H_d > 0.0d0)) 
       ljeans(:) = cs_d(:)*Sqrt(3.0d0*pi/32.0d0/G/(sigma_d(:)/(2.0d0*H_d(:))))
    elsewhere
       ljeans(:) = 0.0d0
    end where 

!    If (gen_embs == 'Y') Then
!      print*, 'calling generate embryos'
  
!      call generate_embryos
 
!      If (nembryo .gt. 0) Then 
!        print*, 'setting up planets'

!        call setup_planets
!      EndIf 
!    EndIf    
  endif

return

END SUBROUTINE midplane_properties_grav_fixedalpha_2
