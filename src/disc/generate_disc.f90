SUBROUTINE generate_disc
! Subroutine creates a self-gravitating disc, constrained by M_d and dot(M_d)
! We fix Q=1, and follow a similar procedure to Clarke (2009)

use stardata
use eosdata
use unitdata
use embryodata, only: fg
use clouddata
use winddata

implicit none

integer :: i, nzeros, ntries
real :: mtry, mdot_try,sigma_old,T,dT,fine
real :: beta

logical :: file_exists

mtry = 0.0
i = 0

! TODO - randomly generate fg, q_disc in some way

nfiles = tmax/t_disc_dump
nzeros = int(log10(nfiles)) + 2
write(zerostring, '(I1)') nzeros
snapshotformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"

print*, tmax, t_disc_dump, nfiles, nzeros, snapshotformat

call random_number(rnum)
fg = 0.02 + rnum*0.01

q_disc = 10.0
do while (q_disc .gt. q_disc1)
  call random_number(rnum)
  q_disc = q_disc0 + rnum*2.0d0*q_disc1
enddo
!q_disc = q_disc0

mdisc = q_disc*mstar

mdotvisc = 1.0d0
do while (mdotvisc .gt. mdotvisc1)
  call random_number(rnum)
  mdotvisc = mdotvisc0 + 2.0d0*(mdotvisc1-mdotvisc0)*rnum
enddo
mdotvisc = 10.0**(mdotvisc)*umass/yr

If (stell_irr .eq. 'n') Then
  Tirr = 1000.0d0
  do while (Tirr .gt. Tirr1)
    call random_number(rnum)
    Tirr = Tirr0 + rnum*2.0d0*Tirr1
  enddo
EndIf

alpha_visc = 1.0d0
do while (alpha_visc .gt. alpha_visc1)
  call random_number(rnum)
  alpha_visc = alpha_visc0 + rnum*2.0d0*alpha_visc1
enddo

print*, 'Generating disc for star ',istar, mstar/umass

sigma_old = 1.0e6
maxstep = 1.0d0

Lx = 10.0d0*Lx_1
do while (Lx .gt. Lx_1)
  call random_number(rnum)
  Lx = Log10(Lx_0) + 3.0d0*Log10(Lx_1)*rnum
  Lx = 10.0d0**Lx
enddo
print*, 'Luminosity of central star if 1 MSun: ', Lx 

If (runmode == 'C1') Then
  call setup_cloud
EndIf

print*, TRIM(prefix)//'.log'

inquire(file=TRIM(prefix)//'.log', exist=file_exists)

print*, file_exists

If (.not. file_exists) Then
  open(itime,file=TRIM(prefix)//'.log',access='append',status='unknown')
  write(itime,*) mstar/umass, mdisc/umass, q_disc, mdotvisc*yr/umass, Tirr, alpha_visc, Lx
  close(itime)
EndIf


IF(allocated(sigma_d)) deallocate(sigma_d)
IF(allocated(sigma_d_collapse)) deallocate(sigma_d_collapse)
IF(allocated(snew)) deallocate(snew)
IF(allocated(dsigma_cloud)) deallocate(dsigma_cloud)
IF(allocated(fM)) deallocate(fM)
IF(allocated(E_acc)) deallocate(E_acc)
IF(allocated(cs_d)) deallocate(cs_d)
IF(allocated(omega_d)) deallocate(omega_d)
IF(allocated(betac_d)) deallocate(betac_d)
IF(allocated(mjeans)) deallocate(mjeans)
IF(allocated(ljeans)) deallocate(ljeans)
IF(allocated(r_d)) deallocate(r_d)
IF(allocated(H_d)) deallocate(H_d)
IF(allocated(alpha_d)) deallocate(alpha_d)
IF(allocated(T_d)) deallocate(T_d)
IF(allocated(Tnew)) deallocate(Tnew)
IF(allocated(kappa_d)) deallocate(kappa_d)
IF(allocated(gamma_d)) deallocate(gamma_d)
IF(allocated(tau_d)) deallocate(tau_d)
IF(allocated(gamma_j)) deallocate(gamma_J)
IF(allocated(rz)) deallocate(rz)
IF(allocated(rzm1)) deallocate(rzm1)
IF(allocated(rz1_2)) deallocate(rz1_2)
IF(allocated(drzm1)) deallocate(drzm1)
IF(allocated(rf)) deallocate(rf)
If(allocated(rf1_2)) deallocate(rf1_2)
IF(allocated(drfm1)) deallocate(drfm1)
IF(allocated(T_source)) deallocate(T_source)
IF(allocated(heatfunc)) deallocate(heatfunc)
IF(allocated(Q)) deallocate(Q)

allocate(sigma_d(nrannuli))
allocate(sigma_d_collapse(nrannuli))
allocate(snew(nrannuli))
allocate(dsigma_cloud(nrannuli))
allocate(fM(nrannuli))
allocate(E_acc(nrannuli))
allocate(cs_d(nrannuli))
allocate(omega_d(nrannuli))
allocate(betac_d(nrannuli))
allocate(mjeans(nrannuli))
allocate(ljeans(nrannuli))
allocate(gamma_j(nrannuli))
allocate(r_d(nrannuli))
allocate(H_d(nrannuli))
allocate(alpha_d(nrannuli))
allocate(T_d(nrannuli))
allocate(Tnew(nrannuli))
allocate(kappa_d(nrannuli))
allocate(gamma_d(nrannuli))
allocate(tau_d(nrannuli))

allocate(T_source(nrannuli))
allocate(heatfunc(nrannuli))
allocate(Q(nrannuli))

! Set up grid arrays
allocate(rz(nrannuli))
allocate(rzm1(nrannuli))
allocate(rz1_2(nrannuli))
allocate(drzm1(nrannuli))

allocate(rf(nrannuli))
allocate(rf1_2(nrannuli))
allocate(drfm1(nrannuli))

rout = 500.0d0*udist

beta = (rout/rin)**(1.0d0/dble(nrannuli))

isr = 2
ier = nrannuli-1

rf(isr) = rin
rz(isr) = rin * sqrt(beta)
rzm1(isr)  = 1.0d0 / rz(isr)
rf1_2(isr) = sqrt(rf(isr))
rz1_2(isr) = sqrt(rz(isr))

rf(1) = rf(isr) / beta
rz(1) = rz(isr) / beta
rzm1(1)  = 1.0d0 / rz(1)
rf1_2(1) = sqrt(rf(1))
rz1_2(1) = sqrt(rz(1))

do i = isr+1,ier+1
   rf(i) = rf(i-1) * beta
   rz(i) = rz(i-1) * beta
   rzm1(i)  = 1.0d0 / rz(i)
   rf1_2(i) = sqrt(rf(i))
   rz1_2(i) = sqrt(rz(i))
enddo

do i = isr, ier+1
   drfm1(i) = 1.0d0 / (rz(i)-rz(i-1))
enddo

do i = isr-1, ier
   drzm1(i) = 1.0d0 / (rf(i+1)-rf(i))
enddo

do i = 1, nrannuli
  sigma_d(i) = 0.0d0
  cs_d(i) = 0.0d0
  betac_d(i) = 0.0d0
  mjeans(i) = 0.0d0
  alpha_d(i) = alpha_visc 
  r_d(i) = rz(i)

  if (stell_irr == 'y') then 
    T_source(i) = (Lstar*3.8d33/4.0d0/pi/stefan/r_d(i)**2.0d0)**0.25d0
  else
    T_source(i) = Tirr
  endif
  T_d(i) = T_source(i)
enddo

sigma_d_collapse(:) = 0.0d0
dsigma_cloud(:) = 0.0d0
fM(:) = 0.0d0
E_acc(:) = 0.0d0

i = 0
mtry = 0.0d0

If (runmode .ne. 'C1') then

  DO WHILE(mtry < mdisc .and. i<nrannuli-1)

    i = i+1
    r_d(i) = rz(i)
    rmax = r_d(i)
    dr = 1.0/drzm1(i)

    omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

    sigma_d(i) = sigma_old

    dT = 1.0e30
    ntries = 0
    fine = 0.01
    DO WHILE(ABS(dT) > tolerance)

!	Calculate sound speed assuming fixed Q

      cs_d(i) = Qcrit*pi*G*sigma_d(i)/(omega_d(i))

!	Calculate scale height

      H_d(i) = cs_d(i)/omega_d(i)

      rhomid = sigma_d(i)/(2.0*H_d(i))

!	Use EoS to calculate tau, gamma, betac

      CALL eos_cs(rhomid, cs_d(i))

      gamma_d(i) = gammamuT(1)
      T_d(i) = gammamuT(3)
      kappa_d(i) = gammamuT(4)

      tau_d(i) = sigma_d(i)*kappa_d(i)
      betac_d(i) = (tau_d(i)+1.0/tau_d(i))*cs_d(i)*cs_d(i)*omega_d(i)*sigma_d(i)*3.0d0/&
                   (16.0d0*stefan*(T_d(i)**4.0-T_source(i)**4.0d0)*gamma_d(i)*(gamma_d(i)-1.0))

!	Calculate alpha from this value --> accretion rate

      alpha_d(i) = 4.0/(9.0*gamma_d(i)*(gamma_d(i)-1.0)*betac_d(i))

!	Compare with imposed accretion rate

      if (alpha_d(i) .lt. alpha_visc) Then
        alpha_d(i) = alpha_visc

        heatfunc(i) = 9.0d0*sigma_d(i)*alpha_d(i)*cs_d(i)**2.0d0*omega_d(i)/4.0d0
        T = heatfunc(i)*(tau_d(i) + 1.0d0/tau_d(i))
        T = T + 16.0d0*stefan/3.0d0*T_source(i)**4.0d0
        T = (T*3.0d0/stefan/16.0d0)**0.25d0   

        T_d(i) = T
      endif

      call eos_T(rhomid,T_d(i))
      gamma_d(i) = gammamuT(1)
      kappa_d(i) = gammamuT(4)
      cs_d(i) = SQRT(gammamuT(1)*k_B*T_d(i)/(gammamuT(2)*m_H))

      Q(i) = cs_d(i)*omega_d(i)/pi/G/sigma_d(i)

      mdot_try = 3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i)*sigma_d(i)/omega_d(i)

      dT = (mdotvisc-mdot_try)/mdotvisc

      sigma_old = sigma_d(i)
      sigma_d(i) = sigma_d(i)*(1.0 + dT/(abs(dT))*fine)

      ntries = ntries + 1

      IF (ntries>500) THEN
        fine = fine/10.0
        if (fine .lt. 1.0d-10) exit 
        ntries = 0
      ENDIF

    ENDDO

!	Check for MRI activation
!	If so, then set alpha=0.01 and readjust sigma to maintain mdotvisc

!IF(T_d(i)>1000.0) THEN
!	print*, 'MRI active ',r/udist, T, alpha, sigma
!alpha_d(i) = 0.01
!sigma_d(i) = mdotvisc*omega_d(i)/(3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i))
!	print*, 'MRI active ', T, alpha, sigma
!ENDIF

!  MJEANS(i) = 4.0*1.412*pi*pi*pi/(3.0*G)
    mjeans(i) = Sqrt(3.0)*pi*pi*pi*Sqrt(Qcrit)/(32.0*G) 
    mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*Sqrt(1.0+4.47*Sqrt(alpha_d(i))))

! old mjeans
!    mjeans(i) = 4.0*sqrt(2.0)*pi*pi*pi*sqrt(Qcrit)/(3.0*G)
!    mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*(1.0+4.47*sqrt(alpha_d(i))))

!    mjeans(i) = pi*cs_d(i)*cs_d(i)*cs_d(i)/6.0d0/G**(3.0d0/2.0d0)
!    mjeans(i) = mjeans(i)/(sigma_d(i)/(2.0d0*H_d(i)))**0.5d0

    ljeans(i) = cs_d(i)*sqrt(3.0d0*pi/32.0d0/G/rhomid)

    IF(alpha_d(i) .lt. alpha_frag) THEN
       gamma_j(i) = 2.0*(9.0*alpha_d(i)*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i)) 
    ELSE
       gamma_j(i) = 2.0*(9.0*alpha_frag*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i))
    ENDIF

    If (gamma_j(i) .ne. 0.0d0) Then
       gamma_j(i) = 1.0/gamma_j(i)
    Else 
       gamma_j(i) = 10.0d0
    EndIf   

!	Calculate enclosed Mass and mass ratio
    If (r_d(i) .gt. 5.0d0*udist) Then
      mtry = mtry + twopi*r_d(i)*sigma_d(i)*dr
    EndIf

    IF(i==nrannuli) THEN
      print*, 'ERROR: nrannuli insufficient'
      print*, 'Either recompile or increase dr in .params file'
      STOP
    ENDIF
  ENDDO

EndIf

If (runmode .ne. 'C1') then
  rout = r_d(i)
  irout = i

  print*, 'Mass Ratio = ', q_disc
  print*, 'Disc Mass = ', q_disc*mstar/umass, mtry/umass
  print*, 'Accretion Rate = ', mdotvisc*yr/umass, mdot_try*yr/umass
  print*, 'Rout = ', rout/AU
  print*, 'fg = ', fg
  print*, 'Viscous alpha = ', alpha_visc
  print*, 'Tirr = ', Tirr

  print*, sigma_d(irout), sigma_d(irout+1), irout

  T_d(ier+1) = T_d(ier)
endif

!
! The next do loop allows you to set up a prescribed disc profile for testing
!

!do i = 1, nrannuli
!  If (r_d(i) .le. 200.0d0*au) then
!    gamma_d(i) = 1.0d0 
!    sigma_d(i) = 0.1d0*Mstar/(2.0d0*pi*au**1.5*2.0*((35.0d0*au)**0.5-au**0.5))
!    sigma_d(i) = sigma_d(i)*(r_d(i)/au)**(-1.5d0)
!    sigma_d(i) = 0.1d0*Mstar/(2.0d0*pi*au*(100.0d0*au - au))
!    sigma_d(i) = sigma_d(i)*(r_d(i)/au)**(-1.0d0)
!    sigma_d(i) = 1.0d3*(r_d(i)/au)**(-1.0d0)
!    T_d(i) = 250.0d0*(r_d(i)/1.496d13)**(-0.5d0)
!    cs_d(i) = Sqrt(gamma_d(i)*k_B*T_d(i)/2.4d0/m_H)
!    alpha_d(i) = 0.01d0
!  else
!    gamma_d(i) = 0.0d0
!    sigma_d(i) = 0.0d0
!    T_d(i) = 0.0d0
!    cs_d(i) = 0.0d0
!    alpha_d(i) = 0.0d0
!  endif
!  write(25,*) i, r_d(i), sigma_d(i), cs_d(i), T_d(i), tau_d(i), kappa_d(i), alpha_d(i)
!enddo

If (runmode .ne. 'C1') Then
  mtry = 0.0d0
  do i = 1, nrannuli-1
    mtry = mtry + 2.0d0*pi*r_d(i)*sigma_d(i)*(r_d(i+1) - r_d(i))
  enddo

  If (mstar .gt. 0.0d0) then
    print*, mtry/mstar,'   disc-to-star mass ratio'
  EndIf
EndIf

END SUBROUTINE generate_disc
