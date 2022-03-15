SUBROUTINE generate_disc
! Subroutine creates a self-gravitating disc, constrained by M_d and dot(M_d)
! We fix Q=1, and follow a similar procedure to Clarke (2009)

use stardata
use eosdata
use unitdata
use embryodata, only: fg

implicit none

integer :: i, nzeros
real :: mtry, mdot_try,sigma_old,heatfunc,T,dT,fine
real :: beta

mtry = 0.0
i=0

! TODO - randomly generate fg, q_disc in some way

nfiles = tmax/t_disc_dump
nzeros = int(log10(nfiles)) + 2
write(zerostring, '(I1)') nzeros
snapshotformat = "(I"//TRIM(zerostring)//"."//TRIM(zerostring)//")"

print*, tmax, t_disc_dump, nfiles, nzeros, snapshotformat

fg = 0.01 + ran2(iseed)*0.01
q_disc = 0.5 + ran2(iseed)*0.25
!q_disc = 1.0d0

mdisc = q_disc*mstar
mdotvisc = -6.0 + 2.0*ran2(iseed)
!mdotvisc = -6.5
mdotvisc = 10.0**(mdotvisc)*umass/yr

print*, 'Generating disc for star ',istar,': ', fg,mstar, q_disc, mdotvisc*yr/umass

print*, mdisc/mdotvisc, 5.0*2.0*pi/Sqrt(G*mstar/(100.0*AU)**3.0), mdisc/mdotvisc/(5.0*2.0*pi/Sqrt(G*mstar/(100.0*AU)**3.0))

sigma_old = 1.0e6
maxstep = 1.0d0

IF(allocated(sigma_d)) deallocate(sigma_d)
IF(allocated(snew)) deallocate(snew)
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

allocate(sigma_d(nrannuli))
allocate(snew(nrannuli))
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

! Set up grid arrays
allocate(rz(nrannuli))
allocate(rzm1(nrannuli))
allocate(rz1_2(nrannuli))
allocate(drzm1(nrannuli))

allocate(rf(nrannuli))
allocate(rf1_2(nrannuli))
allocate(drfm1(nrannuli))

rout = 300.0d0*udist

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

!do i = 1, ier+1
!  print*, rf(i)/udist, rz(i)/udist
!enddo

do i = 1, nrannuli
  sigma_d(i) = 0.0d0
  cs_d(i) = 0.0d0
  betac_d(i) = 0.0d0
  mjeans(i) = 0.0d0
  alpha_d(i) = 0.0d0
  r_d(i) = rz(i)

  T_source(i) = 10.0d0  
enddo

i = 0
mtry = 0.0d0

DO WHILE(mtry < mdisc .and. i<nrannuli-1)
!DO WHILE(i < nrannuli)

  i = i+1
!  r_d(i) = rin + (i-1)*dr
  r_d(i) = rz(i)
  rmax = r_d(i)
  dr = 1.0/drzm1(i)

!  print*, i, r_d(i)

  omega_d(i) = sqrt(G*mstar/(r_d(i)*r_d(i)*r_d(i)))

!  omega_d(i) = sqrt(G*mstar/(rz(i)*rz(i)*rz(i)))

! Iterate to find sigma

  sigma_d(i) = sigma_old

  dT = 1.0e30
  ntries = 0
  fine = 0.01
  DO WHILE(ABS(dT)> tolerance)
	
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
                 (16.0d0*stefan*(T_d(i)**4.0-20.0d0**4.0d0)*gamma_d(i)*(gamma_d(i)-1.0))

!	Calculate alpha from this value --> accretion rate

    alpha_d(i) = 4.0/(9.0*gamma_d(i)*(gamma_d(i)-1.0)*betac_d(i))

!    print*, r_d(i)/1.5d13, betac_d(i), alpha_d(i), sigma_d(i)

!	print*, alpha_d(i), gamma_d(i),betac_d(i),cs_d(i),omega_d(i)		
!	Compare with imposed accretion rate

    if (alpha_d(i) .lt. 5.0d-4) alpha_d(i) = 5.0d-4

    heatfunc = 9.0d0*sigma_d(i)*alpha_d(i)*cs_d(i)**2.0d0*omega_d(i)
    T = heatfunc*(tau_d(i) + 1.0d0/tau_d(i))
    T = T + 16.0d0*stefan/3.0d0*20.0d0**4.0d0
    T = (T*3.0d0/stefan/16.0d0)**0.25d0

!    print*,i, r_d(i)/AU, T     

    if (alpha_d(i) .eq. 5.0d-4) cs_d(i) = SQRT(gamma_d(i)*k_B*T/(2.4*m_H))

    mdot_try = 3.0*pi*alpha_d(i)*cs_d(i)*cs_d(i)*sigma_d(i)/omega_d(i)

    dT = (mdotvisc-mdot_try)/mdotvisc

    sigma_old = sigma_d(i)
!	sigma_d(i) = sigma_d(i)*(1.0 + dT*0.01)
    sigma_d(i) = sigma_d(i)*(1.0 +dT/(abs(dT))*fine)		
	
!	print*,i, sigma_d(i), sigma_old, dT, mdot_try, mdotvisc
    ntries = ntries + 1

    IF (ntries>500) THEN
      fine = fine/10.0
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

  sigma_d(i+1) = sigma_old
!  MJEANS(i) = 4.0*1.412*pi*pi*pi/(3.0*G)
  mjeans(i) = Sqrt(3.0)*pi*pi*pi*Sqrt(Qcrit)/(32.0*G) 

!  mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*Sqrt(1.0+1.0/sqrt(betac_d(i))))
  mjeans(i) = mjeans(i)*cs_d(i)*cs_d(i)*cs_d(i)/(omega_d(i)*Sqrt(1.0+4.47*Sqrt(alpha_d(i))))
!  mjeans(i) = (pi/6.0)*cs_d(i)*cs_d(i)*cs_d(i)/G**1.5d0/Sqrt(rhomid)

  ljeans(i) = cs_d(i)*sqrt(3.0d0*pi/32.0d0/G/rhomid)

  IF(alpha_d(i)<0.1) THEN
     gamma_j(i) = 2.0*(9.0*alpha_d(i)*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i)) 
  ELSE
     gamma_j(i) = 2.0*(9.0*0.1*gamma_d(i)*(gamma_d(i)-1.0)/4.0 - 1/betac_d(i))
  ENDIF

  If (gamma_j(i) .ne. 0.0d0) Then
     gamma_j(i) = 1.0/gamma_j(i)
  Else 
     gamma_j(i) = 10.0d0
  EndIf   

!	Calculate enclosed Mass and mass ratio
  If (r_d(i) .gt. 5.0d0*udist) Then
    mtry = mtry + twopi*r_d(i)*sigma_d(i)*dr

!    print*, i, r_d(i)/udist, dr/udist, mtry/1.989d33, mdisc/1.989d33
  EndIf  

!  print*, mtry, mdisc

  IF(i==nrannuli) THEN
    print*, 'ERROR: nrannuli insufficient'
    print*, 'Either recompile or increase dr in .params file'
    STOP
  ENDIF
ENDDO

rout = r_d(i)
irout =i

do i = 1, nrannuli
  write(25,*) i, r_d(i), sigma_d(i), cs_d(i), T_d(i), tau_d(i), kappa_d(i), alpha_d(i)
!  print*, i, r_d(i)/AU, sigma_d(i), T_d(i), gamma_j(i)
enddo

END SUBROUTINE generate_disc
