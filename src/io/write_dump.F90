
!-----------------------------------------------------------
! Subroutine to write radial dumps 
!-----------------------------------------------------------

subroutine write_dump(t)
#ifdef _OPENMP
  use OMP_LIB
#endif
  use stardata
!  use magdata
  use planetdata
  use embryodata
  use unitdata
  use winddata

  implicit none

  real(kind=8) :: sig_max,mdisk, tot_lumin,tcoolmin
  real :: t
  real :: mdot_grav, mdot_mag, grav_max,mag_max, mmag,mgrav
  real :: mdisk10, mdisk100
  integer :: i,ifirst, iplanet
  integer :: iprofomp, itimeomp, itorqueomp

  write (*,103) 'Output at time t = ',t/yr, '  years, dt = ',dt/yr, ' years'
102 format (12(1PE15.5))
103 format (A,1PE11.3,A,1PE11.3,A)
104 format (4(1PE15.5))
110 format (10(1PE15.5))
111 format (8(1PE15.5))

!if using openmp, make each file id unique to each thread
#ifdef _OPENMP
  iprofomp = iprof + OMP_GET_THREAD_NUM()
  itimeomp = itime + OMP_GET_THREAD_NUM()
  itorqueomp = itorque + OMP_GET_THREAD_NUM()
#else
  iprofomp = iprof
  itimeomp = itime + 1
  itorqueomp = itorque
#endif

  ! Calculate disc properties and spectrum
  call disc_properties(t)
!  call luminosity(tot_lumin,spectrum,Tc,tau)

  write(fileno, snapshotformat) snapshotcounter

  fileno = TRIM(fileno)
! write out disc profiles
  OPEN(iprofomp, file=TRIM(prefix)//'_profile.'//fileno, status='unknown')
  write(iprofomp,*) t/yr, nrannuli
  do i = isr, ier
     write (iprofomp,102) rz(i)/AU,sigma_d(i), cs_d(i), kappa_d(i), gamma_d(i),mu,T_d(i), &
    tau_d(i), nu_tc(i), alpha_d(i),Q(i), 0.0d0
!, torque_term(i)
  enddo
  close(iprofomp)

!  OPEN(ispec, file=TRIM(prefix)//'_spectrum.'//fileno, status='unknown')
!  write(ispec,*) 1000
!  do i = 1, 1000
!     write(ispec,*) spectrum(1,i), spectrum(2,i)
!  enddo

!  close(ispec)

!if (layerchoice=='y') then
!OPEN(iprof,file=TRIM(prefix)//'layer.'//fileno,status='unknown')
!WRITE(iprof,*) t/yr,nrgrid
!do i=isr,ier
!write(iprof,111) rz(i)/AU, sigma(i), sigma_m(i), sigma_tot(i), cs_m(i),&
!    kappa_m(i),gamma_m(i), mu_m(i),tau_m(i),nu_m(i),alpha_m
!enddo
!close(iprof)
!endif

if(nplanet .gt. 0) then
  open(iprofomp, file=TRIM(prefix)//'_planets.'//fileno,status='unknown')
  write(iprofomp,*)t/yr,nplanet,nactive
  do iplanet=1,nplanet
     write(iprofomp,*) alive(iplanet),mp(iplanet)/mjup, &
         ap(iplanet)/AU, tmig(iplanet)/yr, embryo(iplanet)%mcore/mjup, &
         embryo(iplanet)%r/rjup, embryo(iplanet)%rcore/rjup
  enddo
  close(iprofomp)

  open(itorqueomp,file=TRIM(prefix)//'_torque.'//fileno,status='unknown')
  do i=isr,ier
     write(itorqueomp,*) rz(i)/AU, 0.1*H_d(i)/rz(i), lambdaI(1,i), lambdaII(1,i), &
          total_planet_torque(i), torque_term(i), total_planet_torque(i+1)-total_planet_torque(i)
  enddo
  close(itorqueomp)

endif

if((nplanet .eq. 0).and.(runmode .eq. 'C')) then
  open(iprofomp, file=TRIM(prefix)//'_planets.'//fileno,status='unknown')
  write(iprofomp,*)t/yr,nplanet,nactive
  write(iprofomp,*) 0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0
  close(iprofomp)

  open(itorqueomp,file=TRIM(prefix)//'_torque.'//fileno,status='unknown')
  do i=isr,ier
     write(itorqueomp,*) rz(i)/AU, 0.1*H_d(i)/rz(i), 0.0, 0.0, 0.0, 0.0, 0.0
  enddo
  close(itorqueomp)
endif

!if(nplanet .eq. 0) then
!open(iprofomp, file=TRIM(prefix)//'_planets.'//fileno,status='unknown')
!write(iprofomp,*)t/yr,nplanet,nactive
!do iplanet=1,nplanet
!
!    write(iprof,*) alive(iplanet),mp(iplanet)/mjup, &
!         ap(iplanet)/AU, tmig(iplanet)/yr
!enddo
!close(iprofomp)
!endif


! Compute disk mass and maximum surface density
! Also compute radially averaged accretion rates

mdisk = 0.0d0
mdisk10 = 0.0d0
mdisk100 = 0.0d0
mgrav= 0.0d0
mmag = 0.0d0

grav_max = 0.0d0
mag_max = 0.0d0
sig_max = 0.0d0

mdot_wind = 0.0d0

do i = isr, ier
!mmag = mmag + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_m(i)
!mgrav = mgrav + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma(i)
mdisk = mdisk+2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_d(i)

if (rz(i) .le. 10.0d0*AU) then
 mdisk10 = mdisk10 + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_d(i)
endif
if (rz(i) .le. 100.0d0*AU) then
 mdisk100 = mdisk100 + 2.0d0*pi*rz(i)*(rf(i+1)-rf(i))*sigma_d(i)
endif

!grav_max = max(grav_max,sigma(i))
!mag_max = max(mag_max,sigma_m(i))
sig_max = max(sig_max,sigma_d(i))

!mdot_grav = mdot_grav + 3.0*pi*nu_tc(i)*sigma(i)/REAL(nrgrid)
!mdot_mag = mdot_mag + 3.0*pi*nu_m(i)*sigma_m(i)/REAL(nrgrid)

mdot_wind = mdot_wind + sigdot_wind(i)*twopi*rz(i)/drzm1(i)

enddo

mdisk = mdisk/solarmass
mdisk10 = mdisk10/solarmass
mdisk100 = mdisk100/solarmass
!mgrav = mgrav/solarmass
!mmag = mmag/solarmass

!mdot_grav = mdot_grav/msolyr
!mdot_mag = mdot_mag/msolyr
mdot_wind = mdot_wind/msolyr

! Write out snapshot data (disk mass, sigma, total luminosity)
!open(itimeomp,file=TRIM(prefix)//'_disc.'//fileno,status='unknown')

#ifdef _OPENMP
open(itimeomp,file=TRIM(prefix)//'.log',access='append',status='unknown')
write(itimeomp,111) t/yr, dt/yr, mstar/solarmass, mdisk10, mdisk100, mdisk, sig_max, mdot_wind
close(itimeomp)
#else
open(itime,file=TRIM(prefix)//'.log',access='append',status='unknown')
write(itime,111) t/yr, dt/yr, mstar/solarmass, mdisk10, mdisk100, mdisk, sig_max, mdot_wind
print*, t/yr, dt/yr, mstar/solarmass, mdisk10, mdisk100, mdisk, sig_max, mdot_wind
print*,' here log:', t/yr, itime
close(itime)
#endif

return

end  subroutine write_dump
