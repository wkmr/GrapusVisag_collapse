subroutine setup_cloud
  
  use stardata
  use clouddata
  use unitdata
  use winddata

  real :: Mcloud, Rcloud, f_cloud

  logical :: file_exists

  Mcloud = (Mcloud_0 + Mcloud_0/4.0d0 - 0.5d0*Mcloud_0*ran2(iseed))
  Mcloud = Mcloud*solarmass
  Rcloud = (Rcloud_0 + Rcloud_0/2.0d0 - 1.3d0*Rcloud_0*ran2(iseed))
  Rcloud = Rcloud*pc
!  f_cloud = f_cloud_0 + 2.0d0*f_cloud_0 - 2.0d0*f_cloud_0*ran2(iseed)
  f_cloud = f_cloud_0 + 0.1d0 - 0.1d0*ran2(iseed)
!
!  Mcloud = Mcloud_0*solarmass
!  Rcloud = Rcloud_0*pc
!  f_cloud = f_cloud_0 

  print*, 'Cloud props (M,R,f):', Mcloud/solarmass, Rcloud/pc, f_cloud
  print*, 'Cloud mass = ', Mcloud/solarmass, 'Cloud radius = ', Rcloud/pc

  rhocloud = Mcloud/((4.0d0*pi*Rcloud**3.0d0)/3.0d0)

  print*, 'Cloud density = ', rhocloud

  beta_cloud = (f_cloud**2.0d0)/(2.0d0*pi)
  omega_cloud = f_cloud*Sqrt(G*rhocloud)

  amax_cloud = ((3.0d0*f_cloud**2.0d0/(4.0d0*pi))*Rcloud)

  print*, 'Cloud rotation = ', f_cloud, 'max radius = ', amax_cloud/AU
  print*, 'Angular speed = ', omega_cloud

  amax_cloud = omega_cloud**2.0d0*Rcloud**4.0d0/G/Mcloud

  tff = Sqrt((3.0d0*pi)/(32.0d0*G*rhocloud))

  mdot_cloud = Mcloud/tff

  print*, 'infall rate = ',mdot_cloud, 'free fall time = ', tff/yr

  Mcloud = Mcloud/solarmass
  Rcloud = Rcloud/pc

  print*, TRIM(prefix)//'.log'

  inquire(file=TRIM(prefix)//'.log', exist=file_exists)

  print*, file_exists

  If (.not. file_exists) Then
    open(itime,file=TRIM(prefix)//'.log',access='append',status='unknown')
    write(itime,*) Mcloud, Rcloud, f_cloud, amax_cloud/AU, mdot_cloud*yr/solarmass, tff/yr, Tirr, alpha_visc, Lx
    close(itime)
  EndIf

  dMtot = 0.0d0

end subroutine setup_cloud
