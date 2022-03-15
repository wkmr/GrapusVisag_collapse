subroutine setup_cloud
  
  use stardata
  use clouddata
  use unitdata

  real :: Mcloud, Rcloud, f_cloud

  Mcloud = (Mcloud_0 + Mcloud_0/4.0d0 - Mcloud_0*ran2(iseed)/2.0d0)
  Mcloud = Mcloud*solarmass
  Rcloud = (Rcloud_0 + Rcloud_0/4.0d0 - Rcloud_0*ran2(iseed)/2.0d0)
  Rcloud = Rcloud_0*pc
  f_cloud = f_cloud_0 + 2.0d0*f_cloud_0 - f_cloud_0*ran2(iseed)

  print*, 'Cloud props (M,R,f):', Mcloud/solarmass, Rcloud/pc, f_cloud
  print*, 'Cloud mass = ', Mcloud/solarmass, 'Cloud radius = ', Rcloud/AU

  rhocloud = Mcloud/((4.0d0*pi*Rcloud**3.0d0)/3.0d0)

  beta_cloud = (f_cloud**2.0d0)/(2.0d0*pi)
  omega_cloud = f_cloud*Sqrt(G*rhocloud)

  amax_cloud = ((3.0d0*f_cloud**2.0d0/(4.0d0*pi))*Rcloud)

  print*, 'Cloud rotation = ', f_cloud, 'max radius = ', amax_cloud/AU

  amax_cloud = omega_cloud**2.0d0*Rcloud**4.0d0/G/Mcloud

  print*, amax_cloud
  
  tff = Sqrt((3.0d0*pi)/(32.0d0*G*rhocloud))

  mdot_cloud = Mcloud/tff

  print*, 'infall rate = ',mdot_cloud, 'free fall time = ', tff/3.15d7

  Mcloud = Mcloud/solarmass
  Rcloud = Rcloud/pc

end subroutine setup_cloud
