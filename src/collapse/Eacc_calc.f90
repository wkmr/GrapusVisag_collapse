subroutine Eacc_calc(t) 

   use clouddata 
   use stardata
   use unitdata
       
   real(kind=8) :: ro, M_ro, vel, vphi, vperp
   real(kind=8) :: t

   integer i

   IF(allocated(E_acc)) deallocate(E_acc)

   allocate(E_acc(nrannuli))

   E_acc(:) = 0.0d0

   Router = ((3.0d0*dMtot)/(rhocloud*4.0d0*pi) + rin**3.0d0)**(1.0d0/3.0d0)
   MRouter = rhocloud*4.0d0*pi*Router**3.0d0/3.0d0

   If (t .le. tff) Then
     azero = amax_cloud*(t/tff)**(1.0d0/3.0d0)
   Else
     azero = amax_cloud
   EndIf

   ro = 4.0d0*azero*rhocloud*G/3.0d0/omega_cloud**2.0d0

   If (azero .gt. 0.0d0) Then
     M_ro = omega_cloud**2.0d0*ro**2.0d0/G/azero
     M_ro = M_ro*ro**2.0d0
   Else
     M_ro = 0.0d0
   EndIf

   If (azero .gt. 0.0d0) Then
     do i = isr, ier
       omega_d(i) = Sqrt(G*mstar/rz(i)**3.0d0)

       vel = G*M_ro/rz(i)*(1 + rz(i)**2.0/ro**2.0)

       vel = vel - G*M_ro/ro

       vel = Sqrt(2.0d0*vel)  

       vphi = Sqrt(G*M_ro/azero)

       if (vel > vphi) then
         vperp = Sqrt(vel**2.0d0 - vphi**2.0d0)

         E_acc(i) = 0.5d0*vperp**2.0d0
 
         E_acc(i) = E_acc(i) + 0.5d0*(vphi - rz(i)*omega_d(i))**2.0d0
       else
         E_acc(i) = 0.0d0
       endif

     enddo
   EndIf

end subroutine Eacc_calc

