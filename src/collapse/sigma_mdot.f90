subroutine sigma_mdot(t)

   use clouddata 
   use stardata
   use unitdata, only: pi
       
   real :: t

   integer i

   IF(allocated(dsigma_cloud)) deallocate(dsigma_cloud)
   IF(allocated(fM)) deallocate(fM)
  
   allocate(dsigma_cloud(nrannuli))
   allocate(fM(nrannuli))

   call mdotcalc(t)

   dsigma_cloud(:) = 0.0d0
   fM(:) = 0.0d0

   Router = ((3.0d0*dMtot)/(rhocloud*4.0d0*pi) + rin**3.0d0)**(1.0d0/3.0d0)
   MRouter = rhocloud*4.0d0*pi*Router**3.0d0/3.0d0

   If (t .le. tff) Then
     azero = amax_cloud*(t/tff)**(1.0d0/3.0d0)
   Else
     azero = amax_cloud
   EndIf

   do i = isr-1, ier
     if (rf(i) .lt. azero) Then
       fM(i) = 1.0 - Sqrt(1.0 - rf(i)/azero)
     else
       fM(i) = 1.0
     endif
   enddo
   do i = isr, ier-1
     dsigma_cloud(i) = (fM(i+1)-fM(i))*mdot/2.0d0/pi/rz(i)/(rf(i+1)-rf(i))
   enddo

   dmstar_cloud = fM(1)*mdot

end subroutine sigma_mdot

