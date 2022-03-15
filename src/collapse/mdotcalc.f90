subroutine mdotcalc(t)

   use clouddata 
   use stardata
   use unitdata, only: pi
   
   real :: t

   If ((t .gt. 0.0D0).and.(t.lt.tff)) Then
      mdot = mdot_cloud
   ElseIf ((tff.le.t).and.(t.lt.2.0D0*tff)) Then
      mdot = 0.5D0*mdot_cloud*(1.0D0+COS(pi*(t-tff)/tff))
   Else
      mdot = 0.0D0
   EndIf

end subroutine mdotcalc
