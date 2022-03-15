SUBROUTINE generate_star
! 
use stardata
use eosdata
use embryodata

! For now, simply postulate star mass (in solar masses)
! (TODO - stellar mass function?)

If (runmode .ne. 'C') then
  mstar = (mstar0 + mstar0*(ran2(iseed)-0.5)/4.0)*umass
Else
  mstar = 0.0d0
EndIf

END SUBROUTINE generate_star
