SUBROUTINE generate_star
! 
use stardata
use eosdata
use embryodata

! For now, simply postulate star mass (in solar masses)
! (TODO - stellar mass function?)

If (runmode .ne. 'C1') then
  mstar = 10.0d0*umass
  do while (mstar .gt. mstar1*umass) 
!     mstar = (mstar0 + mstar0*(1.0 - ran2(iseed))/1.0)*umass
    call random_number(rnum)
    print*, rnum
    mstar = (mstar0 + 2.0d0*mstar1*rnum)*umass
!    mstar = mstar0*umass
  enddo
!  mstar = mstar0*umass
Else
  mstar = 0.0d0
  mstar_collapse = 0.0d0
EndIf

END SUBROUTINE generate_star
