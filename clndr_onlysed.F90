#if defined(SED_OFFLINE)
subroutine updcln_onlysed()

!-----------------------------------------------------------------------
! Update the calendar
!
! 2018-06   Code based on updcln() from phy/clndr.F, here adjusted for
!           monthly timestepping in the sediment() spin-up routine.
!-----------------------------------------------------------------------

use mod_xc

implicit none

#include "common_clndr.h90"

! get new date
nday=1
nmonth=nmonth+1
if (nmonth > 12) then
   nday_of_year=1
   nmonth=1
   nyear=nyear+1
   if (calendar(1:3) == 'sta') then
      if (mod(nyear,4)   == 0 .and.                                  &
       & (mod(nyear,100) /= 0 .or. mod(nyear,400) == 0)) then
         nd_in_m(2)=29
         nday_in_year=366
      else
         nd_in_m(2)=28
         nday_in_year=365
      endif
   endif
endif

if ( (calendar(1:3) == 'sta'  .or. calendar(1:3) == 'mix' .or.       &
 &    calendar(1:3) == 'gre').and. nyear <= 1582 ) then
   if (mnproc == 1) then
      write (lp,*)                                                   &
         & 'Do not use mixed Julian/Gregorian calendar before Oct 10th 1582!'
   endif
   call xcstop('(updcln)')
endif

return
end
#endif
