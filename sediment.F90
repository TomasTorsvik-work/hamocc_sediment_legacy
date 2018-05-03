#if defined(SED_OFFLINE) || defined(SED_RCLIM) || defined(SED_WCLIM)
subroutine sediment(kpie,kpje,kpke,pglat,                   &
     &              ptho,psao,ppao,prho,pddpo,pdlxp,pdlyp,  &
     &              omask)

!-----------------------------------------------------------------------
!
! *MODULE sediment* - run sediment model for multiple timesteps
!
! Copyright (C) 2017 Marco van Hulten <Marco.Hulten@uib.no>
!
! This module is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!
! Method:
!
! After spinning up the water column, this module takes over to spin up
! the sediment.  During the water column spin-up, the model collects
! monthly averages of bottom water layer particle fluxes and dissolved
! concentrations.  Those monthly averages are used by this module as a
! forcing for the sediment model.
!
! History:
!
! The sediment model is originally based on Heinze, Maier-Reimer,
! Winguth and Archer (1999) doi:10.1029/98GB02812.
!
! 2017-07   Code based on hamocc4bcm.F90.
!
!-----------------------------------------------------------------------

use mod_xc
use mo_bgcmean
use mo_control_bgc, only: io_stdo_bgc, nyear_global, maxyear_sediment
use mo_param1_bgc 
use mo_sedmnt

implicit none

#include "common_clndr.h90"   ! F90 - MvH

! Function arguments; see hamocc4bcm.F90 for variable description
!
integer, intent(in)  :: kpie,kpje,kpke
real, intent(in)     :: pglat  (kpie,kpje)
real, intent(in)     :: ptho   (kpie,kpje,kpke)
real, intent(in)     :: psao   (kpie,kpje,kpke)
real, intent(in)     :: ppao   (kpie,kpje)
real, intent(in)     :: prho   (kpie,kpje,kpke)
real, intent(in)     :: pddpo  (kpie,kpje,kpke)
real, intent(in)     :: pdlxp  (kpie,kpje)
real, intent(in)     :: pdlyp  (kpie,kpje)
real, intent(in)     :: omask  (kpie,kpje)

! Local variables
!
integer :: iyear, imonth      ! timestepping over sediment only

!-----------------------------------------------------------------------
!     Sediment module

IF (mnproc.eq.1) THEN
   WRITE(io_stdo_bgc,*) 'Starting sediment spin-up'
ENDIF

do iyear = 1, maxyear_sediment
   nyear_global = nyear_global + 1
   IF (mnproc.eq.1) WRITE(io_stdo_bgc,*) 'sediment(): nyear_global = ', nyear_global
   do imonth = 1, 12
      CALL powach_onlysed(kpie,kpje,kpke,pdlxp,pdlyp,psao,prho,omask,imonth)

#ifdef PBGC_CK_TIMESTEP
      IF (mnproc.eq.1) THEN
         WRITE(io_stdo_bgc,*)' '
         WRITE(io_stdo_bgc,*)'after powach_onlysed(): call INVENTORY'
      ENDIF
      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif

!     sediment is shifted every sediment() timestep (Euler)
      CALL SEDSHI(kpie,kpje,omask)

!     accumulate sediments
      call accsdm(jpowaic,powtra(1,1,1,ipowaic))
      call accsdm(jpowaal,powtra(1,1,1,ipowaal))
      call accsdm(jpowaph,powtra(1,1,1,ipowaph))
      call accsdm(jpowaox,powtra(1,1,1,ipowaox))
      call accsdm(jpown2 ,powtra(1,1,1,ipown2) )
      call accsdm(jpowno3,powtra(1,1,1,ipowno3))
      call accsdm(jpowasi,powtra(1,1,1,ipowasi))
      call accsdm(jssso12,sedlay(1,1,1,issso12))
      call accsdm(jssssil,sedlay(1,1,1,issssil))
      call accsdm(jsssc12,sedlay(1,1,1,isssc12))
      call accsdm(jssster,sedlay(1,1,1,issster))

!     accumulate sediment burial
      call accbur(jburssso12,burial(1,1,issso12))
      call accbur(jburssssil,burial(1,1,issssil))
      call accbur(jbursssc12,burial(1,1,isssc12))
      call accbur(jburssster,burial(1,1,issster))

   enddo
enddo

!-----------------------------------------------------------------------
return
end
#endif
