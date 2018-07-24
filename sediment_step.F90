subroutine sediment_step(kpie, kpje, kpke, pglat, pddpo, pdlxp, pdlyp,  &
   &                     psao_, prho_, omask,                           &
   &                     ocetra_, bolay_, keqb_,                        &
   &                     prorca_, prcaca_, silpro_, produs_, co3_)

!-----------------------------------------------------------------------
!
! Perform one sediment timestep
!
! Copyright (C) 2018 Marco van Hulten <Marco.Hulten@uib.no> et al.
!                    Geophysical Institute @ University of Bergen
!
! This subroutine is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!
! Method:
!
! Set the sediment timestep, depending if we are in an off-line sediment
! spin-up or full MICOM/HAMOCC.  Then execute the sediment routines
! powach(), sedshi(), accsdm() and accbur().
!
! History:
!
! The sediment model is originally based on Heinze, Maier-Reimer,
! Winguth and Archer (1999) doi:10.1029/98GB02812.
!
! Marco van Hulten <Marco.Hulten@uib.no>           2017-07-12
! - Initial code based on hamocc4bcm.F90.
!
!-----------------------------------------------------------------------

use mod_xc
use mo_bgcmean
use mo_control_bgc!, only: io_stdo_bgc
use mo_param1_bgc
use mo_sedmnt

implicit none

#include "common_clndr.h90"   ! F90 - MvH
#include "common_bgcs.h90"    ! F90 - MvH

! Function arguments; see hamocc4bcm.F90 for variable description
!
integer, intent(in)  :: kpie,kpje,kpke
real, intent(in)     :: pglat  (kpie,kpje)
real, intent(in)     :: pddpo  (kpie,kpje,kpke)
real, intent(in)     :: pdlxp  (kpie,kpje)
real, intent(in)     :: pdlyp  (kpie,kpje)
real, intent(in)     :: psao_  (kpie,kpje)
real, intent(in)     :: prho_  (kpie,kpje)
real, intent(in)     :: omask  (kpie,kpje)
real, intent(inout)  :: ocetra_(kpie,kpje,nocetra)
real, intent(in)     :: bolay_ (kpie,kpje)
real, intent(in)     :: keqb_  (11,kpie,kpje)
real, intent(inout)  :: prorca_(kpie,kpje)
real, intent(inout)  :: prcaca_(kpie,kpje)
real, intent(inout)  :: silpro_(kpie,kpje)
real, intent(inout)  :: produs_(kpie,kpje)
real, intent(in)     :: co3_   (kpie,kpje)

!-----------------------------------------------------------------------
! Sediment module

if (lspinup_sediment) then
   dtsed = dtoff
   rdtsed = dtoff/dtbgc
else
   dtsed = dtbgc
   rdtsed = 1.
endif

call powach(kpie,kpje,kpke,pdlxp,pdlyp,psao_,prho_,                    &
   &        bolay_,ocetra_,keqb_,prorca_,prcaca_,silpro_,produs_,co3_)

#ifdef PBGC_CK_TIMESTEP
IF (mnproc.eq.1) THEN
   WRITE(io_stdo_bgc,*)' '
   WRITE(io_stdo_bgc,*)'after powach(): call INVENTORY'
ENDIF
CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
#endif

! hamocc4bcm() calls CARCHM(), which includes code for sediment C-14 decay.
! This should be added (TODO, e.g., here) when isotopes are included. - MvH

! sediment is shifted every sediment() timestep
call sedshi(kpie,kpje,omask)

! accumulate sediments
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

! accumulate sediment burial
call accbur(jburssso12,burial(1,1,issso12))
call accbur(jburssssil,burial(1,1,issssil))
call accbur(jbursssc12,burial(1,1,isssc12))
call accbur(jburssster,burial(1,1,issster))

!-----------------------------------------------------------------------
return
end subroutine sediment_step

