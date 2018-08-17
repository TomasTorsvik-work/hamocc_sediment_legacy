subroutine powach(kpie,kpje,kpke,pdlxp,pdlyp,psao_,prho_,              &
   &       bolay_,ocetra_,keqb_,prorca_,prcaca_,silpro_,produs_,co3_)

!-----------------------------------------------------------------------
!
! Decomposition of particulate matter and diffusion of pore water
! constituents in the sediment
!
! Copyright (C) 2004 Ernst Maier-Reimer,  *MPI-Met, HH*
! Copyright (C) 2004 S. Legutke,          *MPI-MaD, HH*
! Copyright (C) 2018 Marco van Hulten <Marco.Hulten@uib.no>
!                    Geophysical Institute @ University of Bergen
!
! This module is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
!
! Purpose:
!
! Simulate the decomposition of particulate matter simultaneously with
! the diffusion of pore water consituents (oxygen, silicate and
! dissolved inorganic carbon).
!
! Methods:
!
! A backward approach is used (Heinze and Maier-Reimer, 1999).
!
! History:
!
! The sediment model is originally based on Heinze, Maier-Reimer,
! Winguth and Archer (1999) doi:10.1029/98GB02812.
!
! 2001-04   Changes by S. Legutke
! 2018-06   Refactored by Marco van Hulten
!
! Interface to ocean model (parameter list):
!
!  INTEGER  kpie     - 1st dimension of model grid.
!  INTEGER  kpje     - 2nd dimension of model grid.
!  INTEGER  kpke     - 3rd (vertical) dimension of model grid.
!  REAL     pdlxp    - size of scalar grid cell (1st dimension) [m].
!  REAL     pdlxp    - size of scalar grid cell (1st dimension) [m].
!  REAL     psao     - salinity [psu].
!  REAL     prho     - seawater density [g/cm^3].
!  REAL     bolay    - thickness of the bottom gridbox [m].
!  REAL     ocetra   - bottom layer ocean tracer ocetra(:,:,kbo,:).
!  REAL     keqb     - chemical equilibrium constants.
!  REAL     prorca   - sedimentation of carbon [kmol/m^2/s].
!  REAL     prcaca   - sedimentation of calcium carbonate [kmol/m^2/s].
!  REAL     silpro   - sedimentation of biogenic silica [kmol/m^2/s].
!  REAL     produs   - sedimentation of lithogenic dust [kmol/m^2/s].
!  REAL     co3      - dissolved carbonate in the bottom gridbox [mol/kg].
!
! The actual argument names sometimes end with an underscore, signifying
! that they may be passed in the call as a mutable array, or they may be
! passed as an object that must not be changed, a true Fortran 90 dummy
! variable whose mutability is normally restricted with INTENT(IN).
! We could overload the INTENT, possibly through an optional logical
! argument in an INTERFACE, but it may be more confusing than helpful.
! In the standard (.not. lspinup_sediment) case, prorca etc. will have
! been passed for the prorca_ etc. dummy arguments (prorca_ => prorca),
! so it doesn't matter to which we assign values.
! However, in the standard case we must assign to ocetra(:,:,kbo,:) as
! ocetra_ => ocetra_kbo (newly defined in hamocc4bcm.F90), not ocetra.
!
!-----------------------------------------------------------------------

use mo_carbch, only: sedfluxo, ocetra
use mo_chemcon, only: calcon
use mo_sedmnt
#if defined(SED_OFFLINE)
use mo_sedmnt_offline
#endif
use mo_biomod
use mo_control_bgc
use mo_param1_bgc
use mo_common_bgc, only: omask

implicit none

integer :: i, j, k, l
integer, intent(in)  :: kpie, kpje, kpke

real, intent(in)     :: psao_(kpie,kpje)
real, intent(in)     :: prho_(kpie,kpje)

real, intent(in)     :: pdlxp(kpie,kpje), pdlyp(kpie,kpje)
real, intent(in)     :: bolay_ (kpie,kpje)
real, intent(inout)  :: ocetra_(kpie,kpje,nocetra)
real, intent(in)     :: keqb_  (11,kpie,kpje)
real, intent(inout)  :: prorca_(kpie,kpje)
real, intent(inout)  :: prcaca_(kpie,kpje)
real, intent(inout)  :: silpro_(kpie,kpje)
real, intent(inout)  :: produs_(kpie,kpje)
real, intent(in)     :: co3_   (kpie,kpje)

real :: sedb1(kpie,0:ks), sediso(kpie,0:ks)
real :: solrat(kpie,ks), powcar(kpie,ks)
real :: aerob(kpie,ks), anaerob(kpie,ks)

real :: disso, dissot, undsa, silsat, posol
real :: umfa, denit, saln, rrho, alk, c, sit, pt
real :: K1, K2, Kb, Kw, Ks1, Kf, Ksi, K1p, K2p, K3p
real :: ah1, ac, cu, cb, cc, satlev
real :: ratc13, ratc14, rato13, rato14, poso13, poso14

integer, parameter :: niter = 5 ! number of iterations for carchm_solve

!-----------------------------------------------------------------------
! accelerated sediment
! needed for boundary layer vertilation in fast sediment routine

real :: bolven(kpie)

! A LOOP OVER J
! RJ: This loop must go from 1 to kpje in the parallel version,
!     otherways we had to do a boundary exchange

!$OMP PARALLEL DO                                                &
!$OMP&PRIVATE(sedb1,sediso,solrat,powcar,aerob,anaerob,          &
!$OMP&        disso,dissot,undsa,silsat,posol,                   &
!$OMP&        umfa,denit,saln,rrho,alk,c,sit,pt,                 &
!$OMP&        K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p,                &
!$OMP&        ah1,ac,cu,cb,cc,satlev,bolven,                     &
!$OMP&        ratc13,ratc14,rato13,rato14,poso13,poso14)
DO 8888 j=1,kpje

do k = 1, ks
   do i = 1, kpie
      solrat(i,k) = 0.
      powcar(i,k) = 0.
      anaerob(i,k)= 0.
      aerob(i,k)  = 0.
   enddo
enddo


! calculate bottom ventilation rate for scaling of sediment-water exchange
do i = 1, kpie
   bolven(i) = 1.
enddo

do k = 0, ks
   do i = 1, kpie
      sedb1(i,k) = 0.
      sediso(i,k)= 0.
   enddo
enddo


! Calculate silicate-opal cycle and simultaneous silicate diffusion
!-----------------------------------------------------------------------

! Dissolution rate constant of opal (disso) [1/(kmol Si(OH)4/m3)*1/sec]

!disso=1.e-8
disso  = 1.e-6 ! test vom 03.03.04 half live sil ca. 20.000 yr
dissot = disso*dtsed

! Silicate saturation concentration is 1 mol/m3

silsat = 0.001

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

do i = 1, kpie
   if(omask(i,j) > 0.5) then
      undsa = silsat - powtra(i,j,1,ipowasi)
      sedb1(i,0) = bolay_(i,j) * (silsat - ocetra_(i,j,isilica))  &
         &       * bolven(i)
      solrat(i,1) = ( sedlay(i,j,1,issssil)                       &
         &        + silpro_(i,j)/(porsol(1)*seddw(1)) )           &
         &        * dissot/(1.+dissot*undsa) * porsol(1)/porwat(1)
   endif
enddo


! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         undsa = silsat - powtra(i,j,k,ipowasi)
         sedb1(i,k) = seddw(k)*porwat(k)*(silsat-powtra(i,j,k,ipowasi))
         if ( k > 1 ) solrat(i,k) = sedlay(i,j,k,issssil)         &
            &         * dissot/(1.+dissot*undsa) * porsol(k)/porwat(k)
      endif
   enddo
enddo

! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.

call powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,bolay_)

! Store the flux for budget.
! Add sedimentation to first layer.

do i = 1, kpie
   if(omask(i,j) > 0.5) then
      if ( .not. lspinup_sediment ) then
         sedfluxo(i,j,ipowasi) = sedfluxo(i,j,ipowasi) +                   &
            &    (silsat-sediso(i,0) - ocetra_(i,j,isilica))*bolay_(i,j)

         ! NOTE: If ocetra_(:,:,:) instead of ocetra(:,:,kbo,:) were to be updated,
         !       the latter would need to be updated by the former later on (MvH)!
         !
         ocetra(i,j,kbo(i,j),isilica) = silsat - sediso(i,0)
      endif
      sedlay(i,j,1,issssil) =                                           &
         &  sedlay(i,j,1,issssil) + silpro_(i,j) / (porsol(1)*seddw(1)) &
         &                                       * rdtsed
   endif
enddo

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! Update pore water concentration from new undersaturation.

do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         umfa = porsol(k)/porwat(k)
         solrat(i,k) = sedlay(i,j,k,issssil)                        &
            &        * dissot / (1. + dissot*sediso(i,k))
         posol = sediso(i,k)*solrat(i,k)
         sedlay(i,j,k,issssil) = sedlay(i,j,k,issssil) - posol
         powtra(i,j,k,ipowasi) = silsat - sediso(i,k)
      endif
   enddo
enddo

! Calculate oxygen-POC cycle and simultaneous oxygen diffusion
!-----------------------------------------------------------------------

! Degradation rate constant of POP (disso) [1/(kmol O2/m3)*1/sec]

disso  = 0.01/86400.  !  disso=3.e-5 was quite high
dissot = disso*dtsed

! This scheme is not based on undersaturation, but on O2 itself

! Evaluate boundary conditions for sediment-water column exchange.
! Current concentration of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

do i = 1, kpie
   if(omask(i,j) > 0.5) then
      undsa = powtra(i,j,1,ipowaox)
      sedb1(i,0) = bolay_(i,j)*ocetra_(i,j,ioxygen)         &
         &       * bolven(i)
      solrat(i,1) = (sedlay(i,j,1,issso12)+prorca_(i,j)/(porsol(1)*seddw(1)))  &
         &        * ro2ut * dissot/(1.+dissot*undsa) * porsol(1)/porwat(1)
   endif
enddo

! Evaluate sediment concentration and degradation.
! Current concentration in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

do k = 1, ks
   do i = 1, kpie
      if(bolay_(i,j) > 0.) then
         undsa = powtra(i,j,k,ipowaox)
         sedb1(i,k) = seddw(k)*porwat(k)*powtra(i,j,k,ipowaox)
         if (k > 1) solrat(i,k) = sedlay(i,j,k,issso12)               &
            &       * ro2ut*dissot/(1.+dissot*undsa) * porsol(k)/porwat(k)
      endif
   enddo
enddo

! Solve for new O2 concentration sediso, from current concentration sedb1,
! and first guess of new solid sediment solrat.

call powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,bolay_)

! Update water column oxygen, and store the flux for budget (opwflux).
! Add sedimentation to first layer.

do i = 1, kpie
   if(omask(i,j) > 0.5) then
      if ( .not. lspinup_sediment ) then
         ocetra(i,j,kbo(i,j),ioxygen) = sediso(i,0)
      endif
      sedlay(i,j,1,issso12) = sedlay(i,j,1,issso12)                  &
         &                  + prorca_(i,j) / (porsol(1)*seddw(1))    &
         &                  * rdtsed
#ifdef __c_isotopes
      sedlay(i,j,1,issso13)                                     &
         &      = sedlay(i,j,1,issso13)+pror13(i,j)/(porsol(1)*seddw(1))
      sedlay(i,j,1,issso14)                                     &
         &      = sedlay(i,j,1,issso14)+pror14(i,j)/(porsol(1)*seddw(1))
#endif
      if ( .not. lspinup_sediment ) then
         prorca(i,j) = 0.
#ifdef __c_isotopes
         pror13(i,j) = 0.
         pror14(i,j) = 0.
#endif
      endif
   endif
enddo


! Calculate updated degradation rate from updated concentration.
! Calculate new solid sediment.
! Update pore water concentration.
! Store flux in array aerob, for later computation of DIC and alkalinity.
do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         umfa = porsol(k)/porwat(k)
         solrat(i,k) = sedlay(i,j,k,issso12) * dissot/(1.+dissot*sediso(i,k))
         posol = sediso(i,k)*solrat(i,k)
#ifdef __c_isotopes
         rato13 = sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24)
         rato14 = sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24)
         poso13 = posol*rato13
         poso14 = posol*rato14
#endif
         aerob(i,k) = posol*umfa !this has P units: kmol P/m3 of pore water
         sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12)-posol
         powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph)+posol*umfa
         powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3)+posol*rnit*umfa
         powtra(i,j,k,ipowaox) = sediso(i,k)
#ifdef __c_isotopes
         sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13)-poso13
         sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14)-poso14
         powtra(i,j,k,ipowc13) = powtra(i,j,k,ipowc13)+poso13*umfa
         powtra(i,j,k,ipowc14) = powtra(i,j,k,ipowc14)+poso14*umfa
#endif
      endif
   enddo
enddo

! Calculate nitrate reduction under anaerobic conditions explicitely
!-----------------------------------------------------------------------

! Denitrification rate constant of POP (disso) [1/sec]
! Store flux in array anaerob, for later computation of DIC and alkalinity.

denit = 0.01/86400. *dtsed
do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         if(powtra(i,j,k,ipowaox) < 1.e-6) then
            posol = denit*MIN(0.5*powtra(i,j,k,ipowno3)/114.,          &
               &                  sedlay(i,j,k,issso12))
            umfa = porsol(k)/porwat(k)
            anaerob(i,k) = posol*umfa !this has P units: kmol P/m3 of pore water
            sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12)-posol
            powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph)+posol*umfa
            powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3)-98.*posol*umfa
            powtra(i,j,k,ipown2) = powtra(i,j,k,ipown2)+57.*posol*umfa
#ifdef __c_isotopes
            rato13 = sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24)
            rato14 = sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24)
            poso13 = posol*rato13
            poso14 = posol*rato14
            sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13)-poso13
            sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14)-poso14
            powtra(i,j,k,ipowc13) = powtra(i,j,k,ipowc13)+poso13*umfa
            powtra(i,j,k,ipowc14) = powtra(i,j,k,ipowc14)+poso14*umfa
#endif
         endif
      endif
   enddo
enddo

!    sulphate reduction in sediments
do k = 1, ks
   do i = 1, kpie
      if (omask(i,j) > 0.5) then
      if (powtra(i,j,k,ipowaox) < 3.e-6 .and. powtra(i,j,k,ipowno3) < 3.e-6) then
         posol = denit* sedlay(i,j,k,issso12)        ! remineralization of POC
         umfa = porsol(k)/porwat(k)
         anaerob(i,k) = anaerob(i,k)+posol*umfa !this has P units: kmol P/m3 of pore water
                                              !this overwrites anaerob from denitrification

         sedlay(i,j,k,issso12) = sedlay(i,j,k,issso12)-posol
         powtra(i,j,k,ipowaph) = powtra(i,j,k,ipowaph)+posol*umfa
         powtra(i,j,k,ipowno3) = powtra(i,j,k,ipowno3)+posol*umfa*rno3
#ifdef __c_isotopes
         rato13 = sedlay(i,j,k,issso13)/(sedlay(i,j,k,issso12)+1.e-24)
         rato14 = sedlay(i,j,k,issso14)/(sedlay(i,j,k,issso12)+1.e-24)
         poso13 = posol*rato13
         poso14 = posol*rato14
         sedlay(i,j,k,issso13) = sedlay(i,j,k,issso13)-poso13

         sedlay(i,j,k,issso14) = sedlay(i,j,k,issso14)-poso14
         powtra(i,j,k,ipowc13) = powtra(i,j,k,ipowc13)+poso13*umfa
         powtra(i,j,k,ipowc14) = powtra(i,j,k,ipowc14)+poso14*umfa
#endif
      endif
      endif
   enddo
enddo


! Calculate CaCO3-CO3 cycle and simultaneous CO3-undersaturation diffusion
!-----------------------------------------------------------------------


! Compute new powcar, carbonate ion concentration in the sediment
! from changed alkalinity (nitrate production during remineralisation)
! and DIC gain. Iterate 5 times. This changes pH (sedhpl) of sediment.

do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         saln= psao_(i,j)
         rrho= prho_(i,j)
         alk = (powtra(i,j,k,ipowaal)-(anaerob(i,k)+aerob(i,k))*16.)  / rrho
         c   = (powtra(i,j,k,ipowaic)+(anaerob(i,k)+aerob(i,k))*122.) / rrho
         sit =  powtra(i,j,k,ipowasi) / rrho
         pt  =  powtra(i,j,k,ipowaph) / rrho
         ah1 = sedhpl(i,j,k)
         K1  = keqb_( 1,i,j)
         K2  = keqb_( 2,i,j)
         Kb  = keqb_( 3,i,j)
         Kw  = keqb_( 4,i,j)
         Ks1 = keqb_( 5,i,j)
         Kf  = keqb_( 6,i,j)
         Ksi = keqb_( 7,i,j)
         K1p = keqb_( 8,i,j)
         K2p = keqb_( 9,i,j)
         K3p = keqb_(10,i,j)

         call carchm_solve(saln,c,alk,sit,pt,                  &
                           K1,K2,Kb,Kw,Ks1,Kf,Ksi,K1p,K2p,K3p, &
                           ah1,ac,niter)

         cu = ( 2. * c - ac ) / ( 2. + K1 / ah1 )
         cb = K1 * cu / ah1
         cc = K2 * cb / ah1
         sedhpl(i,j,k) = max(1.e-20,ah1)
         powcar(i,k)   = cc * rrho
      endif
   enddo
enddo


! Dissolution rate constant of CaCO3 (disso) [1/(kmol CO3--/m3)*1/sec]
disso = 1.e-7
dissot = disso*dtsed

! Evaluate boundary conditions for sediment-water column exchange.
! Current undersaturation of bottom water: sedb(i,0) and
! Approximation for new solid sediment, as from sedimentation flux: solrat(i,1)

! CO3 saturation concentration is aksp/calcon as in CARCHM
! (calcon defined in BELEG_BGC with 1.03e-2; 1/calcon =~ 97.)

do i = 1, kpie
   if(omask(i,j) > 0.5) then
      satlev = keqb_(11,i,j)/calcon+2.e-5
      undsa = MAX(satlev-powcar(i,1),0.)
      sedb1(i,0) = bolay_(i,j) * (satlev-co3_(i,j)) * bolven(i)
      solrat(i,1) = (sedlay(i,j,1,isssc12)+prcaca_(i,j) / (porsol(1)*seddw(1)))  &
         &        * dissot / (1.+dissot*undsa) * porsol(1)/porwat(1)
   endif
enddo

! Evaluate sediment undersaturation and degradation.
! Current undersaturation in pore water: sedb(i,k) and
! Approximation for new solid sediment, as from degradation: solrat(i,k)

do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         undsa=MAX(keqb_(11,i,j)/calcon-powcar(i,k),0.)
         sedb1(i,k) = seddw(k)*porwat(k)*undsa
         if (k > 1) solrat(i,k) = sedlay(i,j,k,isssc12)                 &
            &                   * dissot/(1.+dissot*undsa) * porsol(k)/porwat(k)
         if (undsa <= 0.) solrat(i,k) = 0.
      endif
   enddo
enddo

! Solve for new undersaturation sediso, from current undersaturation sedb1,
! and first guess of new solid sediment solrat.

call powadi(j,kpie,kpje,solrat,sedb1,sediso,bolven,bolay_)

! There is no exchange between water and sediment with respect to co3 so far.
! Add sedimentation to first layer.
do i = 1, kpie
   if(omask(i,j) > 0.5) then
      sedlay(i,j,1,isssc12) = sedlay(i,j,1,isssc12)                  &
         &                  + prcaca_(i,j) / (porsol(1)*seddw(1))    &
         &                  * rdtsed
#ifdef __c_isotopes
      sedlay(i,j,1,isssc13) = sedlay(i,j,1,isssc13)+prca13(i,j)/(porsol(1)*seddw(1))
      sedlay(i,j,1,isssc14) = sedlay(i,j,1,isssc14)+prca14(i,j)/(porsol(1)*seddw(1))
#endif
      if ( .not. lspinup_sediment ) then
         prcaca(i,j) = 0.
#ifdef __c_isotopes
         prca13(i,j) = 0.
         prca14(i,j) = 0.
#endif
      endif
   endif
enddo

! Calculate updated degradation rate from updated undersaturation.
! Calculate new solid sediment.
! No update of powcar pore water concentration from new undersaturation so far.
! Instead, only update DIC, and, of course, alkalinity.
! This also includes gains from aerobic and anaerobic decomposition.

do k = 1, ks
   do i = 1, kpie
      if(omask(i,j) > 0.5) then
         umfa = porsol(k)/porwat(k)
         solrat(i,k) = sedlay(i,j,k,isssc12)                           &
            &        * dissot / (1.+dissot*sediso(i,k))
         posol = sediso(i,k)*solrat(i,k)
         sedlay(i,j,k,isssc12) = sedlay(i,j,k,isssc12)-posol
         powtra(i,j,k,ipowaic) = powtra(i,j,k,ipowaic)                 &
            &                  + posol*umfa+(aerob(i,k)+anaerob(i,k))*122.
         powtra(i,j,k,ipowaal) = powtra(i,j,k,ipowaal)                 &
            &                  + 2.*posol*umfa-16.*(aerob(i,k)+anaerob(i,k))
#ifdef __c_isotopes
         ratc13 = sedlay(i,j,k,isssc13)/(sedlay(i,j,k,isssc12)+1.e-24)
         ratc14 = sedlay(i,j,k,isssc14)/(sedlay(i,j,k,isssc12)+1.e-24)
         poso13 = posol*ratc13
         poso14 = posol*ratc14
         sedlay(i,j,k,isssc13) = sedlay(i,j,k,isssc13)-poso13
         sedlay(i,j,k,isssc14) = sedlay(i,j,k,isssc14)-poso14
         powtra(i,j,k,ipowc13) = powtra(i,j,k,ipowc13)+poso13*umfa
         powtra(i,j,k,ipowc14) = powtra(i,j,k,ipowc14)+poso14*umfa
#endif
      endif
   enddo
enddo

8888  CONTINUE
!$OMP END PARALLEL DO

call dipowa(kpie,kpje,kpke,pdlxp,pdlyp,omask,                        &
   &        bolay_)

!ik add clay sedimentation onto sediment
!ik this is currently assumed to depend on total and corg sedimentation:
!ik f(POC) [kg C] / f(total) [kg] = 0.05
!ik thus it is
!$OMP PARALLEL DO
do j = 1, kpje
   do i = 1, kpie
      sedlay(i,j,1,issster) = sedlay(i,j,1,issster)                 &
         &                  + produs_(i,j) / (porsol(1)*seddw(1))   &
         &                  * rdtsed
   enddo
enddo
!$OMP END PARALLEL DO

if ( .not. lspinup_sediment ) then
!$OMP PARALLEL DO
do j = 1, kpje
   do i = 1, kpie
      silpro_(i,j) = 0.
      prorca_(i,j) = 0.
#ifdef __c_isotopes
      pror13(i,j) = 0.
      pror14(i,j) = 0.
      prca13(i,j) = 0.
      prca14(i,j) = 0.
#endif
      prcaca_(i,j) = 0.
      produs_(i,j) = 0.
   enddo
enddo
!$OMP END PARALLEL DO
endif

return
end
