#if defined(SED_OFFLINE)
module mo_sedmnt_offline

!-----------------------------------------------------------------------
!
! Routines for offline sediment spin-up
!
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
! Routines for reading, calculating and writing a bottom sediment water
! climatology, and/or spinning up the sediment.
!
! Description:
!
! Depending on the truth values of lsed_*, this subroutine will
!     - run the sediment routines, just as hamocc4bcm() does, but by using
!       a different timestep and different bottom-seawater properties;
!     - read a bottom-seawater climatology from file (otherwise full HAMOCC
!       will run for one year to get a climatology);
!     - write a bottom-seawater climatology based on a preceding one-year
!       full HAMOCC simulation.
!
! - subroutine read_clim()
!     Read a bottom water climatology from file.
!
! - subroutine write_clim()
!     Write a bottom water climatology to file.  This is done based on
!     one-year MICOM/HAMOCC run over which monthly averages are
!     accumulated.
!
! - subroutine sedmnt_offline()
!     Main off-line sediment routine.  Typically, after spinning up the
!     water column, this subprogram should be used to spin up the
!     sediment.  If lsed_rclim is set (namelist), lread_clim will be set
!     and it will read_clim().  This is only done one time.
!     In case no climatology was read, it will first collect monthly
!     averages of bottom water layer particle fluxes and dissolved
!     concentrations by running HAMOCC for at least one year.  Then, if
!     lsed_wclim is set, it will write the monthly averages to file
!     through write_clim().  Finally, it will spin up the sediment by
!     looping sediment_step() over maxyear model years.
!
! - subroutine updcln_onlysed()
!     Update the calendar based on updcln() from phy/clndr.F, here
!     adjusted for monthly timestepping in the sediment() spin-up
!     routine.
!     NOTE: At the end of sedmnt_offline(), the calendar is set back to
!     MICOM's last known date.
!
! - subroutine alloc_mem_sedmnt_offline()
!     Allocate memory for bottom-water climatology, analogous to
!     mo_sedmnt's alloc_mem_sedmnt().
!
! Refer to the flowchart (version 1.1) in the documentation.
!
! History:
!
! 2018-06   Initial version by MvH.
! 2018-06   updcln_onlysed() based on updcln() from phy/clndr.F.
! 2018-07   Rewrite.
!
!-----------------------------------------------------------------------

use mo_control_bgc
use mo_carbch,      only: ocetra, keqb, co3
use mo_param1_bgc,  only: nocetra, idet, icalc, iopal, ifdust        &
   &                    , nsedtra, npowtra
use mo_biomod,      only: kbo, bolay, wpoc, wmin, wmax, wdust        &
   &                    , wlin, wopal, wcal
use mod_xc
use mo_common_bgc
use mo_bgcmean


implicit none

#include "common_clndr.h90"
#include "common_bgcs.h90"

public

real, dimension (:,:,:),   allocatable :: ocetra_kbo_avg
real, dimension (:,:,:,:), allocatable :: ocetra_kbo_clim
real, dimension (:,:),     allocatable :: prorca_avg
real, dimension (:,:,:),   allocatable :: prorca_clim
real, dimension (:,:),     allocatable :: prcaca_avg
real, dimension (:,:,:),   allocatable :: prcaca_clim
real, dimension (:,:),     allocatable :: silpro_avg
real, dimension (:,:,:),   allocatable :: silpro_clim
real, dimension (:,:),     allocatable :: produs_avg
real, dimension (:,:,:),   allocatable :: produs_clim
real, dimension (:,:,:),   allocatable :: keqb_avg
real, dimension (:,:,:,:), allocatable :: keqb_clim
real, dimension (:,:),     allocatable :: bolay_avg
real, dimension (:,:,:),   allocatable :: bolay_clim
real, dimension (:,:),     allocatable :: bgc_t_kbo_avg
real, dimension (:,:,:),   allocatable :: bgc_t_kbo_clim
real, dimension (:,:),     allocatable :: bgc_s_kbo_avg
real, dimension (:,:,:),   allocatable :: bgc_s_kbo_clim
real, dimension (:,:),     allocatable :: bgc_rho_kbo_avg
real, dimension (:,:,:),   allocatable :: bgc_rho_kbo_clim
real, dimension (:,:),     allocatable :: co3_kbo_avg
real, dimension (:,:,:),   allocatable :: co3_kbo_clim

! averaging and writing frequencies for diagnostic output
integer, save      :: nsed
integer, parameter :: nsedmax = 5
real,    dimension(nsedmax), save :: diagfq_sed,filefq_sed
integer, dimension(nsedmax), save :: nacc_sed
logical, dimension(nsedmax), save :: diagmon_sed, diagann_sed,    &
   &                                 diagdec_sed, diagcen_sed,    &
   &                                 diagmil_sed,                 &
   &                                 filemon_sed, fileann_sed,    &
   &                                 filedec_sed, filecen_sed,    &
   &                                 filemil_sed, sedwrt

! namelist for diagnostic output
!
integer, dimension(nsedmax), save ::                                &
   & SDM_POWAIC    =0    ,SDM_POWAAL    =0    ,SDM_POWAPH    =0  ,  &
   & SDM_POWAOX    =0    ,SDM_POWN2     =0    ,SDM_POWNO3    =0  ,  &
   & SDM_POWASI    =0    ,SDM_SSSO12    =0    ,SDM_SSSSIL    =0  ,  &
   & SDM_SSSC12    =0    ,SDM_SSSTER    =0                       ,  &
   & BUR_SSSO12    =0    ,BUR_SSSC12    =0    ,BUR_SSSSIL    =0  ,  &
   & BUR_SSSTER    =0                                            ,  &
   & GLB_AVEPERIO  =0    ,GLB_FILEFREQ  =0    ,GLB_COMPFLAG  =0  ,  &
   & GLB_NCFORMAT  =0    ,GLB_INVENTORY =0
character(len=10), dimension(nsedmax), save :: GLB_FNAMETAG
namelist /DIASED/                                                   &
   & SDM_POWAIC        ,SDM_POWAAL        ,SDM_POWAPH        ,      &
   & SDM_POWAOX        ,SDM_POWN2         ,SDM_POWNO3        ,      &
   & SDM_POWASI        ,SDM_SSSO12        ,SDM_SSSSIL        ,      &
   & SDM_SSSC12        ,SDM_SSSTER                           ,      &
   & BUR_SSSO12        ,BUR_SSSC12        ,BUR_SSSSIL        ,      &
   & BUR_SSSTER                                              ,      &
   & GLB_AVEPERIO      ,GLB_FILEFREQ      ,GLB_COMPFLAG      ,      &
   & GLB_NCFORMAT      ,GLB_FNAMETAG      ,GLB_INVENTORY

! private variables and subprograms
!
integer, private :: i,j,iocetra,n
character(len = *), parameter, private :: bscfnmbase = "bottom_seawater_clim"
character(len =30),            private :: bscfnm
character(len = 4),            private :: seqstring
integer, private :: nday_of_month ! current day of month

private :: read_clim
private :: updcln_onlysed

contains

subroutine read_clim(nstep)
   integer, intent(in) :: nstep
   integer             :: imonth

   ! Read the bottom seawater climatology from the netCDF forcing file
   if (mnproc == 1) write(io_stdo_bgc,*)                             &
      &        'hamocc_step(): starting bottom seawater climatology read loop'
   do imonth = 1,12
      bscfnm = bscfnmbase//".read.nc"
      call aufr_bgc_onlysed(idm,jdm,kdm,imonth,bscfnm)
      ocetra_kbo_clim(:,:,:,imonth) = ocetra_kbo_avg
      bolay_clim(:,:,imonth) = bolay_avg
      keqb_clim(:,:,:,imonth) = keqb_avg
      prorca_clim(:,:,imonth) = prorca_avg
      prcaca_clim(:,:,imonth) = prcaca_avg
      silpro_clim(:,:,imonth) = silpro_avg
      produs_clim(:,:,imonth) = produs_avg
      bgc_t_kbo_clim(:,:,imonth) = bgc_t_kbo_avg
      bgc_s_kbo_clim(:,:,imonth) = bgc_s_kbo_avg
      bgc_rho_kbo_clim(:,:,imonth) = bgc_rho_kbo_avg
      co3_kbo_clim(:,:,imonth) = co3_kbo_avg
   enddo
   lcompleted_clim = .true.

   ocetra_kbo_avg = 0.0
   bolay_avg = 4000.0 ! maximum of lowest bottom layer thickness
   keqb_avg  = 0.0
   prorca_avg = 0.0
   prcaca_avg = 0.0
   silpro_avg = 0.0
   produs_avg = 0.0
   bgc_t_kbo_avg = 0.0
   bgc_s_kbo_avg = 0.0
   bgc_rho_kbo_avg = 0.0
   co3_kbo_avg = 0.0
end subroutine read_clim

subroutine prepare_clim(kpie, kpje, kpke, maxyear, nstep)
   ! Subprogram arguments
   !
   integer, intent(in)                          :: kpie,kpje,kpke
   integer, intent(in)                          :: maxyear
   integer, intent(in)                          :: nstep

   ! Accumulate the bottom seawater fields from HAMOCC
   nstep_in_month  = nstep_in_month + 1
   do iocetra = 1, nocetra
      do j = 1, jdm
         do i = 1, idm
            ocetra_kbo_avg(i,j,iocetra) = ocetra_kbo_avg(i,j, iocetra)  &
               &                        + ocetra(i,j,kbo(i,j),iocetra)
         enddo
      enddo
   enddo
   do j = 1, jdm
      do i = 1, idm
         bgc_t_kbo_avg(i,j) = bgc_t_kbo_avg(i,j) + bgc_t(i,j,kbo(i,j))
         bgc_s_kbo_avg(i,j) = bgc_s_kbo_avg(i,j) + bgc_s(i,j,kbo(i,j))
         bgc_rho_kbo_avg(i,j) = bgc_rho_kbo_avg(i,j) + bgc_rho(i,j,kbo(i,j))
         co3_kbo_avg(i,j) = co3_kbo_avg(i,j) + co3(i,j,kbo(i,j))
         bolay_avg(i,j) = min(bolay_avg(i,j), bolay(i,j))
      enddo
   enddo
   keqb_avg  = keqb_avg  + keqb
   !bolay_avg = bolay_avg + bolay

   ! Calculate the day of month
   nday_of_month = nday_of_year
   do i = 1, nmonth - 1
      nday_of_month = nday_of_month - nd_in_m(i)
   enddo

   if ( nday_of_month==nd_in_m(nmonth) .and. is_end_of_day ) then

      ! Calculate tracer monthly average
      if (mnproc == 1) write(io_stdo_bgc,*)                             &
         &  'hamocc_step(): end of month, set tracer avg for last month'
      ocetra_kbo_avg = ocetra_kbo_avg / nstep_in_month
      !bolay_avg = bolay_avg / nstep_in_month
      keqb_avg = keqb_avg / nstep_in_month
      bgc_t_kbo_avg = bgc_t_kbo_avg / nstep_in_month
      bgc_s_kbo_avg = bgc_s_kbo_avg / nstep_in_month
      bgc_rho_kbo_avg = bgc_rho_kbo_avg / nstep_in_month
      co3_kbo_avg = co3_kbo_avg / nstep_in_month
      do j = 1, jdm     ! Here we assume the WLIN case.
         do i = 1, idm
            wpoc = min(wmin+wlin*bolay_avg(i,j),wmax)
            prorca_avg(i,j)=ocetra_kbo_avg(i,j,idet  )*wpoc
            prcaca_avg(i,j)=ocetra_kbo_avg(i,j,icalc )*wcal
            silpro_avg(i,j)=ocetra_kbo_avg(i,j,iopal )*wopal
            produs_avg(i,j)=ocetra_kbo_avg(i,j,ifdust)*wdust
         enddo
      enddo
      nstep_in_month = 0

      ! Write tracer monthly averages to netCDF file
      if (lsed_wclim) then
         write (seqstring,'(I0.4)') nyear
         bscfnm = bscfnmbase//"."//seqstring//".nc"
         call aufw_bgc_onlysed(idm,jdm,kdm,nyear,nmonth,nday,nstep,bscfnm)
      endif

      if (lsed_spinup) then
         ! Save averages of last month in climatology matrix
         ocetra_kbo_clim(:,:,:,nmonth) = ocetra_kbo_avg
         bolay_clim(:,:,nmonth) = bolay_avg
         keqb_clim(:,:,:,nmonth) = keqb_avg
         prorca_clim(:,:,nmonth) = prorca_avg
         prcaca_clim(:,:,nmonth) = prcaca_avg
         silpro_clim(:,:,nmonth) = silpro_avg
         produs_clim(:,:,nmonth) = produs_avg
         bgc_t_kbo_clim(:,:,nmonth) = bgc_t_kbo_avg
         bgc_s_kbo_clim(:,:,nmonth) = bgc_s_kbo_avg
         bgc_rho_kbo_clim(:,:,nmonth) = bgc_rho_kbo_avg
         co3_kbo_clim(:,:,nmonth) = co3_kbo_avg

         if ( nmonth == 12 ) then ! we have a complete climatology
            lcompleted_clim = .true.
         endif
      endif

      ocetra_kbo_avg = 0.0
      bolay_avg = 4000.0
      keqb_avg  = 0.0
      prorca_avg = 0.0
      prcaca_avg = 0.0
      silpro_avg = 0.0
      produs_avg = 0.0
      bgc_t_kbo_avg = 0.0
      bgc_s_kbo_avg = 0.0
      bgc_rho_kbo_avg = 0.0
      co3_kbo_avg = 0.0
   endif ! end of month routines
end subroutine prepare_clim

subroutine updcln_onlysed()

!-----------------------------------------------------------------------
! Update the calendar
!-----------------------------------------------------------------------

! get new date
nday=1
if (nmonth == 1) then
   nday_of_year=1
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
   &  calendar(1:3) == 'gre').and. nyear <= 1582 ) then
   if (mnproc == 1) then
      write (lp,*)                                                   &
         & 'Do not use mixed Julian/Gregorian calendar before Oct 10th 1582!'
   endif
   call xcstop('(updcln)')
endif

return
end subroutine updcln_onlysed

subroutine sedmnt_offline(kpie, kpje, kpke, maxyear, nstep,            &
   &                      pglat, pddpo, pdlxp, pdlyp, omask)

   ! Arrays needed to update two time levels after the sediment spin-up
   !
   use mo_sedmnt, only: sedlay, powtra, burial

   ! Subprogram arguments
   !
   integer, intent(in)                          :: kpie,kpje,kpke
   integer, intent(in)                          :: maxyear
   integer, intent(in)                          :: nstep
   real, dimension(kpie,kpje), intent(in)       :: pglat
   real, dimension(kpie,kpje,kpke), intent(in)  :: pddpo
   real, dimension(kpie,kpje), intent(in)       :: pdlxp, pdlyp
   real, dimension(kpie,kpje), intent(in)       :: omask

   integer  :: nday_save, nmonth_save   ! calendar variables outside sedmnt_offline()
   integer  :: nyear_save, nday_of_year_save, nd_in_m_2, nday_in_year_save
   integer  :: nacc_bgc_1, nacc_bgc_2

   ! Read the bottom water climatology
   !
   if (lread_clim) then
      call read_clim(nstep)      ! Read the bottom water climatology,
      lread_clim = .false.       !  but only one time.
   endif

   if ( lsed_spinup .and. lcompleted_clim ) then
      ! increase the off-line sediment integration counter
      nburst = nburst + 1

      lspinning_up_sed = .true.
      if (mnproc == 1) write(io_stdo_bgc,'(a,i4)')                         &
         &     'sedmnt_offline(): stand-alone sediment spin-up starting, nburst =', nburst

      ! save calendar variables
      nday_save         = nday
      nmonth_save       = nmonth
      nyear_save        = nyear
      nday_of_year_save = nday_of_year
      nd_in_m_2         = nd_in_m(2)
      nday_in_year_save = nday_in_year
      nacc_bgc_1        = nacc_bgc(1)
      nacc_bgc_2        = nacc_bgc(2)

      do nyear = 1, maxyear
         if (mnproc == 1) write(io_stdo_bgc,'(a,i6)')                      &
               &         'sedmnt_offline(): nyear_global = ', nyear_global
         do nmonth = 1, 12
            ! do a monthly calendar update only for sediment_step()
            call updcln_onlysed()

            ! set up sediment layers (much higher diffusion rate; correct timestep)
            dtoff = 3600*24*nd_in_m(nmonth)
            call bodensed(kpie,kpje,kpke,pddpo)

            call sediment_step(idm,jdm,kdm,pglat, bgc_dp,bgc_dx,bgc_dy,    &
               & bgc_s_kbo_clim(:,:,nmonth), bgc_rho_kbo_clim(:,:,nmonth), &
               & omask,                                                    &
               & ocetra_kbo_clim(:,:,nmonth,:), bolay_clim(:,:,nmonth),    &
               & keqb_clim(:,:,:,nmonth),                                  &
               & prorca_clim(:,:,nmonth), prcaca_clim(:,:,nmonth),         &
               & silpro_clim(:,:,nmonth), produs_clim(:,:,nmonth),         &
               & co3_kbo_clim(:,:,nmonth))

!TODO: 1. update sedwrt() instead of setting..? (inner/outer loop)
!      2. inner/outer loop
!      3. no nstep, nday_... but use nyear and nmonth counters.

            ! write monthly output fields
            do n = 1, nsed
               nacc_sed(n) = 1 ! mo (#timesteps)
               if diagmon_sed(n) then
                  if (GLB_INVENTORY(n) /= 0) then                          &
                     call INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
                  endif
                  call ncwrt_onlysed(n) ! ncout_onlysed.F90
                  nacc_sed(n)=0
               endif
            enddo

         enddo

         !FIXME: Following not very useful because no special code in ncout_onlysed.F90.
         !       Moreover, now all other cases are excluded.  Better: something with fq.
         !TODO:  Can also combine below in one loop / if .. or .. or ... -> DONE.

         ! write yearly and multi-yearly output fields
         do n = 1, nsed
            nacc_sed(n) = nacc_sed(n) + 12 ! mo (#timesteps)
            if ( (diagann_sed(n)                              .or.      &
               & (diagdec_sed(n) .and. mod(nyear,   10) == 0) .or.      &
               & (diagcen_sed(n) .and. mod(nyear,  100) == 0) .or.      &
               & (diagmil_sed(n) .and. mod(nyear, 1000) == 0) ) then
               if (GLB_INVENTORY(n) /= 0) then                          &
                  call INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
               endif
               call ncwrt_onlysed(n) ! ncout_onlysed.F90
               nacc_sed(n) = 0
            endif
         enddo

         nyear_global = nyear_global + 1
      enddo

      lspinning_up_sed = .false.
      if (mnproc == 1) write(io_stdo_bgc,*)                             &
         &     'sedmnt_offline(): stand-alone sediment spin-up ended'

      ! set up sediment layers for normal MICOM/HAMOCC use
      call bodensed(kpie,kpje,kpke,pddpo)

      ! set calendar variables back to original values
      nday         = nday_save
      nmonth       = nmonth_save
      nyear        = nyear_save
      nday_of_year = nday_of_year_save
      nd_in_m(2)   = nd_in_m_2
      nday_in_year = nday_in_year_save
      nacc_bgc(1)  = nacc_bgc_1
      nacc_bgc(2)  = nacc_bgc_2

      ! update time levels, such that no old time levels will be used later
      sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
      sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:)
      powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
      powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:)
      burial2(:,:,1,:)         = burial(:,:,:)
      burial2(:,:,2,:)         = burial(:,:,:)
   endif
   lcompleted_clim = .false.  ! don't reuse the climatology
end subroutine sedmnt_offline

subroutine alloc_mem_sedmnt_offline(kpie, kpje)
   integer :: kpie, kpje
   integer :: errstat

   ! determine number of output groups
   nsed = 0
   do n = 1, nsedmax
      if (GLB_AVEPERIO(n) /= 0) then
         nsed = nsed + 1
         nacc_sed(n) = 0
      endif
   enddo

   do n = 1, nsed
      GLB_FILEFREQ(n) = max(GLB_AVEPERIO(n), GLB_FILEFREQ(n))
      !diagfq_sed(n) = GLB_AVEPERIO(n)/30 !FIXME: Do we want to use it?
      diagmon_sed(n) = .false.
      diagann_sed(n) = .false.
      diagdec_sed(n) = .false. ! isn't actually used now
      diagcen_sed(n) = .false. ! isn't actually used now
      diagmil_sed(n) = .false. ! isn't actually used now
      if (GLB_AVEPERIO(n) == 30) then
         diagmon_sed(n) = .true.
      elseif (GLB_AVEPERIO(n) == 365) then
         diagann_sed(n) = .true.
      elseif (GLB_AVEPERIO(n) == 3650) then
         diagdec_sed(n) = .true.
      elseif (GLB_AVEPERIO(n) == 36500) then
         diagcen_sed(n) = .true.
      elseif (GLB_AVEPERIO(n) == 365000) then
         diagmil_sed(n) = .true.
      endif
      if (GLB_FILEFREQ(n) < 0) then
         filefq_sed(n)=-real(nstepinday)/GLB_FILEFREQ(n)
      else
         filefq_sed(n)=nstepinday*max(1,GLB_FILEFREQ(n))
      endif
      filemon_sed(n) = .false.
      fileann_sed(n) = .false.
      filedec_sed(n) = .false.
      filecen_sed(n) = .false.
      filemil_sed(n) = .false.
      if (GLB_FILEFREQ(n) == 30) then
         filemon_sed(n) = .true.
      elseif (GLB_FILEFREQ(n) == 365) then
         fileann_sed(n) = .true.
      elseif (GLB_FILEPERIO(n) == 3650) then
         filedec_sed(n) = .true.
      elseif (GLB_FILEPERIO(n) == 36500) then
         filecen_sed(n) = .true.
      elseif (GLB_FILEPERIO(n) == 365000) then
         filemil_sed(n) = .true.
      endif
   enddo

   i_bsc_sed = 0
   do n = 1,nsed
      if (SDM_POWAIC(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowaic(n) = i_bsc_sed*min(1,SDM_POWAIC(n))
      if (SDM_POWAAL(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowaal(n) = i_bsc_sed*min(1,SDM_POWAAL(n))
      if (SDM_POWAPH(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowaph(n) = i_bsc_sed*min(1,SDM_POWAPH(n))
      if (SDM_POWAOX(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowaox(n) = i_bsc_sed*min(1,SDM_POWAOX(n))
      if (SDM_POWN2(n)  > 0) i_bsc_sed = i_bsc_sed+1
      jpown2(n)  = i_bsc_sed*min(1,SDM_POWN2(n))
      if (SDM_POWNO3(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowno3(n) = i_bsc_sed*min(1,SDM_POWNO3(n))
      if (SDM_POWASI(n) > 0) i_bsc_sed = i_bsc_sed+1
      jpowasi(n) = i_bsc_sed*min(1,SDM_POWASI(n))
      if (SDM_SSSO12(n) > 0) i_bsc_sed = i_bsc_sed+1
      jssso12(n) = i_bsc_sed*min(1,SDM_SSSO12(n))
      if (SDM_SSSSIL(n) > 0) i_bsc_sed = i_bsc_sed+1
      jssssil(n) = i_bsc_sed*min(1,SDM_SSSSIL(n))
      if (SDM_SSSC12(n) > 0) i_bsc_sed = i_bsc_sed+1
      jsssc12(n) = i_bsc_sed*min(1,SDM_SSSC12(n))
      if (SDM_SSSTER(n) > 0) i_bsc_sed = i_bsc_sed+1
      jssster(n) = i_bsc_sed*min(1,SDM_SSSTER(n))
   enddo

   i_bsc_bur = 0
   do n = 1,nsed
      if (BUR_SSSO12(n) > 0) i_bsc_bur = i_bsc_bur+1
      jburssso12(n) = i_bsc_bur*min(1,BUR_SSSO12(n))
      if (BUR_SSSC12(n) > 0) i_bsc_bur = i_bsc_bur+1
      jbursssc12(n) = i_bsc_bur*min(1,BUR_SSSC12(n))
      if (BUR_SSSSIL(n) > 0) i_bsc_bur = i_bsc_bur+1
      jburssssil(n) = i_bsc_bur*min(1,BUR_SSSSIL(n))
      if (BUR_SSSTER(n) > 0) i_bsc_bur = i_bsc_bur+1
      jburssster(n) = i_bsc_bur*min(1,BUR_SSSTER(n))
   enddo

   nsedt_sed  = i_bsc_sed
   nsedt_bur  = i_bsc_bur

   if (mnproc == 1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable sedtsed ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',ks
        write(io_stdo_bgc,*)'Forth dimension    : ',nsedt_sed
   endif

   allocate (sedt_sed(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,ks,          &
      &      nsedt_sed),stat=errstat)
   if (errstat /= 0) STOP 'not enough memory sedt_sed'
   if (nsedt_sed /= 0) sedt_sed = 0.

   if (mnproc == 1) THEN
        write(io_stdo_bgc,*)'Memory allocation for variable sedtbur ...'
        write(io_stdo_bgc,*)'First dimension    : ',kpie
        write(io_stdo_bgc,*)'Second dimension   : ',kpje
        write(io_stdo_bgc,*)'Third dimension    : ',nsedt_bur
   endif

   allocate (sedt_bur(1-nbdy:kpie+nbdy,1-nbdy:kpje+nbdy,             &
      &      nsedt_bur),stat=errstat)
   if (errstat /= 0) STOP 'not enough memory sedt_sed'
   if (nsedt_bur /= 0) sedt_bur = 0.

   !-- In HAMOCC, above is done by ALLOC_MEM_BGCMEAN() below by ALLOC_MEM_SEDMNT()

   if (mnproc == 1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'Memory allocation for spin-up-specific sediment modules:'
      write(io_stdo_bgc,*)' '
   endif

   if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variables ocetra_kbo_avg and ocetra_kbo_clim ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',nocetra
      write(io_stdo_bgc,*)'Fourth dimension   : ',12
   endif

   allocate (ocetra_kbo_avg(kpie,kpje,nocetra),stat=errstat)
   if(errstat /= 0) stop 'not enough memory ocetra_kbo_avg'
   allocate (ocetra_kbo_clim(kpie,kpje,nocetra,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory ocetra_kbo_clim'
   ocetra_kbo_avg(:,:,:) = 0.0
   ocetra_kbo_clim(:,:,:,:) = 0.0

   if (mnproc == 1) then
      WRITE(io_stdo_bgc,*)'Memory allocation for flux climatologies (*_clim) ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',12
   endif

   allocate (prorca_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory prorca_avg'
   allocate (prcaca_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory prcaca_avg'
   allocate (silpro_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory silpro_avg'
   allocate (produs_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory produs_avg'
   prorca_avg(:,:) = 0.0
   prcaca_avg(:,:) = 0.0
   silpro_avg(:,:) = 0.0
   produs_avg(:,:) = 0.0

   allocate (prorca_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory prorca_clim'
   allocate (prcaca_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory prcaca_clim'
   allocate (silpro_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory silpro_clim'
   allocate (produs_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory produs_clim'
   prorca_clim(:,:,:) = 0.0
   prcaca_clim(:,:,:) = 0.0
   silpro_clim(:,:,:) = 0.0
   produs_clim(:,:,:) = 0.0

   if (mnproc == 1) then
      WRITE(io_stdo_bgc,*)'Memory allocation for variables bolay_* ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',12
   endif

   allocate (bolay_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bolay_avg'
   allocate (bolay_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bolay_clim'
   bolay_avg(:,:) = 4000.0 ! maximum of lowest bottom layer thickness
   bolay_clim(:,:,:) = 0.0

   ALLOCATE (bgc_t_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_t_kbo_avg'
   ALLOCATE (bgc_t_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_t_kbo_clim'
   bgc_t_kbo_avg(:,:) = 0.0
   bgc_t_kbo_clim(:,:,:) = 0.0

   allocate (bgc_s_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_s_kbo_avg'
   allocate (bgc_s_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_s_kbo_clim'
   bgc_s_kbo_avg(:,:) = 0.0
   bgc_s_kbo_clim(:,:,:) = 0.0

   allocate (bgc_rho_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_rho_kbo_avg'
   allocate (bgc_rho_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory bgc_rho_kbo_clim'
   bgc_rho_kbo_avg(:,:) = 0.0
   bgc_rho_kbo_clim(:,:,:) = 0.0

   allocate (co3_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory co3_kbo_avg'
   allocate (co3_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory co3_kbo_clim'
   co3_kbo_avg(:,:) = 0.0
   co3_kbo_clim(:,:,:) = 0.0

   if (mnproc == 1) then
      write(io_stdo_bgc,*)'Memory allocation for variables keqb_* ...'
      write(io_stdo_bgc,*)'First dimension    : ',11
      write(io_stdo_bgc,*)'Second dimension   : ',kpie
      write(io_stdo_bgc,*)'Third dimension    : ',kpje
      write(io_stdo_bgc,*)'Fourth dimension   : ',12
   endif

   allocate (keqb_avg(11,kpie,kpje),stat=errstat)
   if(errstat /= 0) stop 'not enough memory keqb_avg'
   allocate (keqb_clim(11,kpie,kpje,12),stat=errstat)
   if(errstat /= 0) stop 'not enough memory keqb_clim'
   keqb_avg(:,:,:) = 0.0
   keqb_clim(:,:,:,:) = 0.0

end subroutine alloc_mem_sedmnt_offline

end module mo_sedmnt_offline
#endif
