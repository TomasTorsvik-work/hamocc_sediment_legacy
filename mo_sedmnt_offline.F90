#if defined(SED_OFFLINE)
module mo_sedmnt_offline

!-----------------------------------------------------------------------
!
! Routines for offline sediment spin-up
!
! Copyright (C) 2018 Marco van Hulten <Marco.Hulten@uib.no>
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
!  Depending on the truth values of lsed_*, this subroutine will
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
! - subroutine sed_offline()
!     Main off-line sediment routine.
!     If lsed_rclim is set (namelist), it will read_clim(), and then
!     spin up the sediment by looping sediment_step() over maxyear_sediment.
!     In case no climatology was read, it will first collect monthly
!     averages of bottom water layer particle fluxes and dissolved
!     concentrations by running HAMOCC for one year, writing the monthly
!     averages to file through write_clim().
!
! - subroutine alloc_mem_sedmnt_offline()
!     Allocate memory for bottom-water climatology, analogous to mo_sedmnt's
!     alloc_mem_sedmnt().
!
! Refer to the flowchart in the documentation.
!
! History:
!
! 2018-06   Initial version by MvH.
! 2018-07   Rewrite.
!
!-----------------------------------------------------------------------

use mo_control_bgc, only: nyear_global, is_end_of_day, io_stdo_bgc   &
   &                    , maxyear_ocean                              &
   &                    , rmasko, nstep_in_month                     &
   &                    , lsed_rclim, lsed_wclim, lsed_spinup        &
   &                    , lread_clim, lwrite_clim, ldo_spinup        &
   &                    , imonth, iyear
use mo_carbch,      only: ocetra, keqb, co3
use mo_param1_bgc,  only: nocetra, idet, icalc, iopal, ifdust
use mo_biomod,      only: kbo, bolay, wpoc, wmin, wmax, wdust        &
   &                    , wlin, wopal, wcal
use mod_xc  !, only: ii,jj,kk,idm,jdm,kdm,nbdy,ifp,isp,ilp
use mo_common_bgc
use mo_common_bgcs
use mo_bgcmean!,      only: nacc_bgc


implicit none

!#include "common_bgc.h90"
!#include "common_blocks.h"
#include "common_clndr.h90"
!#include "common_bgcs.h90"

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

! private variables
!
integer, private :: i,j,iocetra
character(len = *), parameter, private :: bscfnm = "bottom_seawater_clim.nc"
integer, private :: nday_of_month ! current day of month

contains

subroutine sed_init()
   ! Default action is not to do the sediment spin-up
   ldo_spinup = .false.
end subroutine

subroutine read_clim()
   ! Read the bottom seawater climatology from the netCDF forcing file
   if (mnproc.eq.1) write(io_stdo_bgc,*)                             &
      &        'hamocc_step(): starting bottom seawater climatology read loop'
   do imonth = 1,12
      call aufr_bgc_onlysed(idm,jdm,kdm,nyear,imonth,nday,nstep,bscfnm)
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

subroutine write_clim()
   ! Write tracer monthly averages to netCDF file
   call aufw_bgc_onlysed(idm,jdm,kdm,nyear,nmonth,nday,nstep,bscfnm)
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
end subroutine write_clim

subroutine sed_offline(kpie, kpje, kpke, maxyear,                    &
                &      pglat, pddpo, pdlxp, pdlyp, omask)
   ! FIXME: bit ugly:
   use mo_sedmnt, only: sedlay, powtra, burial

   ! Subprogram arguments
   !
   integer, intent(in)                          :: kpie,kpje,kpke
   integer, intent(in)                          :: maxyear
   real, dimension(kpie,kpje), intent(in)       :: pglat
   real, dimension(kpie,kpje,kpke), intent(in)  :: pddpo
   real, dimension(kpie,kpje), intent(in)       :: pdlxp, pdlyp
   real, dimension(kpie,kpje), intent(in)       :: omask

   integer  :: nday_save, nmonth_save   ! calendar variables outside sediment()
   integer  :: nyear_save, nday_of_year_save, nd_in_m_2, nday_in_year_save

   if (lread_clim) then
      call read_clim()
      lread_clim = .false.
   else
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
      do i = 1, nmonth-1
         nday_of_month = nday_of_month - nd_in_m(i)
      enddo

      ! Calculate tracer monthly average
      if ( nday_of_month==nd_in_m(nmonth) .and. is_end_of_day ) then
         if (mnproc.eq.1) write(io_stdo_bgc,*)                             &
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
         if (lsed_wclim) call write_clim()
      endif ! end of month routines
   endif ! (if .not. lsed_rclim)

   if (lsed_spinup) then
   ! The time loop over the sediment routines is done as soon as the needed bottom
   ! water variables are available

   ! start immediately with sediment() if bottom water climatology is read,
   ! else wait until the end of the year when we have bottom-water fields
   if (lsed_rclim) then
      ldo_spinup = .true.
   elseif (maxyear_ocean<=0) then
      ldo_spinup = .false.
   else
      ldo_spinup = mod(nyear,maxyear_ocean)==0 .and. nday_of_year==nday_in_year  &
         &      .and. is_end_of_day
   endif

   if ( ldo_spinup ) then
      ! set up sediment layers (mainly for much higher diffusion rate)
      call bodensed_onlysed(kpie,kpje,kpke,pddpo)

      ! save calendar variables
      nday_save         = nday
      nmonth_save       = nmonth
      nyear_save        = nyear
      nday_of_year_save = nday_of_year
      nd_in_m_2         = nd_in_m(2)
      nday_in_year_save = nday_in_year

      do iyear = 1, maxyear
         nyear_global = nyear_global + 1
         if (mnproc.eq.1) write(io_stdo_bgc,*) 'sediment(): nyear_global = ', nyear_global
         do imonth = 1, 12
            if (mnproc.eq.1) write(io_stdo_bgc,*)                             &
               &     'hamocc_step(): sediment spin-up starting'
            call updcln_onlysed() ! do a monthly calendar update only for sediment_step()
            call sediment_step(idm,jdm,kdm,pglat, bgc_dp,bgc_dx,bgc_dy,    &
               & bgc_s_kbo_clim(:,:,imonth), bgc_rho_kbo_clim(:,:,imonth), &
               & ocetra_kbo_clim(:,:,imonth,:), bolay_clim(:,:,imonth),    &
               & keqb_clim(:,:,:,imonth),                                  &
               & prorca_clim(:,:,imonth), prcaca_clim(:,:,imonth),         &
               & silpro_clim(:,:,imonth), produs_clim(:,:,imonth),         &
               & co3_kbo_clim(:,:,imonth))
            ! write monthly outputs (assuming first index is mo)
            if (maxyear <= 20) then
               nacc_bgc(1) = 1
               if (GLB_INVENTORY(1).ne.0)                                        &
                  &  CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
               call ncwrt_bgc(1) ! ncout_hamocc.F
               nacc_bgc(1) = 0
            endif
         enddo
         ! write yearly outputs (assuming second index is yr)
         nacc_bgc(2) = 1
         if (GLB_INVENTORY(2).ne.0)                                        &
            &  CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,0)
         imonth = 0 ! triggers writing of yearly averages
         call ncwrt_bgc(2) ! ncout_hamocc.F
         nacc_bgc(2) = 0
      enddo

      ! update both time levels before returning to full model (leapfrog scheme)
      sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
      sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:)
      powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
      powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:)
      burial2(:,:,1,:)         = burial(:,:,:)
      burial2(:,:,2,:)         = burial(:,:,:)

      ! set calendar variables back to original values
      nday         = nday_save
      nmonth       = nmonth_save
      nyear        = nyear_save
      nday_of_year = nday_of_year_save
      nd_in_m(2)   = nd_in_m_2
      nday_in_year = nday_in_year_save

      if (mnproc.eq.1) write(io_stdo_bgc,*)                             &
         &     'hamocc_step(): sediment spin-up ended'
      endif
      ldo_spinup = .false.
   endif ! spin-up
end subroutine sed_offline

subroutine alloc_mem_sedmnt_offline(kpie, kpje)
   use mod_xc
   use mo_control_bgc
   use mo_param1_bgc

   integer :: kpie, kpje
   integer :: errstat

   if (mnproc.eq.1) then
      write(io_stdo_bgc,*)' '
      write(io_stdo_bgc,*)'Memory allocation for spin-up-specific sediment modules:'
      write(io_stdo_bgc,*)' '
   endif

   if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variables ocetra_kbo_avg and ocetra_kbo_clim ...'
      write(io_stdo_bgc,*)'First dimension    : ',kpie
      write(io_stdo_bgc,*)'Second dimension   : ',kpje
      write(io_stdo_bgc,*)'Third dimension    : ',nocetra
      write(io_stdo_bgc,*)'Fourth dimension   : ',12
   endif

   allocate (ocetra_kbo_avg(kpie,kpje,nocetra),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory ocetra_kbo_avg'
   allocate (ocetra_kbo_clim(kpie,kpje,nocetra,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory ocetra_kbo_clim'
   ocetra_kbo_avg(:,:,:) = 0.0
   ocetra_kbo_clim(:,:,:,:) = 0.0

   if (mnproc.eq.1) then
      WRITE(io_stdo_bgc,*)'Memory allocation for flux climatologies (*_clim) ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',12
   endif

   allocate (prorca_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory prorca_avg'
   allocate (prcaca_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory prcaca_avg'
   allocate (silpro_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory silpro_avg'
   allocate (produs_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory produs_avg'
   prorca_avg(:,:) = 0.0
   prcaca_avg(:,:) = 0.0
   silpro_avg(:,:) = 0.0
   produs_avg(:,:) = 0.0

   allocate (prorca_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory prorca_clim'
   allocate (prcaca_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory prcaca_clim'
   allocate (silpro_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory silpro_clim'
   allocate (produs_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory produs_clim'
   prorca_clim(:,:,:) = 0.0
   prcaca_clim(:,:,:) = 0.0
   silpro_clim(:,:,:) = 0.0
   produs_clim(:,:,:) = 0.0

   if (mnproc.eq.1) then
      WRITE(io_stdo_bgc,*)'Memory allocation for variables bolay_* ...'
      WRITE(io_stdo_bgc,*)'First dimension    : ',kpie
      WRITE(io_stdo_bgc,*)'Second dimension   : ',kpje
      WRITE(io_stdo_bgc,*)'Third dimension    : ',12
   endif

   allocate (bolay_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bolay_avg'
   allocate (bolay_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bolay_clim'
   bolay_avg(:,:) = 4000.0 ! maximum of lowest bottom layer thickness
   bolay_clim(:,:,:) = 0.0

   ALLOCATE (bgc_t_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_t_kbo_avg'
   ALLOCATE (bgc_t_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_t_kbo_clim'
   bgc_t_kbo_avg(:,:) = 0.0
   bgc_t_kbo_clim(:,:,:) = 0.0

   allocate (bgc_s_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_s_kbo_avg'
   allocate (bgc_s_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_s_kbo_clim'
   bgc_s_kbo_avg(:,:) = 0.0
   bgc_s_kbo_clim(:,:,:) = 0.0

   allocate (bgc_rho_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_rho_kbo_avg'
   allocate (bgc_rho_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory bgc_rho_kbo_clim'
   bgc_rho_kbo_avg(:,:) = 0.0
   bgc_rho_kbo_clim(:,:,:) = 0.0

   allocate (co3_kbo_avg(kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory co3_kbo_avg'
   allocate (co3_kbo_clim(kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory co3_kbo_clim'
   co3_kbo_avg(:,:) = 0.0
   co3_kbo_clim(:,:,:) = 0.0

   if (mnproc.eq.1) then
      write(io_stdo_bgc,*)'Memory allocation for variables keqb_* ...'
      write(io_stdo_bgc,*)'First dimension    : ',11
      write(io_stdo_bgc,*)'Second dimension   : ',kpie
      write(io_stdo_bgc,*)'Third dimension    : ',kpje
      write(io_stdo_bgc,*)'Fourth dimension   : ',12
   endif

   allocate (keqb_avg(11,kpie,kpje),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory keqb_avg'
   allocate (keqb_clim(11,kpie,kpje,12),stat=errstat)
   if(errstat.ne.0) stop 'not enough memory keqb_clim'
   keqb_avg(:,:,:) = 0.0
   keqb_clim(:,:,:,:) = 0.0

end subroutine alloc_mem_sedmnt_offline

end module mo_sedmnt_offline
#endif
