      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- perform one HAMOCC step     
c --- depending on preprocessor directives, run a sediment spin-up (MvH)
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mo_bgcmean
      use mo_control_bgc, only: nyear_global, is_end_of_day, io_stdo_bgc
#if defined(SED_OFFLINE)
     .                        , maxyear_ocean
#endif
#if defined(SED_OFFLINE) || defined(SED_RCLIM) || defined(SED_WCLIM)
     .                        , rmasko, nstep_in_month
     .                        , must_read_clim
     .                        , lrunsed, imonth
      use mo_carbch,      only: ocetra, keqb, co3
      use mo_param1_bgc,  only: nocetra, idet, icalc, iopal, ifdust
      use mo_biomod,      only: kbo, bolay, wpoc, wmin, wmax, wdust
     .                        , wlin, wopal, wcal
      use mo_sedmnt
#endif
c
      implicit none
c
#include "common_bgc.h"
#include "common_blocks.h"
#include "common_clndr.h"
c
      integer, intent(in) :: m,n,mm,nn,k1m,k1n
      integer :: l
#if defined(SED_OFFLINE) || defined(SED_RCLIM) || defined(SED_WCLIM)
      integer :: i,j,iocetra
      character (len = *), parameter :: bscfnm = "bottom_seawater_clim.nc"
      integer :: nday_of_month                   ! current day of month
#endif
c
      call trc_limitc(nn)
c
      call micom2hamocc(m,n,mm,nn)
c
      do l=1,nbgc
        bgcwrt(l)=0.
        if (((diagann_bgc(l).and.nday_of_year.eq.1.or.diagmon_bgc(l)
     .    .and.nday.eq.1).and.mod(nstep,nstep_in_day).eq.0).or..not.
     .    (diagann_bgc(l).or.diagmon_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.1.)
     .    bgcwrt(l)=1.
      enddo
c
c --- Apparently, nstep=0 (mod 18) at the start of a day!
      is_end_of_day = mod(nstep+1,nstep_in_day)==0
      if ( nday_of_year==nday_in_year .and. is_end_of_day ) then
         nyear_global = nyear_global + 1
         IF (mnproc.eq.1) WRITE(io_stdo_bgc,*)
     .            'hamocc_step(): nyear_global = ', nyear_global
      endif
c

#if defined(SED_OFFLINE) || defined(SED_RCLIM) || defined(SED_WCLIM)
c --- ------------------------------------------------------------------
c --- Sediment spin-up code - MvH
c --- 
c --- Depending on the pre-processing keys, read or calculate (and
c --- write to disk) a bottom seawater array needed as a forcing for the
c --- subsequent sediment() spin-up routine.
c --- ------------------------------------------------------------------

c --- Default action is not to run sediment()
      lrunsed = .false.

# if defined(SED_RCLIM) && defined(SED_WCLIM)
      IF( mnproc==1 ) WRITE(io_stdo_bgc,*)
     .   'hamocc_step(): WARNING: Bottom-water climatology will be read
     . from and written to the same file.'
# endif

c# if defined(SED_RCLIM)
c --- Read the bottom seawater climatology from the netCDF forcing file
      if (must_read_clim) then
      IF (mnproc.eq.1) WRITE(io_stdo_bgc,*)
     .         'hamocc_step(): starting bottom seawater climatology read loop'
      do imonth = 1,12
         call aufr_bgc_onlysed(idm,jdm,kdm,nyear,imonth,nday,nstep
     .                        ,bscfnm)
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
      else
c# else
c --- Accumulate the bottom seawater fields from HAMOCC
      nstep_in_month  = nstep_in_month + 1
      do iocetra = 1, nocetra
         do j = 1, jdm
            do i = 1, idm
               ocetra_kbo_avg(i,j,iocetra) = ocetra_kbo_avg(i,j, iocetra)
     .                                     + ocetra(i,j,kbo(i,j),iocetra)
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

c --- Calculate the day of month
      nday_of_month = nday_of_year
      do i = 1, nmonth-1
         nday_of_month = nday_of_month - nd_in_m(i)
      enddo

c --- Calculate tracer monthly average
      if ( nday_of_month==nd_in_m(nmonth) .and. is_end_of_day ) then
         IF (mnproc.eq.1) WRITE(io_stdo_bgc,*)
     .      'hamocc_step(): end of month, set tracer avg for last month'
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
#  if defined(SED_WCLIM)
c --- -- Write tracer monthly averages to netCDF file
         call aufw_bgc_onlysed(idm,jdm,kdm,nyear,nmonth,nday,nstep
     .                        ,bscfnm)
#  endif
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
      endif ! end of month routines
      endif
c# endif

# if defined(SED_OFFLINE)
c --- The whole sediment() routine is done as soon as the needed bottom
c --- water variables are available, at two subsequent MICOM timesteps,
c --- following the leapfrog scheme

c --- start immediately with sediment() if bottom water climatology is read,
c --- else wait until the end of the year when we have bottom-water fields
      if (must_read_clim) then
         lrunsed = .true.
      elseif (maxyear_ocean<=0) then
         lrunsed = .false.
      else
         lrunsed = mod(nyear,maxyear_ocean)==0 .and. nday_of_year==nday_in_year
     .             .and. is_end_of_day
      endif

      if ( lrunsed ) then
         IF (mnproc.eq.1) WRITE(io_stdo_bgc,*)
     .            'hamocc_step(): sediment spin-up starting'
         call sediment(idm,jdm,kdm,pglat,
     .            bgc_dp,bgc_dx,bgc_dy,
     .            omask)
         IF (mnproc.eq.1) WRITE(io_stdo_bgc,*)
     .            'hamocc_step(): sediment spin-up ended'
      endif
      lrunsed = .false.
# endif
      must_read_clim = .false.

#endif
      call hamocc4bcm(idm,jdm,kdm,pglat,bgc_swr,bgc_fice,
     .            bgc_t,bgc_s,bgc_slp,bgc_rho,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,nday_in_year,bgc_flxdms)

      call hamocc2micom(m,n,mm,nn)
c
      return
      end
