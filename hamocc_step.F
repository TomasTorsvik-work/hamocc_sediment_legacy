      subroutine hamocc_step(m,n,mm,nn,k1m,k1n)
c
c --- ------------------------------------------------------------------
c --- perform one HAMOCC step     
c --- depending on preprocessor directives, run a sediment spin-up (MvH)
c --- ------------------------------------------------------------------
c
      use mod_xc
      use mo_bgcmean
      use mo_control_bgc
#ifdef SED_OFFLINE
      use mo_sedmnt_offline, only: sedmnt_offline, prepare_clim ! subroutines
#endif
      use mo_common_bgc
c
      implicit none
c
#include "common_blocks.h"
#include "common_clndr.h"
c
      integer, intent(in) :: m,n,mm,nn,k1m,k1n
      integer :: l
c
      call trc_limitc(nn)
c
      call micom2hamocc(m,n,mm,nn)
c
      do l=1,nbgc
        bgcwrt(l) = .false.
        if (((diagann_bgc(l).and.nday_of_year.eq.1.or.diagmon_bgc(l)
     .    .and.nday.eq.1).and.mod(nstep,nstep_in_day).eq.0).or..not.
     .    (diagann_bgc(l).or.diagmon_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.1.)
     .    bgcwrt(l) = .true.
      enddo
c
      is_start_of_day = mod(nstep-1,nstep_in_day)==0
#ifdef SED_OFFLINE
      call sedmnt_offline(idm, jdm, kdm, maxyear_sediment, nstep,
     .                    bgc_dp, bgc_dx, bgc_dy, omask)
#endif

      ! January 2 to make it work with MICOM clock for a startup run
      if ( nday_of_year==2 .and. is_start_of_day ) then
         nyear_global = nyear_global + 1
         if (mnproc.eq.1) write(io_stdo_bgc,'(a,i6)')
     .            'hamocc_step(): nyear_global = ', nyear_global
      endif
c
      call hamocc4bcm(idm,jdm,kdm,pglat,bgc_swr,bgc_fice,
     .            bgc_t,bgc_s,bgc_slp,bgc_rho,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,nday_in_year,bgc_flxdms)

      is_end_of_day = mod(nstep,nstep_in_day)==0
#ifdef SED_OFFLINE
      ! accumulate bottom water tracers when in
      ! the last year of MICOM/HAMOCC simulation
      if ( (lsed_wclim .or. lsed_spinup) .and. mod(nyear,maxyear_ocean)==0
     .     .and. .not. nyear_global==0 ) then
         call prepare_clim(nstep)
      endif
#endif
c
      call hamocc2micom(m,n,mm,nn)
c
      return
      end
