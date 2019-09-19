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
        if ((((diagann_bgc(l).or.diagdec_bgc(l).and.mod(nyear,10)==0.or.
     .    diagcen_bgc(l).and.mod(nyear,100)==0.or.diagmil_bgc(l).and.
     .    mod(nyear,1000)==0).and.nday_of_year.eq.nday_in_year.or.diagmon_bgc(l)
     .    .and.nday.eq.nd_in_m(nmonth)).and.mod(nstep,nstep_in_day).eq.0).or..not.
     .    (diagann_bgc(l).or.diagmon_bgc(l).or.diagdec_bgc(l).or.
     .    diagcen_bgc(l).or.diagmil_bgc(l)).and.
     .    mod(nstep+.5,diagfq_bgc(l)).lt.1.)
     .    bgcwrt(l) = .true.
      enddo
c
c     At the start of every year, and when not reading climatology (or
c     at the start of the simulation for startup which begins at day 2),
c     increase the simulation year nyear_global.
      if ( (nday_of_year==1 .and. mod(ldtrunbgc,nstep_in_day)==0 .or.
     .   ldtrunbgc==0)
#ifdef SED_OFFLINE
     .   .and. .not. lread_clim
#endif
     .   ) then
         nyear_global = nyear_global + 1
         if (mnproc == 1) write(io_stdo_bgc,'(a,i6)')
     .            'hamocc_step(): nyear        =', nyear
         if (mnproc.eq.1) write(io_stdo_bgc,'(a,i6)')
     .            'hamocc_step(): nyear_global =', nyear_global
      endif
c
#ifdef SED_OFFLINE
c     Before running full HAMOCC, the offline routine must be called
c     (including optional reading of climatology and check if ready).
      call sedmnt_offline(idm, jdm, kdm, maxyear_sediment, nstep,
     .                    bgc_dp, bgc_dx, bgc_dy, omask)
#endif
c
      call hamocc4bcm(idm,jdm,kdm,pglat,bgc_swr,bgc_fice,
     .            bgc_t,bgc_s,bgc_slp,bgc_rho,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,nday_in_year,bgc_flxdms)
c
#ifdef SED_OFFLINE
c     Accumulate bottom water tracers (obtained from HAMOCC integration)
c     when in the last year of MICOM/HAMOCC simulation.
      if ( (lsed_wclim .or. lsed_spinup) .and. mod(nyear,maxyear_ocean)==0 ) then
         call prepare_clim(nstep)
      endif
#endif
c
      call hamocc2micom(m,n,mm,nn)
c
      return
      end
