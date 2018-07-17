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
      use mo_sedmnt_offline, only: sed_offline ! subroutine
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
#ifdef SED_OFFLINE
      call sed_offline(idm, jdm, kdm, maxyear_sediment,
     .                 pglat, bgc_dp, bgc_dx, bgc_dy, omask)
#endif
c
      call hamocc4bcm(idm,jdm,kdm,pglat,bgc_swr,bgc_fice,
     .            bgc_t,bgc_s,bgc_slp,bgc_rho,bgc_dp,bgc_dx,bgc_dy,
     .            bgc_pu,bgc_pw,bgc_dpio,bgc_awnd,bgc_atmco2,bgc_flxco2,
     .            nyear,nmonth,nday,nd_in_m(nmonth),ldtmonth,ldtday,
     .            omask,nday_in_year,bgc_flxdms)

      call hamocc2micom(m,n,mm,nn)
c
      return
      end
