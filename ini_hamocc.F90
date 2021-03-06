      SUBROUTINE INI_HAMOCC(kpaufr,kpicycli,pdt,kpndtrun,kpie,kpje,kpke  &
     &            ,kpbe,pddpo,ptho,psao,prho,pdlxp,pdlyp,ptiestu,ptiestw &
     &            ,kplyear,kplmonth,kplday,kpldtoce                      &
     &            ,pglon,pglat,omask,ntr,ntrbgc,itrbgc                   &
     &            ,trc,sedlay2,powtra2,burial2                           &    
     &            ,rstfnm_ocn,path,path2)

!****************************************************************
!
!**** *INI_BGC* - initialize marine bio-geo-chemistry module.
!
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Modified
!     --------
!     J.Schwinger       *GFI, Bergen*    2013-10-21
!     - added GNEWS2 option for riverine input of carbon and nutrients
!     - code cleanup
!     J.Schwinger,      *GFI, Bergen*    2014-05-21
!     - adapted code for use with two time level tracer field in MICOM
!     
!     Purpose
!     -------
!     - allocate and initialize bgc variables
!     - allocate and initialize sediment layering
!     - initialize bgc constants (beleg.F90)
!     - read restart fields if kpaufr=1
!
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpaufr*     - 1/0 for read / do not read restart file
!     *INTEGER* *kpicycli*   - flag for cyclicity.
!     *REAL*    *pdt*        - ocean model time step [sec].
!     *INTEGER* *kpndtrun*   - total no. of time steps of run.
!     *INTEGER* *kpie*       - zonal dimension of model grid.
!     *INTEGER* *kpje*       - meridional dimension of model grid.
!     *INTEGER* *kpke*       - vertical dimension of model grid.
!     *REAL*    *pddpo*      - size of scalar grid cell (vertical dimension) [m].
!     *REAL*    *ptho*       - potential temperature [deg C].
!     *REAL*    *psao*       - salinity [psu].
!     *REAL*    *prho*       - density [g/cm^3].
!     *REAL*    *pdlxp*      - size of scalar grid cell (zonal) [m].
!     *REAL*    *pdlyp*      - size of scalar grid cell (meridional) [m].
!     *REAL*    *ptiestu*    - depth of level [m].
!     *REAL*    *ptiestw*    - depth of level interface [m].
!     *INTEGER* *kplyear*    - year  in ocean restart date
!     *INTEGER* *kplmonth*   - month in ocean restart date
!     *INTEGER* *kplday*     - day   in ocean restart date
!     *INTEGER* *kpldtoce*   - step  in ocean restart date
!     *REAL*    *pglon*      - geographical longitude of grid points [degree E].
!     *REAL*    *pglat*      - geographical latitude  of grid points [degree N].
!     *REAL*    *omask*      - land/ocean mask
!     *INTEGER* *ntr*        - number of tracers in tracer field
!     *INTEGER* *ntrbgc*     - number of biogechemical tracers in tracer field
!     *INTEGER* *itrbgc*     - start index for biogeochemical tracers in tracer field
!     *REAL*    *trc*        - initial/restart tracer field to be passed to the 
!                              ocean model [mol/kg]
!     *REAL*    *sedlay2*    - initial/restart sediment (two time levels) field
!     *REAL*    *powtra2*    - initial/restart pore water tracer (two time levels) field
!     *REAL*    *burial2*    - initial/restart sediment burial (two time levels) field
!     *CHAR*    *rstfnm_ocn* - restart file name-informations
!     *CHAR*    *path*       - path to data files
!     *CHAR*    *path2*      - path to restart files
!
!**********************************************************************

      USE mo_carbch
      USE mo_sedmnt
#if defined(SED_OFFLINE)
      USE mo_sedmnt_offline, only: alloc_mem_sedmnt_offline
#endif
      USE mo_biomod
      USE mo_control_bgc
      use mo_param1_bgc 
      use mod_xc, only: mnproc,lp,nfu
      use mo_bgcmean
#ifdef RIV_GNEWS
      use mo_riverinpt
#endif
 
      implicit none

#include "common_clndr.h90"

      INTEGER :: kpie,kpje,kpke,kpbe,ntr,ntrbgc,itrbgc
      INTEGER :: kplyear,kplmonth,kplday,kpldtoce
      INTEGER :: kpaufr,kpicycli,kpndtrun,k,l
      INTEGER :: i,j
      
      REAL    :: pddpo(kpie,kpje,kpke)
      REAL    :: ptho (kpie,kpje,kpke)
      REAL    :: psao (kpie,kpje,kpke)
      REAL    :: prho (kpie,kpje,kpke)
      REAL    :: pdlxp(kpie,kpje),pdlyp(kpie,kpje)
      REAL    :: pglon(kpie,kpje),pglat(kpie,kpje)
      REAL    :: ptiestu(kpie,kpje,kpke+1),ptiestw(kpie,kpje,kpke+1)
      REAL    :: omask(kpie,kpje)
      REAL    :: trc(1-kpbe:kpie+kpbe,1-kpbe:kpje+kpbe,2*kpke,ntr)
      REAL    :: sedlay2(kpie,kpje,2*ks,nsedtra)
      REAL    :: powtra2(kpie,kpje,2*ks,npowtra)
      REAL    :: burial2(kpie,kpje,2,   nsedtra)
      REAL    :: pdt
      character(len=*) :: rstfnm_ocn,path,path2

! Define io units

      io_stdo_bgc = lp      !  standard out.
      io_stdi_bgc = 5       !  standard in.
      io_rsti_bgc = nfu     !  restart in. 
      io_rsto_bgc = nfu     !  restart out.
      io_nml = nfu          !  namelist

      if (mnproc.eq.1) then
      write(io_stdo_bgc,*) 'HAMOCC initialisation'
      write(io_stdo_bgc,*) 'restart',kpaufr,kpicycli,kpndtrun
      write(io_stdo_bgc,*) 'dims',kpie,kpje,kpke
      write(io_stdo_bgc,*) 'time',kplyear,kplmonth,kplday,kpldtoce
      write(io_stdo_bgc,*) 'time step',pdt
      endif
!                    
! Set control constants ( mo_control_bgc )
!
      dtbgc = pdt                   !  time step length [sec].
      ndtdaybgc=NINT(86400./dtbgc)  !  time steps per day [No].
      dtb=1./ndtdaybgc              !  time step length [days].
      dtsed = dtbgc                 !  time step length for sediment [sec].

#if defined(SED_OFFLINE)
      if (lsed_rclim .and. .not. (kplmonth==1 .and. kplday<=2)) then
         write(io_stdo_bgc,*) 'WARNING: Not at start of year!  The transition between'
         write(io_stdo_bgc,*) 'stand-alone and coupled sediment will be inconsistent!'
      endif
#endif

      icyclibgc = kpicycli
      ndtrunbgc = kpndtrun

!
! Initialize some namelist parameters
!
      isac = 1

!
! Initialize time step counter of run.
!
      ldtrunbgc = 0
#if defined(SED_OFFLINE)
      nstep_in_month = 0
#endif
      nyear_global = 0
      if (mnproc == 1) write(io_stdo_bgc,'(a,i6)')                      &
            &         'ini_hamocc(): nyear_global = ', nyear_global

      CALL ALLOC_MEM_BGCMEAN(kpie,kpje,kpke)

!                        
! Allocate memory : biology
!
      CALL ALLOC_MEM_BIOMOD(kpie,kpje,kpke)

!                        
! Allocate memory : sediment
!
      CALL ALLOC_MEM_SEDMNT(kpie,kpje)
#if defined(SED_OFFLINE)
      call alloc_mem_sedmnt_offline(kpie, kpje)
#endif

!                        
! Allocate memory : inorganic carbon cycle
!
      CALL ALLOC_MEM_CARBCH(kpie,kpje,kpke)

!                        
! Initialize sediment layering
!
      CALL BODENSED(kpie,kpje,kpke,pddpo)

!                        
! Initialize sediment and ocean tracer.
! 
      CALL BELEG_BGC(kpie,kpje,kpke,pddpo,ptiestw,prho,                  &
     &               omask,pglon,pglat,path)
     
!
! Initialise river input
!
#ifdef RIV_GNEWS
      call INI_RIVERINPT(path)
#endif

!                        
! Read restart fields from restart file if requested, otherwise 
! (at first start-up) copy ocetra and sediment arrays (which are
! initialised in BELEG) to both timelevels of their respective
! two-time-level counterpart
!
      IF(kpaufr.eq.1) THEN
         CALL AUFR_BGC(kpie,kpje,kpke,ntr,ntrbgc,itrbgc,                 &
     &                 trc,sedlay2,powtra2,burial2,                      &
     &                 kplyear,kplmonth,kplday,kpldtoce,                 &
     &                 rstfnm_ocn,path2)
      ELSE
         trc(1:kpie,1:kpje,1:kpke,       itrbgc:itrbgc+ntrbgc-1) =       &
     &     ocetra(:,:,:,:)
         trc(1:kpie,1:kpje,kpke+1:2*kpke,itrbgc:itrbgc+ntrbgc-1) =       &
     &     ocetra(:,:,:,:)
         sedlay2(:,:,1:ks,:)      = sedlay(:,:,:,:)
         sedlay2(:,:,ks+1:2*ks,:) = sedlay(:,:,:,:) 
         powtra2(:,:,1:ks,:)      = powtra(:,:,:,:)
         powtra2(:,:,ks+1:2*ks,:) = powtra(:,:,:,:) 
         burial2(:,:,1,:)         = burial(:,:,:)
         burial2(:,:,2,:)         = burial(:,:,:) 
      ENDIF

! aufsetz! (for initialization of 14C)
!      call c14_correction(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,psao,         &
!     &                    ptho,omask)

!#ifdef DIFFAT
! correction of alkalinity during spin-up of pco2 
!
!      CALL SPINUP_BGC(kpie,kpje,kpke,omask,pdlxp,pdlyp)
!#endif

!
! Global inventory of all tracers
!      
!      CALL INVENTORY_BGC(kpie,kpje,kpke,pdlxp,pdlyp,pddpo,omask,1)


      RETURN
      END
