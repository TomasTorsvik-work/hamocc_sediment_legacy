      SUBROUTINE BELEG_BGC(kpie,kpje,kpke,pddpo,ptiestw,prho,  &
     &                     omask,pglon,pglat,path)
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/beleg_bgc.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\

!****************************************************************
!
!**** *BELEG_BGC* - initialize bgc variables.
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*     10.04.01
!
!     J.Schwinger,      *GFI, Bergen*     2013-10-21
!     - corrected units of tracer fields at initialisation to
!       mol/kg as tracers are passed to ocean model first,
!       unit conversion is done using in situ density
!     - code clean up
!
!     I. Kriest, GEOMAR, 11.08.2016
!     - included T-dependence of cyanobacteria growth
!     - modified stoichiometry for denitrification
! 
!     Purpose
!     -------
!     - set start values for  bgc variables.
!
!     Method
!     -------
!     - bgc variables are initialized. They might be overwritten
!       if read from restart by call of AUFR_BGC.
!     - physical constants are defined
!     - fields used for budget calculations (should be in extra SBR!)
!       and as boundary conditions are set to zero.
!     
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!     *INTEGER*   *kpie*    - 1st dimension of model grid.
!     *INTEGER*   *kpje*    - 2nd dimension of model grid.
!     *INTEGER*   *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*      *pddpo*   - size of grid cell (3rd dimension) [m].
!     *REAL*      *ptiestw* - depth of layer interfaces [m].
!     *REAL*      *prho*    - density [g/cm^3].
!     *REAL*      *omask*   - ocean mask.
!     *REAL*      *pglon*   - longitude of grid cell [deg].
!     *REAL*      *pglat*   - latitude  of grid cell [deg].
!     *CHARACTER* *path*    - path to input data files.
!
!
!**********************************************************************

      USE mo_carbch
      USE mo_biomod
      USE mo_sedmnt
      USE mo_control_bgc
      use mo_param1_bgc 
      USE mod_xc, only: mnproc


      implicit none      

      REAL :: pddpo(kpie,kpje,kpke)
      REAL :: ptiestw(kpie,kpje,kpke+1)
      REAL :: prho (kpie,kpje,kpke)
      REAL :: omask(kpie,kpje)
      REAL :: pglon(kpie,kpje)
      REAL :: pglat(kpie,kpje)
      REAL :: north, south
      INTEGER :: i,j,k,l,ii,jj,m,kpie,kpje,kpke
      character(len=*) :: path

#ifdef AGG
      REAL :: shear,snow
#endif 

#ifndef AGG
      REAL :: dustd1, dustd2, dustsink
#endif

      integer :: p_joff,p_ioff

#ifdef SED_OFFLINE
      lcompleted_clim  = .false.
      namelist /bgcnml/ atm_co2, maxyear_sediment, maxyear_ocean, nburst_last, &
         &              lsed_rclim, lsed_wclim, lsed_spinup
#else
      namelist /bgcnml/ atm_co2
#endif
      lspinning_up_sed = .false.

!
! Initialize overall time step counter.
!
      ldtbgc = 0
!
!
#ifndef DIFFAT            
!
! Obtain the atmospheric CO2 concentration and sediment spin-up parameters.
!
      open (unit=io_nml,file='ocn_in',status='old',action='read',      &
     &      recl=80)
      read (unit=io_nml,nml=BGCNML)
      close (unit=io_nml)
      IF (mnproc.eq.1) THEN
        write(io_stdo_bgc,*) 'HAMOCC: atmospheric co2:',atm_co2
#ifdef SED_OFFLINE
        write(io_stdo_bgc,*) 'HAMOCC: maxyear_sediment =',maxyear_sediment
        write(io_stdo_bgc,*) 'HAMOCC: maxyear_ocean =',maxyear_ocean
        write(io_stdo_bgc,*) 'HAMOCC: nburst_last =',nburst_last
        write(io_stdo_bgc,*) 'HAMOCC: lsed_rclim =',lsed_rclim
        write(io_stdo_bgc,*) 'HAMOCC: lsed_wclim =',lsed_wclim
        write(io_stdo_bgc,*) 'HAMOCC: lsed_spinup =',lsed_spinup
#endif
      ENDIF
#ifdef SED_OFFLINE
      if (lsed_rclim) lread_clim = .true.
      nburst = nburst_last ! set to 0 in namelist for a startup simulation
#endif
      atm_o2  = 196800.
      atm_n2  = 802000.
#ifdef natDIC
      atm_co2_nat = 284.32 ! CMIP6 pre-industrial reference
#endif   
#ifdef __c_isotopes
      atm_c13 = atm_co2
      atm_c14 = atm_co2

      ! Calculation of calibration factors
      PDB            = 0.0112372            ! The Pee Dee Belemnite reference 13C/12C ratio Keeling (1981)
      ref14c  	     = 1.176e-12            ! 14C/12C reference (pre-industrial) value Keeling (1981)
      prei_d13C_atm  = -6.5                 ! Preindustrial atmospheric d13C value in promille
      prei_dd14C_atm = 0.                   ! Preindustrial atmospheric dd14C value in promille
      atm_c13_cal    = (prei_d13C_atm/1000. + 1.) * PDB * atm_co2 / &   ! calibrated 13C model value for atmosphere [ppm]
        (1. + (prei_d13C_atm/1000. + 1.) * PDB)
      atm_dc14_cal   = 2.*(prei_d13C_atm+25.)/ &                        ! calibrated d14C model value for atmosphere [ppm]
        (1.-(2.*(prei_d13C_atm+25.)/1000.)) 
      atm_c14_cal    = (atm_dc14_cal/1000. + 1.) * ref14c * atm_co2     ! calibrated dd14C model value for atmosphere [ppm]
      factor_13c     = atm_c13_cal / atm_c13                            ! calibration factor 13C [-]
      factor_14c     = atm_c14_cal / atm_c14                            ! calibration factor 14C [-]
#endif
#endif   
!
! Biology
!
!ik note that the unit is kmol/m3!
      phytomi=1.e-11		!i.e. 1e-5 mmol P/m3 minimum concentration of phyto plankton (?js)
      grami=1.e-10              !i.e. 1e-5 mmol P/m3 minimum concentration of zooplankton

!ik addded parameter definition; taken from OCPROD.F
      remido=0.004*dtb      !1/d -remineralization rate (of DOM)
      dyphy=0.004*dtb       !1/d -mortality rate of phytoplankton 
      grazra=1.2*dtb        !1/d -grazing rate
      spemor=3.*1.e6*dtb    !1/d -mortality rate
      gammap=0.04*dtb       !1/d -exudation rate
      gammaz=0.06*dtb       !1/d -excretion rate
      ecan=0.95             ! fraction of mortality as PO_4
      pi_alpha=0.02*0.4     ! initial slope of production vs irradiance curve (alpha) (0.002 for 10 steps per day)
#ifdef AGG
      zinges = 0.5          !dimensionless fraction -assimilation efficiency
      epsher = 0.9          !dimensionless fraction -fraction of grazing egested
#elif defined(WLIN)
      zinges = 0.7          !dimensionless fraction -assimilation efficiency
      epsher = 0.9          !dimensionless fraction -fraction of grazing egested
#else
      zinges = 0.6          !dimensionless fraction -assimilation efficiency
      epsher = 0.8          !dimensionless fraction -fraction of grazing egest      
#endif

#ifdef __c_isotopes
! for avoiding too many tracers, surface gross rates work with reduced
! values bifr13 and bifr14
      bifr13=0.98                ! biogenic fractionation used in ocprod.F90
      bifr14=0.98

! decay parameter for sco214, HalfLive = 5730 years
      c14_t_half 	= 5730.*365.                 ! Half life of 14C [days]	
      c14dec 		= (log(2.)/c14_t_half)*dtb   ! The decay constant; labda [1/day]; c14dec[-]
#endif

! half sat. constants, note that the units are kmol/m3 !
      bkphy  = 4.e-8    !i.e. 0.04 mmol P/m3
      bkzoo  = 8.e-8    !i.e. 0.08 mmol P/m3
      bkopal = 5.e-6    !i.e. 4.0 mmol Si/m3

!sinking speed
      wpoc  =  5.*dtb       !m/d  iris : 5.
      wcal  = 30.*dtb       !m/d 
      wopal = 30.*dtb       !m/d  iris : 60
#ifdef WLIN
      wmin  =  1.*dtb       !m/d   minimum sinking speed
      wmax  = 60.*dtb       !m/d   maximum sinking speed
      wlin  = 60./2400.*dtb !m/d/m constant describing incr. with depth, r/a=1.0
#endif

      
! deep see remineralisation constants
      drempoc  = 0.025*dtb    !1/d
      dremdoc  = 0.004*dtb    !1/d
      dphymor  = 0.004*dtb    !1/d
      dzoomor  = 3.*1.e6*dtb  !1/d -mortality rate
      dremopal = 0.003*dtb    !1/d
      dremn2o  = 0.01*dtb     !1/d
      dremsul  = 0.005*dtb    ! remineralization rate for sulphate reduction 
      

! nirogen fixation by blue green algae
      bluefix=0.005*dtb       !1/d

! N2-Fixation following the parameterization in Kriest and Oschlies, 2015.
! Factors tf2, tf1 and tf0 are a polynomial (2nd order) 
! approximation to the functional relationship by Breitbarth et al. (2007),
! for temperature dependence of Trichodesmium growth, their eq. (2)
! The relation will be scaled to their max. growth rate, tff.
! Note that the second order approx. is basically similar to their
! function 2 for T-dependent nitrogen fixation multiplied by 4 
! (2 [N atoms per mole] * 12 [light hrs per day]/6 [C-atoms per N-atoms])
      tf2 = -0.0042
      tf1 = 0.2253
      tf0 = -2.7819
      tff = 0.2395  

! extended redfield ratio declaration
! Note: stoichiometric ratios are based on Takahashi etal. (1985)
! P:N:C:-O2 + 1:16:122:172
      ro2ut=172. 
      rcar=122.
      rnit=16.
      rnoi=1./rnit

! stoichiometric ratios for denitrification from Paulmier et al. 2009, Table 1 and
! equation 18. Note that their R_0=ro2ut-2*rnit.
      rdnit0=0.8*ro2ut           ! moles nitrate lost for remineralisation of 1 mole P
      rdnit1=0.8*ro2ut-rnit      ! moles nitrate net  for remineralisation of 1 mole P
      rdnit2=0.4*ro2ut           ! moles N2 released  for remineralisation of 1 mole P

! stoichiometric ratios for N2O loss by "intermediate dinitrification". Note that there
! is no nitrate created by this process, organic N is released as N2
      rdn2o1=2*ro2ut-2.5*rnit    ! moles N2O used for remineralisation of 1 mole P
      rdn2o2=2*ro2ut-2*rnit      ! moles N2 released  for remineralisation of 1 mole P


#ifdef AGG
      rcalc = 14.  ! calcium carbonate to organic phosphorous production ratio
      ropal = 10.5 ! opal to organic phosphorous production ratio      
      calmax= 0.20
#elif defined(WLIN)
      rcalc = 38.  ! calcium carbonate to organic phosphorous production ratio
      ropal = 45.  ! opal to organic phosphorous production ratio      
#else
      rcalc = 40.  ! iris 40 !calcium carbonate to organic phosphorous production ratio
      ropal = 30.  ! iris 25 !opal to organic phosphorous production ratio      
#endif

! parameters for sw-radiation attenuation
! Analog to Moore et al., Deep-Sea Research II 49 (2002), 403-462
! 1 kmolP = (122*12/60)*10^6 mg[Chlorophyl] 

      ctochl  = 60.        ! C to Chlorophyl ratio
      atten_w = 0.04       ! yellow substances attenuation in 1/m
      atten_c = 0.03*rcar*(12./ctochl)*1.e6  ! phytoplankton attenuation in 1/m 
      atten_f = 0.4        ! fraction of sw-radiation directly absorbed in surface layer 
                           ! (only if FB_BGC_OCE) [feedback bgc-ocean]
      	
!ik for interaction with sediment module
      o2ut=172.
      rno3=16.

!ik weight percent iron in dust deposition times Fe solubility
! the latter three values come from Johnson et al., 1997
      fetune=0.6                  ! factor introduced to tune deposistion/solubility
      perc_diron = fetune * 0.035 * 0.01 / 55.85
      riron= 5.*rcar*1.e-6        ! fe to P ratio in organic matter
      fesoly=0.5*1.e-9            ! max. diss. iron concentration in deep water 
      relaxfe = 0.05/365.*dtb

!                        
! Set constants for calculation of dms ( mo_carbch )
! Parameters are a result from kettle optimisation 02.03.04

       dmspar(6)=0.100000000E-07  !0 half saturation microbial
       dmspar(5)=1.25*0.109784522E-01  !2*0.02   production with delsil
       dmspar(4)=1.25*0.107638502E+00  !2*1.3e-5 production with delcar
       dmspar(3)=0.0864 ! following Kloster et al., 06 Table 1 with 50% reduction to reduce bacterial removal and increase dms emissions
       dmspar(2)=0.0011 ! following Kloster et al., 06 Table 1
       dmspar(1)=10.              !2*5.     production with temp


      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      WRITE(io_stdo_bgc,*)                                             &
     &'* '
      WRITE(io_stdo_bgc,*)                                             &
     &'* Values of BELEG_BGC variables : '
#ifndef DIFFAT      
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_co2      = ',atm_co2      
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_o2       = ',atm_o2           
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atm_n2       = ',atm_n2 
#endif  
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              phytomi      = ',phytomi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grami        = ',grami
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              remido       = ',remido
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dyphy        = ',dyphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              zinges       = ',zinges
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              epsher       = ',epsher
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              grazra       = ',grazra
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              spemor       = ',spemor
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammap       = ',gammap
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              gammaz       = ',gammaz
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ecan         = ',ecan    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkphy        = ',bkphy
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkzoo        = ',bkzoo    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bkopal       = ',bkopal    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wpoc         = ',wpoc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wcal         = ',wcal    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wopal        = ',wopal   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              drempoc      = ',drempoc    
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremdoc      = ',dremdoc   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dremopal     = ',dremopal   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dphymor      = ',dphymor   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dzoomor      = ',dzoomor   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              bluefix      = ',bluefix   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ro2ut        = ',ro2ut   
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcar         = ',rcar 
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnit         = ',rnit
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rnoi         = ',rnoi
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rdnit1       = ',rdnit1
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rdnit2       = ',rdnit2
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rcalc        = ',rcalc
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ropal        = ',ropal
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              ctochl       = ',ctochl
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_w      = ',atten_w
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_c      = ',atten_c
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              atten_f      = ',atten_f
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              o2ut         = ',o2ut
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              rno3         = ',rno3
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              perc_diron   = ',perc_diron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              riron        = ',riron
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              fesoly       = ',fesoly
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              relaxfe      = ',relaxfe
#ifdef __c_isotopes
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              c14dec       = ',c14dec
#endif
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(1)    = ',dmspar(1)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(2)    = ',dmspar(2)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(3)    = ',dmspar(3)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(4)    = ',dmspar(4)
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dmspar(5)    = ',dmspar(5)	    
      ENDIF

#ifndef AGG
      dustd1 = 0.0001 !cm = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustsink = (9.81 * 86400. / 18.                  &  ! g * sec per day / 18.
     &         * (claydens - 1025.) / 1.567 * 1000.    &  !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
      wdust = dustsink

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd1       = ',dustd1	
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustd2       = ',dustd2	 
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              dustsink     = ',dustsink	
      WRITE(io_stdo_bgc,*)                                             &
     &'*                              wdust        = ',wdust	              
      ENDIF
#endif

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      ENDIF

#ifdef AGG
! parameters needed for the aggregation module

      SinkExp = 0.62
      FractDim = 1.62
!      Stick = 0.40
!      Stick = 0.25
      Stick = 0.15
      cellmass = 0.012 / rnit ![nmol P]
!ik      cellmass = 0.0039/ rnit ![nmol P] for a 10 um diameter
      cellsink = 1.40 *dtb! [m/d]
!ik      cellsink = 0.911 *dtb! [m/d]  for a 10 um diameter
!      shear = 86400. !shear in the mixed layer,   1.0  d-1
!      shear = 64800. !shear in the mixed layer,   0.75 d-1
      shear = 43200. !shear in the mixed layer,   0.5  d-1
      fsh = 0.163 * shear *dtb
      fse = 0.125 * 3.1415927 * cellsink * 100.
      alow1 = 0.002 !diameter of smallest particle [cm]
!ik      alow1 = 0.001 !diameter of smallest particle [cm]
      alow2 = alow1 * alow1
      alow3 = alow2 * alow1
!      alar1 = 1.0 !diameter of largest particle for size dependend aggregation and sinking [cm]
!      alar1 = 0.75 !diameter of largest particle for size dependend aggregation and sinking [cm]
      alar1 = 0.5 !diameter of largest particle for size dependend aggregation and sinking [cm]
      vsmall = 1.e-9 
      safe = 1.e-6     
      pupper = safe/((FractDim+safe)*cellmass)
      plower = 1./(1.1*cellmass)
      zdis = 0.01 / ((FractDim + 0.01)*cellmass)
      nmldmin = 0.1 ! minimum particle number in mixed layer

! Determine maximum sinking speed
      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*)
      write(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      write(io_stdo_bgc,*) 'HAMOCC aggregate sinking scheme:'
      write(io_stdo_bgc,*) ' Maximum sinking speed for aggregates of '
      write(io_stdo_bgc,*) ' maximum size ', alar1, ' cm is '
      write(io_stdo_bgc,*)   cellsink/dtb*(alar1/alow1)**SinkExp, ' m/day'
      endif

      alar2 = alar1 * alar1
      alar3 = alar2 * alar1
      TSFac = (alar1/alow1)**SinkExp
      TMFac = (alar1/alow1)**FractDim

! for shear aggregation of dust:
      dustd1 = 0.0001 !cm = 1 um, boundary between clay and silt
      dustd2=dustd1*dustd1
      dustd3=dustd2*dustd1
      dustsink = (9.81 * 86400. / 18.                & ! g * sec per day / 18.                 
     &         * (claydens - 1025.) / 1.567 * 1000.  & !excess density / dyn. visc.
     &         * dustd2 * 1.e-4)*dtb
      IF (mnproc.eq.1) THEN
      write(io_stdo_bgc,*) ' dust diameter (cm)', dustd1
      write(io_stdo_bgc,*) ' dust sinking speed (m/d)', dustsink / dtb
      if(dustsink.gt.cellsink) then 
        write(io_stdo_bgc,*) ' dust sinking speed greater than cellsink'
        dustsink=cellsink
        write(io_stdo_bgc,*) ' set dust sinking speed to cellsink'
      endif
      write(io_stdo_bgc,*)                                             &
     &'****************************************************************'
      endif

#endif /*AGG*/  


      DO  i=1,kpie
      DO  j=1,kpje
         keqb(:,i,j)=rmasko
      ENDDO
      ENDDO

      DO  j=1,kpje
      DO  i=1,kpie
      DO  k=1,ks
      DO  l=1,nsedtra
         sedlay(i,j,k,l)=0.
	 burial(i,j,l)=0.   ! last and final sediment layer
      ENDDO
      ENDDO
      ENDDO
      ENDDO

      DO  j=1,kpje
      DO  i=1,kpie
         expoor(i,j)=0.
         expoca(i,j)=0.
         exposi(i,j)=0.
      ENDDO
      ENDDO


! Initialise ocean tracers with WOA and GLODAP data
      call profile_gd(kpie,kpje,kpke,pglon,pglat,ptiestw,omask,TRIM(path))

      do k=1,kpke
      do j=1,kpje
      do i=1,kpie
        if (omask(i,j) .gt. 0.5 ) then
          ! convert WOA tracers kmol/m^3 -> mol/kg; GLODAP dic and alk
          ! are already in mol/kg. We need these units here, since after 
          ! initialisation the tracer field is passed to the ocean model
          ! first where units are mol/kg.
          ocetra(i,j,k,iphosph) = ocetra(i,j,k,iphosph)/prho(i,j,k)
          ocetra(i,j,k,ioxygen) = ocetra(i,j,k,ioxygen)/prho(i,j,k)
          ocetra(i,j,k,iano3)   = ocetra(i,j,k,iano3)  /prho(i,j,k)
          ocetra(i,j,k,isilica) = ocetra(i,j,k,isilica)/prho(i,j,k)
        endif
      enddo
      enddo
      enddo


! Initialise remaining ocean tracers 
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j) .GT. 0.5) THEN
         ocetra(i,j,k,igasnit)=1.e-10
         ocetra(i,j,k,idoc)   =1.e-8
         ocetra(i,j,k,iphy)   =1.e-8 
         ocetra(i,j,k,izoo)   =1.e-8 
         ocetra(i,j,k,idet)   =1.e-8 
         ocetra(i,j,k,icalc)  =0. 
         ocetra(i,j,k,iopal)  =1.e-8 
!#ifdef __c_isotopes
!	 ocetra(i,j,k,isco214)=0.75*2.27e-3 !Paddy: oldest water: 1600y --> X0.83 
!#endif
         ocetra(i,j,k,ian2o)  =0. 
         ocetra(i,j,k,idms)   =0. 
         ocetra(i,j,k,ifdust) =0. 
         ocetra(i,j,k,iiron)  =fesoly
         ocetra(i,j,k,iprefo2)=0.
         ocetra(i,j,k,iprefpo4)=0.
         ocetra(i,j,k,iprefalk)=0.
         hi(i,j,k)            =1.e-8
         co3(i,j,k)           =0.
#ifdef __c_isotopes
         ocetra(i,j,k,isco213)=ocetra(i,j,k,isco212)     ! 12C because 12-C-deviation is wanted
         ocetra(i,j,k,isco214)=ocetra(i,j,k,isco212)
!         ocetra(i,j,k,isco213)=2.27e-3     ! adjusted to reference ratio 13C/12C=1 (*100)!
!         ocetra(i,j,k,isco214)=2.27e-3
         ocetra(i,j,k,idet13) =1.e-8
         ocetra(i,j,k,idet14) =1.e-8
         ocetra(i,j,k,icalc13)=0.
         ocetra(i,j,k,icalc14)=0.
#endif

#ifdef AGG
! calculate initial numbers from mass, to start with appropriate size distribution
         snow = (ocetra(i,j,k,iphy)+ocetra(i,j,k,idet))*1.e+6
         ocetra(i,j,k,inos)   = snow / cellmass / (FractDim+1.)
         ocetra(i,j,k,iadust) =0. 
#endif /*AGG*/

#ifdef ANTC14
         ocetra(i,j,k,iantc14)=ocetra(i,j,k,isco214)
#endif
#ifdef CFC
         ocetra(i,j,k,icfc11)=0.
         ocetra(i,j,k,icfc12)=0.
         ocetra(i,j,k,isf6)=0.
#endif
#ifdef natDIC
         nathi(i,j,k)            =1.e-8
         natco3(i,j,k)           =0.
         ocetra(i,j,k,inatcalc)  =0. 
#endif
      ENDIF ! omask > 0.5

      ENDDO
      ENDDO
      ENDDO

! Initialise preformed tracers in the mixed layer; note that the 
! whole field has been initialised to zero above
      DO k=1,kmle
      DO j=1,kpje
      DO i=1,kpie

      IF(omask(i,j) .GT. 0.5) THEN
         ocetra(i,j,k,iprefo2) =ocetra(i,j,k,ioxygen)
         ocetra(i,j,k,iprefpo4)=ocetra(i,j,k,iphosph)
         ocetra(i,j,k,iprefalk)=ocetra(i,j,k,ialkali)
      ENDIF

      ENDDO
      ENDDO
      ENDDO

! read in dust fields
     CALL GET_DUST(kpie,kpje,kpke,omask,path)
!
!read in pi_ph fields
!     CALL GET_PI_PH(kpie,kpje,kpke,omask,path)
!  Initial values for sediment pore water tracer.
      DO  k=1,ks
      DO  j=1,kpje
      DO  i=1,kpie 
      IF(omask(i,j) .GT. 0.5) THEN
         powtra(i,j,k,ipowaic)=ocetra(i,j,kbo(i,j),isco212)
         powtra(i,j,k,ipowaal)=ocetra(i,j,kbo(i,j),ialkali)
         powtra(i,j,k,ipowaph)=ocetra(i,j,kbo(i,j),iphosph)
         powtra(i,j,k,ipowaox)=ocetra(i,j,kbo(i,j),ioxygen)
         powtra(i,j,k,ipown2) =0.
         powtra(i,j,k,ipowno3)=ocetra(i,j,kbo(i,j),iano3)
         powtra(i,j,k,ipowasi)=ocetra(i,j,kbo(i,j),isilica)      
         sedlay(i,j,k,issso12)=1.e-8
         sedlay(i,j,k,isssc12)=1.e-8
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=1.e-8
         sedlay(i,j,k,isssc13)=1.e-8
         sedlay(i,j,k,issso14)=1.e-8
         sedlay(i,j,k,isssc14)=1.e-8
	 powtra(i,j,k,ipowc13)=1.e-8
	 powtra(i,j,k,ipowc14)=1.e-8
#endif
         sedlay(i,j,k,issster)=30.
         sedlay(i,j,k,issssil)=1.e-8
         sedhpl(i,j,k)        =hi(i,j,kbo(i,j))
      ELSE
         powtra(i,j,k,ipowno3)=rmasks
         powtra(i,j,k,ipown2) =rmasks
         powtra(i,j,k,ipowaic)=rmasks
         powtra(i,j,k,ipowaal)=rmasks
         powtra(i,j,k,ipowaph)=rmasks
         powtra(i,j,k,ipowaox)=rmasks
         powtra(i,j,k,ipowasi)=rmasks
         sedlay(i,j,k,issso12)=rmasks
         sedlay(i,j,k,isssc12)=rmasks
#ifdef __c_isotopes
         sedlay(i,j,k,issso13)=rmasks
         sedlay(i,j,k,isssc13)=rmasks
         sedlay(i,j,k,issso14)=rmasks
         sedlay(i,j,k,isssc14)=rmasks
	 powtra(i,j,k,ipowc13)=rmasks
	 powtra(i,j,k,ipowc14)=rmasks
#endif
         sedlay(i,j,k,issssil)=rmasks
         sedlay(i,j,k,issster)=rmasks
         sedlay(i,j,k,issssil)=rmasks
         sedhpl(i,j,k)        =rmasks
      ENDIF
      ENDDO
      ENDDO
      ENDDO

!
! values for sediment fluxes
!
      
      DO j=1,kpje
      DO i=1,kpie
      DO k=1,npowtra
        sedfluxo(i,j,k)=0.     
      ENDDO
      ENDDO
      ENDDO
!
! values for sedimentation
!
      
      DO j=1,kpje
      DO i=1,kpie
        prorca(i,j)=0.
        prcaca(i,j)=0.
        silpro(i,j)=0.
        produs(i,j)=0.
      ENDDO
      ENDDO
      
!
!  atmospheric diffusion parameters.
!    
#ifdef DIFFAT
      DO  j=1,kpje
      DO  i=1,kpie 
         atm(i,j,iatmco2) = 278.
         atm(i,j,iatmo2)  = 196800.
         atm(i,j,iatmn2)  = 802000.
#ifdef __c_isotopes
         atm(i,j,iatmc13) = 278.
         atm(i,j,iatmc14) = 278.
#endif
         atdifv(i,j)=1.
      ENDDO
      ENDDO
#endif
   
#ifdef ANTC14
      DO  i=1,kpie
        ii=1+(i+p_ioff-1)   ! global i-index for giph_g
        north=1.
        DO  j=1,kpje
          jj=1+(j+p_joff-1) ! global j-index for giph_g
          north=pglat(ii,jj)
          if (north .gt. 20.) then
            Rbomb(i,j) = D14C_north
          endif
          if (north .lt. -20.) then
            Rbomb(i,j) = D14C_south
          endif
          if ((north .le. 20.).and.(north .ge. -20.)) then
            Rbomb(i,j) = D14C_equator
          endif
!         WRITE(io_stdo_bgc,*)'Rbomb: ',i,j,north,Rbomb(i,j)
        ENDDO
      ENDDO
#endif


#ifdef FB_BGC_OCE
      DO k=1,kpke
      DO j=1,kpje
      DO i=1,kpie
        abs_oce(i,j,k)=1.
      ENDDO
      ENDDO
      ENDDO
#endif

      RETURN
      END
