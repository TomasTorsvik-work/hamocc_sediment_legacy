      SUBROUTINE bodensed(kpie,kpje,kpke,pddpo)

!
!$Source: /server/cvs/mpiom1/mpi-om/src_hamocc/bodensed.f90,v $\\
!$Revision: 1.2 $\\
!$Date: 2004/11/12 15:37:21 $\\
!$Name:  $\\
!
!**********************************************************************
!
!**** *BODENSED* - .
!
!     Ernst Maier-Reimer,    *MPI-Met, HH*    10.04.01
!
!     Modified
!     --------
!     S.Legutke,        *MPI-MaD, HH*    10.04.01
!
!     Purpose
!     -------
!     set up of sediment layer.
!
!     Method:
!     ------
!     -
!
!     *CALL*       *BODENSED
!
!     *PARAMETER*  *PARAM1.h*     - grid size parameters for ocean model.
!     *COMMON*     *PARAM1_BGC.h* - declaration of ocean/sediment tracer.
!     *COMMON*     *COMMO1_BGC.h* - ocean/sediment tracer arrays.
!     *COMMON*     *UNITS_BGC.h*  - std I/O logical units.
!
!**   Interface to ocean model (parameter list):
!     -----------------------------------------
!
!     *INTEGER* *kpie*    - 1st dimension of model grid.
!     *INTEGER* *kpje*    - 2nd dimension of model grid.
!     *INTEGER* *kpke*    - 3rd (vertical) dimension of model grid.
!     *REAL*    *pddpo*   - size of scalar grid cell (3rd dimension) [m].
!     *REAL*    *dtsed*   - sediment model time step [sec].
!
!     Externals
!     ---------
!     none.
!**********************************************************************
! evaluates the min depth of all layers, to be checked vs. sinking rate
!   in routine BELEG_BGC, to assure that flux rate is < 1 per time step

      USE mo_carbch
      USE mo_sedmnt
      USE mo_biomod

      USE mo_control_bgc
      USE mo_param1_bgc 
      USE mod_xc

      implicit none

      integer, intent(in) :: kpie,kpje,kpke
      real, intent(in) :: pddpo(kpie,kpje,kpke)
      integer :: i,j,k
      real :: sumsed

      dzs(1) = 0.001
      dzs(2) = 0.003
      dzs(3) = 0.005
      dzs(4) = 0.007
      dzs(5) = 0.009
      dzs(6) = 0.011
      dzs(7) = 0.013
      dzs(8) = 0.015
      dzs(9) = 0.017
      dzs(10) = 0.019
      dzs(11) = 0.021
      dzs(12) = 0.023
      dzs(13) = 0.025

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)  ' '
      WRITE(io_stdo_bgc,*)  'Sediment layer thickness [m] : '
      WRITE(io_stdo_bgc,'(5F9.3)') dzs
      WRITE(io_stdo_bgc,*)  ' '
      ENDIF
    
      porwat(1) = 0.85
      porwat(2) = 0.83
      porwat(3) = 0.8
      porwat(4) = 0.79
      porwat(5) = 0.77
      porwat(6) = 0.75
      porwat(7) = 0.73
      porwat(8) = 0.7
      porwat(9) = 0.68
      porwat(10) = 0.66
      porwat(11) = 0.64
      porwat(12) = 0.62

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)  'Pore water in sediment: ',porwat
      ENDIF

      
      seddzi(1)=500.
      DO 1131 k=1,ks
         porsol(k)=1.-porwat(k)
         IF(k.GE.2) porwah(k)=0.5*(porwat(k)+porwat(k-1))
         IF(k.EQ.1) porwah(k)=0.5*(1.+porwat(1))
         seddzi(k+1)=1./dzs(k+1)
         seddw(k)=0.5*(dzs(k)+dzs(k+1))
1131  CONTINUE

! ******************************************************************
! the following section is to include the SEDIMENT ACCELERATION
! mechanism to accelerate the sediment:

      sedac = 1./real(isac)

! determine total solid sediment thickness sumsed
! and reduced sediment volume
      sumsed = 0.
      do k = 1, ks
        porwat(k) = porwat(k) * sedac
        porsol(k) = porsol(k) * sedac
        sumsed = sumsed + seddw(k)
      enddo

! determine reduced diffusion rate sedict,
! and scaling for bottom layer ventilation, sedifti

      sedict=1.e-9*dtsed
      sedifti = sedict / (sumsed**2)
      sedict=sedict*sedac

      IF (mnproc.eq.1) THEN
      WRITE(io_stdo_bgc,*)  'sediment acceleration factor: ',isac
      WRITE(io_stdo_bgc,*)  'sediment area reduction: ',sedac
      WRITE(io_stdo_bgc,*)  'new diffusion is: ',sedict
      WRITE(io_stdo_bgc,*)  'total sediment thickness: ',sumsed
      WRITE(io_stdo_bgc,*)  'sedict / sumsed**2: ',sedifti
      WRITE(io_stdo_bgc,*)  'new pore water fraction: ',porwat
      ENDIF


! ******************************************************************

! ******************************************************************
! densities etc. for SEDIMENT SHIFTING

! define weight of calcium carbonate, opal, and poc [kg/kmol]
      calcwei = 100. ! 40+12+3*16 kg/kmol C
      opalwei = 60.  ! 28 + 2*16  kg/kmol Si
      orgwei  = 30.  ! from 12 kg/kmol * 2.5 POC[kg]/DW[kg]
                     ! after Alldredge, 1998: POC(g)/DW(g) = 0.4 of diatom marine snow, size 1mm3

! define densities of opal, caco3, poc [kg/m3]
      calcdens = 2600. 
      opaldens = 2200. 
      orgdens  = 1000. 
      claydens = 2600. !quartz
      
! define volumes occupied by solid constituents [m3/kmol]
      calfa = calcwei/calcdens
      oplfa = opalwei/opaldens
      orgfa = orgwei/orgdens
      clafa = 1./claydens !clay is calculated in kg/m3

! determine total solid sediment volume
      solfu = 0.
      DO k=1,ks
         solfu = solfu + seddw(k)*porsol(k)
      ENDDO

! ******************************************************************
! find lowest mass containing layer and its thickness

      if (.not. lspinning_up_sed) CALL calc_bot(kpie,kpje,kpke,pddpo)

      RETURN
      END
