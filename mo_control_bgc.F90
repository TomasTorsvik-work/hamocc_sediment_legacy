      MODULE mo_control_bgc

!$Source: /scratch/local1/m212047/patrick/SRC_MPI/src_hamocc/RCS/mo_control_bgc.f90,v $\\
!$Revision: 1.1 $\\
!$Date: 2005/01/28 08:37:45 $\\

!***********************************************************************
!
!**** *MODULE mo_control_bgc* - control variables for bgc modules.
!
!     S.Legutke         *MPI-MaD, HH*  2002-02-28
!
!     Modified
!     --------
!     Marco van Hulten  *GFI, Bergen*  2018-04-18
!     - added variables specific for sediment spin-up
!     
!     Purpose
!     -------
!     - declaration
!
!
!**********************************************************************
      implicit none

! Logical unit number for I/O.

      INTEGER :: io_stdo_bgc           !  standard out.
      INTEGER :: io_stdi_bgc           !  standard in.

      INTEGER :: io_rsti_bgc           !  restart in. 
      INTEGER :: io_rsto_bgc           !  restart out. 

      INTEGER :: io_nml                !  namelist

! Control variables

      REAL    :: dtbgc            !  HAMOCC time step length [sec].
      REAL    :: dtb              !  HAMOCC time step length [days].
      INTEGER :: ndtdaybgc        !  time steps per day.
      REAL    :: dtoff            !  off-line sediment time step length [sec].
      REAL    :: dto              !  off-line sediment time step length [days].
      REAL    :: dtsed            !  sediment time step length [sec].
      REAL    :: dts              !  sediment time step length [days].
      REAL    :: rdtsed           !  ratio of sediment over biogeochemistry time step.

      INTEGER :: ldtbgc           !  time step number from bgc restart file
      INTEGER :: ldtrunbgc        !  actual time steps of run.
#if defined(SED_OFFLINE)
      INTEGER :: nstep_in_month   !  accumulation counter for SED_OFFLINE.
      INTEGER :: maxyear_sediment !  number of years for off-line sediment integration.
      INTEGER :: maxyear_ocean    !  number of years for full MICOM-HAMOCC integration.
      INTEGER :: nburst_last      !  nburst from the end of the previous simulation (startup: 0).
      INTEGER :: nburst           !  counter of running sediment off-line.
      LOGICAL :: lsed_rclim       !  whether to read bottom seawater climatology from file (nml).
      LOGICAL :: lsed_wclim       !  whether to write bottom seawater climatology to file (nml).
      LOGICAL :: lsed_spinup      !  whether to spin up the sediment (nml).
      LOGICAL :: lread_clim       !  whether reading the climatology now.
      LOGICAL :: lwrite_clim      !  whether writing the climatology now.
      LOGICAL :: lcompleted_clim  !  whether we have a recent climatology available.
#endif
      LOGICAL :: lspinning_up_sed !  whether spinning up the sediment now.
      INTEGER :: nyear_global     !  ocean model year number, including sediment().


      INTEGER :: icyclibgc        !  switch for cyclicity.
      INTEGER :: ndtrunbgc        !  total no. of time steps of run.


      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL    :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL    :: rmasko = 99999.00  !  value at wet cells in ocean.


      INTEGER :: kchck = 0          !  switch for extended print control (0/1). 
      
      END MODULE mo_control_bgc
