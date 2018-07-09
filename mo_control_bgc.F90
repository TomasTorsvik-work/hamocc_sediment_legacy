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
      REAL    :: dtsed            !  sediment time step length [sec].
      REAL    :: dts              !  sediment time step length [days].
      INTEGER :: ndtdaybgc        !  time steps per day.
      REAL    :: ndtdaysed        !  time steps per day for sediment.

      INTEGER :: ldtbgc           !  time step number from bgc restart file
      INTEGER :: ldtrunbgc        !  actual time steps of run.
!#if defined(SEDI_OFFLINE)
      LOGICAL :: lrunsed = .false.!  whether to run the sediment spin-up code.
      INTEGER :: nstep_in_month   !  accumulation counter for SED_WCLIM.
      INTEGER :: imonth, iyear    !  counters that must be available to ncwrt_bgc().
      INTEGER :: maxyear_sediment !  number of years for off-line sediment integration.
      INTEGER :: maxyear_ocean    !  number of years for full MICOM-HAMOCC integration.
!#endif
      INTEGER :: nyear_global     !  ocean model year number, including sediment().
      LOGICAL :: must_read_clim   !  whether to read bottom seawater climatology from file.
      LOGICAL :: is_end_of_day    !  whether we are at the last timestep of the day.
      LOGICAL :: lsed_rclim, lsed_wclim, lsed_spinup


      INTEGER :: icyclibgc        !  switch for cyclicity.
      INTEGER :: ndtrunbgc        !  total no. of time steps of run.


      INTEGER :: isac             !  acceleration factor for sediment, read from namelist


      REAL    :: rmasks = 0.0       !  value at wet cells in sediment.
      REAL    :: rmasko = 99999.00  !  value at wet cells in ocean.


      INTEGER :: kchck = 0          !  switch for extended print control (0/1). 
      
      END MODULE mo_control_bgc
